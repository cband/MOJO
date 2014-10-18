
#include "DiscordantClusterFinder.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >, 
		(keywords::channel = "DiscordantClusterFinder"));

	int DiscordantClusterFinder::NextClusterID = 1000000;

	int DiscordantCluster::GetNumDiscordReads(bool uniqueOnly) 
	{
		int count = 0;
		for (auto al : alignments) {
			if (uniqueOnly && al.isUnique)
				count++;
			else if ( !uniqueOnly )
				count++;
		}
		return count;
	}

	vector<DiscordantCluster *> DiscordantClusterFinder::LoadDiscordantClusters()
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Loading discordant clusters...";

		vector<DiscordantCluster *> clusters_unfiltered, clusters;
		unordered_map<int, unordered_map<int, DiscordantCluster *> > clustersMap;
		ifstream discFlat(GetFlatFilename().c_str());
		for (std::string str; getline(discFlat, str);) {
			auto sp = Utils::SplitToVector(str, "\t");
			Gene *gA = gm->GetGeneFromMap(sp[0]);
			Gene *gB = gm->GetGeneFromMap(sp[5]);

			// abParts is the only black listed locus (~4,500 exons); including 
			// this works without any issues but increases the downstream 
			// processing by as much as 150% if discordant read pairs are found 
			// involving this gene
			if (gA->name.find("abParts") != string::npos ||
				gB->name.find("abParts") != string::npos)
				continue;

			if (gA == 0 || gB == 0) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error loading clusters: "
					<< " gene(s) not found for line: " << endl << str << endl;
				exit(1);
			}
			if (gA == gB || gA->chr == "chrM" || gB->chr == "chrM") 
				continue;
			DiscordantClusterAlignment dca(sp[1], sp[2], sp[3], sp[6], sp[7], 
				sp[8], sp[11]);


			//require at least one end not mapping to repetitive region
			if (gA->DoesRegionContainRepeat(dca.startA, dca.endA) &&
				gB->DoesRegionContainRepeat(dca.startB, dca.endB))
				continue;

			auto exonsA = gA->FindExonsWithinPoints(dca.startA, dca.endA);
			auto exonsB = gB->FindExonsWithinPoints(dca.startB, dca.endB);

			//require each end of the discordant reads to map to at least one exon
			//of each of the genes; 
			if (exonsA.size() == 0 || exonsB.size() == 0) 
				continue;

			//Uses a precomputed megablast index to determine if either end of the 
			//discordant read maps to a region of shared homology between the two genes
			bool blastHomologyFound = false;
			if (gA->blastHits.find(gB->geneId) != gB->blastHits.end()) {
				vector<CoordinatePair> coords = gA->blastHits[gB->geneId];
				for (auto cp : coords) 
					if ((cp.startA < dca.startA && cp.endA > dca.endA) ||
						(cp.startB < dca.startB && cp.endB > dca.endB))
						blastHomologyFound = true;
			}
			//If homology found, then ignore this discordant read
			if (blastHomologyFound)
				continue;

			if (gA->geneId > gB->geneId) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "Incorrect format for: "
					<< GetFlatFilename() << endl << str << endl;
				exit(1);
			}

			//Check if a cluster is already constructed that represents the fusion 
			//between the two genes
			DiscordantCluster *cluster = 0;
			if ((clustersMap[gA->geneId]).find(gB->geneId) == 
				(clustersMap[gB->geneId]).end()) 
			{
				cluster = new DiscordantCluster(gA, gB);
				clustersMap[gA->geneId][gB->geneId] = cluster;
				clusters_unfiltered.push_back(cluster);
			}
			else
				cluster = clustersMap[gA->geneId][gB->geneId];


			// Merge both vectors
			exonsA.insert(exonsA.end(), exonsB.begin(), exonsB.end());	
			//Count the number discordant reads mapping to each of the exons 
			//between the  two genes
			for (auto ex : exonsA) 
			{
				if (cluster->exonicReadCount.find(ex->exonId) ==
					cluster->exonicReadCount.end())
				{
					cluster->exonicReadCount[ex->exonId] = 0;
				}
				cluster->exonicReadCount[ex->exonId]++;
			}

			cluster->alignments.push_back(dca);
		}

		for (auto cluster : clusters_unfiltered) {
			
			if (cluster->alignments.size() < c->minSpanCount)
				continue;

			//Check if there are at least 2 read-pairs that span a potential breakpoint; 
			//ie, check if on a 2-d plane, if the two reads cluster to map to a 
			//potential breakpoint region.  Not the proper way of doing this, but want 
			//to contain the sensitivity issues with doing more k-means clustering.
			bool twoPairFound = false;
			Gene *gA = cluster->gA, *gB = cluster->gB;

			for (int j = 0; j < cluster->alignments.size() && !twoPairFound; j++) {
				DiscordantClusterAlignment X = cluster->alignments[j];
				for (int k = j + 1; 
					k < cluster->alignments.size() && !twoPairFound; k++) 
				{
					DiscordantClusterAlignment Y = cluster->alignments[k];
					int distA, distB;
					distA = gA->GetExonicDistanceBetweenCoords(X.startA, Y.startA);
					distB = gB->GetExonicDistanceBetweenCoords(X.startB, Y.startB);
					if (min(distA, abs(X.startA - Y.startA)) < 500 &&
						min(distB, abs(X.startB - Y.startB)) < 500 &&
						X.startA != Y.startA && X.startB != Y.startB)
						twoPairFound = true;
				}
			}
			if (!twoPairFound)
				continue;

			//Counting the number of discordant reads with each end mapping 
			//to exonic regions only
			int numExonicDiscordantReads = 0;
			vector<DiscordantClusterAlignment> filteredAligns;
			for (auto dca : cluster->alignments) {
				auto eWA = gA->FindExonsWithinPoints(dca.startA, dca.endA, false);
				auto eWB = gB->FindExonsWithinPoints(dca.startB, dca.endB, false);
				if (numExonicDiscordantReads < 2 
					&& eWA.size() != 0 && eWB.size() != 0)
				{
					numExonicDiscordantReads++;
				}
			}
			if (numExonicDiscordantReads == 0)
				continue;

			//Require a minimum number of reads to map to at least one full exon
			int minCountPerExon = 0;
			for (auto exonReadCount : cluster->exonicReadCount) 
				if (minCountPerExon < exonReadCount.second)
					minCountPerExon = exonReadCount.second;

			if (minCountPerExon >= (c->minSpanCount-1))
				clusters.push_back(cluster);
		}
		std::sort(clusters.begin(), clusters.end(), 
			DiscordantCluster::compareClusters);

		//Output discordant clusters passing the above heuristic filters
		ofstream discGenes(
			(c->workingDir + "discordants.genes").c_str(), ios::out);
		for (auto cl : clusters){
			string exonsStr = "", countsStr = "";
			for (auto exonReadCount : cl->exonicReadCount) {
				exonsStr += lexical_cast<string>(exonReadCount.first) + ",";
				countsStr += lexical_cast<string>(exonReadCount.second) + ",";
			}
			discGenes << cl->clusterId << "\t" << cl->gA->geneId << "\t" 
				<< cl->gA->name << "\t" << cl->gB->geneId << "\t" << cl->gB->name
				<< "\t" << cl->GetNumDiscordReads() << "\t" << exonsStr << "\t" 
				<< countsStr << "\t" << cl->GetNumDiscordReads(true) << endl;
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\t" << clusters.size() << " clusters found";
		return clusters;
	}

	string DiscordantClusterFinder::GetFlatFilename()
	{
		Config *c = Config::GetConfig();
		return c->workingDir + "discordant.flat";
	}
}