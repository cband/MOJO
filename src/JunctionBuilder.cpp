
#include "JunctionBuilder.h"

namespace MOJO {

	int JunctionBuilder::BuildJunctions(vector<DiscordantCluster *> clusters)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		int neighbor_cnt = 3, junctionCount = 0;
		string junctsFileName = c->workingDir + "/junctions.fa";
		ofstream junctionsFileStream(junctsFileName.c_str(), ios::out);

		BamTools::BamReader reader;
		BamTools::BamAlignment aln;
		string alignmentBam = c->workingDir + "/splits.alignments.sorted.bam";
		reader.Open(alignmentBam, alignmentBam + ".bai", true, true);
		BamTools::RefVector refvect = reader.GetReferenceData();
		unordered_map<int, int> geneToRefID;
		for (int i = 0; i < refvect.size(); i++) {
			string nm = refvect[i].RefName;
			int geneId = (nm.substr(0, 1) == "g" ? 
				stoi(nm.substr(1, nm.length() - 1)) : stoi(nm));
			geneToRefID[geneId] = i;
		}
		int counter = 0;
		for (auto cluster : clusters) {
			Gene *gA = cluster->gA, *gB = cluster->gB;
			bool isReadThruAB = Gene::IsReadThru(gA, gB, c->readThruDist);
			bool isReadThruBA = Gene::IsReadThru(gB, gA, c->readThruDist);

			unordered_map<int, Exon*> candExonsA, candExonsB;
			for (auto iso : gA->isoforms) {
				for (int i = 0; i < iso->exons.size(); i++) {
					Exon *ex = iso->exons[i];
					if (cluster->exonicReadCount.find(ex->exonId) !=
						cluster->exonicReadCount.end()) 
					{
						for (int k = max(i - neighbor_cnt, 0);
							k < min((int)iso->exons.size(), i + neighbor_cnt + 1); k++)
							candExonsA[iso->exons[k]->exonId] = iso->exons[k];
					}
				}
			}
			for (auto iso : gB->isoforms) {
				for (int i = 0; i < iso->exons.size(); i++) {
					Exon *ex = iso->exons[i];
					if (cluster->exonicReadCount.find(ex->exonId) !=
						cluster->exonicReadCount.end()) 
					{
						for (int k = max(i - neighbor_cnt, 0);
							k < min((int)iso->exons.size(), i + neighbor_cnt + 1); k++)
							candExonsB[iso->exons[k]->exonId] = iso->exons[k];
					}
				}
			}

			vector<Exon *> exonsA, exonsB;
			for (auto candExonA : candExonsA)
				exonsA.push_back(candExonA.second);
			for (auto candExonB : candExonsB)
				exonsB.push_back(candExonB.second);

			//std::sort(exonsA.begin(), exonsA.end(), Exon::compareExons);
			//std::sort(exonsB.begin(), exonsB.end(), Exon::compareExons);

			//int minA = exonsA[0]->exStart, maxA = exonsA[exonsA.size() - 1]->exEnd;
			//int minB = exonsB[0]->exStart, maxB = exonsB[exonsB.size() - 1]->exEnd;

			string seqA = gA->GetTranscribedSequence();
			string seqB = gB->GetTranscribedSequence();

			unordered_map<string, pair<string, string> > fetchedReads;
			for (auto ex : exonsA) {
				int idA = geneToRefID[gA->geneId];
				if (!reader.SetRegion(idA, ex->exStart, idA, ex->exEnd))
					continue;
				while (reader.GetNextAlignment(aln)){
					string tagValue;
					if (!aln.GetTag("XJ", tagValue))
						continue;
					auto sp = Utils::SplitToVector(tagValue, ",");
					pair<string, string> p(sp[2], sp[3]);
					fetchedReads[sp[0] + "_" + sp[1]] = p;
				}
			}

			for (auto ex : exonsB) {
				int idB = geneToRefID[gB->geneId];
				if (!reader.SetRegion(idB, ex->exStart, idB, ex->exEnd))
					continue;
				while (reader.GetNextAlignment(aln)){
					string tagValue;
					if (!aln.GetTag("XJ", tagValue))
						continue;
					auto sp = Utils::SplitToVector(tagValue, ",");
					pair<string, string> p(sp[2], sp[3]);
					fetchedReads[sp[0] + "_" + sp[1]] = p;
				}
			}

			//for each fetched read, throw out if full length matches to gA or 
			char *gAseq = Utils::StringCopy(seqA.c_str());
			char *gBseq = Utils::StringCopy(seqB.c_str());

			for (auto fetch : fetchedReads)
			{
				counter++;
				string rd1 = fetch.second.first, rd2 = fetch.second.second;
				char *a = Utils::StringCopy(rd1.c_str());
				char *b = Utils::StringCopy(rd2.c_str());
				SWseqpair *prob, *res;
				prob = new SWseqpair(gAseq, a);
				res = SWalign::SmithWaterman(prob, false);
				delete res, prob;
				prob = new SWseqpair(gBseq, b);
				res = SWalign::SmithWaterman(prob, false);
				delete res, prob;
				prob = new SWseqpair(gAseq, a);
				res = SWalign::SmithWaterman(prob, false);
				delete res, prob;
				prob = new SWseqpair(gBseq, b);
				res = SWalign::SmithWaterman(prob, false);
				delete res, prob;
				
				/*			SWseqpair *res;
				res = SWalign::SmithWaterman(rd1, seqA, false);
				res = SWalign::SmithWaterman(rd2, seqA, false);
				res = SWalign::SmithWaterman(rd1, seqB, false);
				res = SWalign::SmithWaterman(rd2, seqB, false);*/
				if (++counter % 100 == 0)
					cout << Utils::GetTimeStamp() << "\t" << counter << endl;
			}

			//Build exon-exon junction database;
			for (auto candExonA : candExonsA){
				for (auto candExonB : candExonsB) {
					junctionCount++;
					Exon *exA = candExonA.second, *exB = candExonB.second;
					stringstream ss;
					if (!isReadThruAB) {
						vector<string> exA_5p_seq = exA->GetJunctionSeqAs5pExon();
						vector<string> exB_3p_seq = exB->GetJunctionSeqAs3pExon();
						int idx = 1;
						for (auto seq5p : exA_5p_seq) {
							for (auto seq3p : exB_3p_seq)
							{
								if (seq5p.length() < 20 || seq3p.length() < 20)
									continue;
								Junction *junct = new Junction(cluster->clusterId,
									cluster, exA, exB, "F", idx++, seq5p, seq3p);
								ss << junct->GetJunctionFa() << endl;
								cluster->junctions.push_back(junct);
							}
						}
					}
					if (!isReadThruBA) {
						vector<string> exB_5p_seq = exB->GetJunctionSeqAs5pExon();
						vector<string> exA_3p_seq = exA->GetJunctionSeqAs3pExon();
						int idx = 1;
						for (auto seq5p : exB_5p_seq)
							for (auto seq3p : exA_3p_seq)
							{
							if (seq5p.length() < 20 || seq3p.length() < 20)
								continue;
							Junction *junct = new Junction(cluster->clusterId,
								cluster, exB, exA, "R", idx++, seq5p, seq3p);
							ss << junct->GetJunctionFa() << endl;
							cluster->junctions.push_back(junct);
							}
					}
					junctionsFileStream << ss.str();
				}
			}
		}
		return junctionCount;
	}


}