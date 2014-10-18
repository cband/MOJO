
#include "FusionCompiler.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = ""));

	// FusionCompiler::Run is the main function in this singleton class.  It 
	// processes the discordant clusters with the the junctions and filtered 
	// anchor reads loaded.  Outputs the .fusions and .fusions.pileup file
	void FusionCompiler::Run(vector <DiscordantCluster *> clusters)
	{
		Config *c = Config::GetConfig();

		unordered_map<string, bool> dupJunctions_DEBUG;
		ofstream pileup(
			(c->outputDir + "/" + c->sampleName + ".fusions.pileup").c_str(), 
			ios::out);
		ofstream fusions(
			(c->outputDir + "/" + c->sampleName + ".fusions").c_str(), 
			ios::out);

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Compiling fusion results...";
		vector<Junction *> resultsJunctions;
		unordered_map<int, DiscordantCluster*> clustersMap;
		for (auto clust : clusters)
		{
			clustersMap[clust->clusterId] = clust;
			for (auto j : clust->junctions)
			{
				if (j->IsSpurious())
					continue;

				//NOTE: Check for special cases where junction.SeqA/B is < 50 bps
				string dupKey = lexical_cast<string>(clust->clusterId) + "_" +
					lexical_cast<string>(j->ex5p->exEnd) + "_" + 
					lexical_cast<string>(j->ex3p->exStart);
				
				//For smaller exons involved in the fusion, multiple isoforms
				//are considered to generate extended sequence; therefore, this is OK
				if (dupJunctions_DEBUG.find(dupKey) != dupJunctions_DEBUG.end())
					continue;   
				dupJunctions_DEBUG[dupKey] = true;

				//If the total number anchor reads is zero, don't proceed further
				if (j->GetTotalARCount() == 0) 
					continue;

				int totalARsWithMismatchesInAnchor = 0;
				for (auto ar : j->anchorReads){
					if (ar->isSpurious || !ar->IsAnchorRead())
						continue;
					if (ar->splitRead.mismatchInAnchor)
						totalARsWithMismatchesInAnchor++;
				}

				if (totalARsWithMismatchesInAnchor == j->GetTotalARCount() &&
					j->GetTotalARCount(true) == 0)
					continue;

				string s5 = j->ex5p->sequence, s3 = j->ex3p->sequence;
				if (s5.size() > 20 && s3.size() > 20) {
					double entropy5p =
						GetDinucleotideEntropy(s5.substr(s5.length() - 20, 20));
					double entropy3p = GetDinucleotideEntropy(s3.substr(0, 20));

					//Based on simulation of ~10,000 random junctions, need to expand this
					//to the entire transcriptome;
					if (entropy5p < 2 || entropy3p < 2)
						continue;

					if (entropy5p < 2.5 && entropy3p < 2.5 &&
						j->GetTotalARCount(true) == 0)
						continue;
				}

				//Mark a junction as spurious if there are no unique spanning reads 
				//supporting the fusion, and no high confidence anchor reads
				if (clust->GetNumDiscordReads(true) == 0 && 
					j->GetTotalARCount(true) == 0)
					continue;

				resultsJunctions.push_back(j);
			}
		}

		//Some anchor read are assigned to multiple fusion junctions; this
		//is mainly a consequence of homology at the fusion junction between 
		//genes A::B and genes A::C.  C and B can be homologous.  To avoid 
		//calling multiple fusions, assign the anchor read to the fusion 
		//junction with highest number of unqiue discordant reads, # of high 
		//conf anchor reads or just # of anchor reads
		try {
			unordered_map<string, vector<Junction *> > arToJ;
			for (auto j : resultsJunctions) {
				for (auto a : j->anchorReads) {
					if (a->isSpurious)
						continue;
					string rnm = a->splitRead.GetTrimmedReadName();
					arToJ[a->splitRead.GetTrimmedReadName()].push_back(j);
				}
			}
			for (auto arj : arToJ){
				auto vectJ = arj.second;
				if (vectJ.size() < 2)
					continue;
				std::sort(vectJ.begin(), vectJ.end(), Junction::compareJunctions);
				Junction *topJ = vectJ[0];
				for (int j = 1; j < vectJ.size(); j++)
					if (topJ != vectJ[j])
						vectJ[j]->MarkARspurious(arj.first);
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error while collapsing Junctions: " << e.what() << endl;
			exit(1);
		}

		int cnt = 0;
		for (auto r : resultsJunctions) {
			if (r->IsSpurious())
				continue;
			cnt++;
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tIdentified "
			<< cnt << " fusions passing filters";
		
		if (cnt == 0)
			return;

		if (resultsJunctions.size() > 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Compiling read mapping "
				<< "statistics for genes involved in fusions";
			FusionQuant::FQ.QuantifyMappingStats(resultsJunctions);
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tCompleted";
		}
		
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Preparing output...";
		
		std::sort(resultsJunctions.begin(), resultsJunctions.end(), Junction::compareJunctions);

		for (auto j : resultsJunctions) {
			if (j->IsSpurious())
				continue;
			JunctionStats stats = FusionQuant::FQ.CompileStatsForJunction(j);
			Gene *g5 = j->ex5p->gene, *g3 = j->ex3p->gene;
			Exon *e5 = j->ex5p, *e3 = j->ex3p;
			string seq5(j->c_seq5p), seq3(j->c_seq3p);
			//NOTE: Check for special cases where junction.SeqA/B is < 50 bps
			// ....

			double entropy5p =
				GetDinucleotideEntropy(seq5.substr(seq5.length() - 20, 20));
			double entropy3p = GetDinucleotideEntropy(seq3.substr(0, 20));

			// Write out fusion junction;
			stringstream lSS, rSS, outputSS;
			lSS << g5->name << ", (exon id:" << e5->exonId << "), "
				<< g5->chr << ":" << e5->exEnd_genomic;
			rSS << g3->name << ", (exon id:" << e3->exonId << "), "
				<< g3->chr << ":" << e3->exStart_genomic;
			outputSS << format("%80s--%-80s") % lSS.str() % rSS.str() << endl;
			outputSS << format("%80s  %-80s") % seq5 % seq3 << endl;

			bool mismatchesInAnchor = false;
			for (auto ar : j->anchorReads){
				if (ar->isSpurious)
					continue;
				SplitRead *sR = &ar->splitRead;
				OtherRead *oR = &ar->otherRead;

				string aS = sR->alignedSequence;
				int breakpoint = ar->junction->breakpoint - sR->position + 1;
				outputSS << format("%80s") % aS.substr(0, breakpoint) << "  ";
				outputSS << format("%-80s") %
					aS.substr(breakpoint, aS.length() - breakpoint);
				outputSS << " <-- SPLIT: (" << sR->ReadName << ", mismatches: "
					<< sR->mismatches << ") ";
				outputSS << " Other: (seq: " << oR->Sequence
					<< ", alignments: ";

				unordered_map<string, bool> avoidDups;
				for (auto orAl : oR->alignments) {
					string name = ((Gene*)orAl.gene)->name;
					if (avoidDups.find(name) != avoidDups.end())
						continue;
					outputSS << "[" << ((Gene*)orAl.gene)->name
						<< "," << orAl.mismatches << "], ";
					avoidDups[name] = true;
				}
				outputSS << ")" << endl;
			}
			pileup << outputSS.str() << endl;

			typedef FusionFragment::BreakPos BRK;
			BRK breakPos5p = BRK::UNKNOWN, breakPos3p = BRK::UNKNOWN;

			string codingSeqsCsv, transcriptSeqsCsv;
			stringstream isoAexpr, isoBexpr, isoABexpr;
			bool inFrame = false;
			for (auto ft : j->GetFusionTranscripts()){
				if (ft->IsInFrame())
					inFrame = true;
				string cds = ft->GetCodingSequence();
				string ts = ft->GetTranscriptSequence();
				string signature = "[" +
					ft->fiveP->fIsoform->isoformIdStr + "-" +
					ft->threeP->fIsoform->isoformIdStr + "]:";

				if (cds != "")
					codingSeqsCsv += signature + cds + ",";
				if (ts != "")
					transcriptSeqsCsv += signature + ts + ",";

				if (ft->fiveP->breakPosition > breakPos5p)
					breakPos5p = ft->fiveP->breakPosition;
				if (ft->threeP->breakPosition > breakPos3p)
					breakPos3p = ft->threeP->breakPosition;				
			}

			int ar10bp = j->GetTotalARCount(false, 10);
			int ar15bp = j->GetTotalARCount(false, 15);
			int ar20bp = j->GetTotalARCount(false, 20);

			auto clust = clustersMap[j->clusterId];
			fusions << g5->name << "_" << g3->name << "\t"
				<< clust->GetNumDiscordReads() << "\t"
				<< clust->GetNumDiscordReads(true) << "\t"
				<< j->GetTotalARCount() << "\t"
				<< j->GetTotalARCount(true) << "\t"
				<< g5->name << "\t"
				<< g5->chr << "\t"
				<< g5->strand << "\t"
				<< e5->exonId << "\t"
				<< (g5->strand == "+" ? e5->exEnd_genomic : e5->exStart_genomic) << "\t"
				<< FusionFragment::GetStringForBreakPos(breakPos5p) << "\t"
				<< g3->name << "\t"
				<< g3->chr << "\t"
				<< g3->strand << "\t"
				<< e3->exonId << "\t"
				<< (g3->strand == "+" ? e3->exStart_genomic : e3->exEnd_genomic) << "\t"
				<< FusionFragment::GetStringForBreakPos(breakPos3p) << "\t"
				<< (inFrame ? "1" : "0") << "\t"
				<< ar10bp - ar15bp << "\t"
				<< ar15bp - ar20bp << "\t"
				<< ar20bp << "\t"
				<< j->GetTotalARCountUniqueBy53(true) << "\t"
				<< j->GetTotalARCountUniqueBy53(false) << "\t"
				<< Gene::DistanceBetweenGenes(g5, g3) << "\t"
				<< entropy5p << "\t"
				<< entropy3p << "\t"
				<< stats.rpkmA << "\t"
				<< stats.rpkmB << "\t" 
				<< stats.rpkmA_5p << "\t"
				<< stats.rpkmA_3p << "\t"
				<< stats.rpkmB_5p << "\t"
				<< stats.rpkmB_3p << "\t"
				<< stats.numConcordsA_5p << "\t"
				<< stats.numConcordsA_3p << "\t"
				<< stats.numConcordsB_5p << "\t"
				<< stats.numConcordsB_3p << "\t"
				<< stats.numConcordsAA_span << "\t"
				<< stats.numConcordsAA_junct << "\t"
				<< stats.numConcordsBB_span << "\t"
				<< stats.numConcordsBB_junct << "\t"
				<< stats.numDiscords_A5p_B3p << "\t"
				<< stats.numDiscords_A3p_B5p << "\t"
				<< stats.numDiscords_A5p_B5p << "\t"
				<< stats.numDiscords_A3p_B3p << "\t"
				;

			fusions << codingSeqsCsv << "\t"
				<< transcriptSeqsCsv << "\t";
			fusions << endl;
		}
		pileup.close();
		fusions.close();
	}

	double FusionCompiler::GetDinucleotideEntropy(string sequence)
	{
		boost::to_upper(sequence);
		unordered_map<string, int> dinuc_counts;
		const char *bps_array[] = { "A", "T", "G", "C" };
		vector<string> bps(begin(bps_array), end(bps_array));
		for (vector<string>::size_type i = 0; i < bps.size(); i++) {
			for (vector<string>::size_type j = i + 1; j < bps.size(); j++) {
				dinuc_counts[bps[i] + bps[j]] = 0;
				dinuc_counts[bps[j] + bps[i]] = 0;
			}
		}
		
		for (string::size_type i = 0; i < sequence.length() - 1; i++)
			dinuc_counts[sequence.substr(i, 2)]++;

		double cum_prob = 0;
		for (auto iter : dinuc_counts){
			double prob = ((double)iter.second) / (sequence.length() - 1);
			if (prob > 0)
				cum_prob += prob * std::log(prob) / std::log((double)2);
		}
		return -1 * cum_prob;
	}
}
