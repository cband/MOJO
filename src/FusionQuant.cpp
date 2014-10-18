
#include "FusionQuant.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = ""));

	FusionQuant FusionQuant::FQ;

	PairedBamAlignment::PairedBamAlignment(BamAlignment a1, BamAlignment a2)
	{
		fAln.position = a1.Position;
		fAln.alignedLength = a1.Qualities.size();
		fAln.readId = stoi(Read::TrimReadName(a1.Name));
		fAln.end = 0;
		if (a1.Name[a1.Name.length() - 2] == '/')
			fAln.end = (a1.Name.back() == '1') ? 1 : 2;

		sAln.position = a2.Position;
		sAln.alignedLength = a2.Qualities.size();
		sAln.readId = stoi(Read::TrimReadName(a2.Name));
		sAln.end = 0;
		if (a2.Name[a2.Name.length() - 2] == '/')
			sAln.end = (a2.Name.back() == '1') ? 1 : 2;
	}

	FusionQuant::FusionQuant()
	{
		
	}

	string FusionQuant::GetPaddingSeq()
	{
		if (FusionQuant::paddingSeq.size() > 0)
			return FusionQuant::paddingSeq;
		for (int i = 0; i < PAD_LENGTH; i++)
			FusionQuant::paddingSeq += "G";
		return FusionQuant::paddingSeq;
	}

	// FusionQuant::QuantifyMappingStats is the main function of the FusionQuant
	// singleton class.  Given an input of all high confidence fusion junctions, 
	// this function constructs the transcriptomes for the select genes involved
	// in the fusions and determines mapping statistics for the fusion as well as 
	// the individual genes.  
	void FusionQuant::QuantifyMappingStats(vector<Junction* > junctions)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();
		Module modMappingStatus =
			Module::LoadStatusFile(c->statusFile_FusionStats);
		
		string mapFa = c->workingDir + "/mapstats/fusion.map.transcripts.fa";
		if (!modMappingStatus.IsComplete) {
			//Build fasta files for all isoforms of all genes in the junctions
			ofstream faStream(mapFa.c_str(), ios::out);
			unordered_map<int, bool> processed;
			for (auto j : junctions){
				Gene *g5p = j->ex5p->gene, *g3p = j->ex3p->gene;
				if (processed.find(g5p->geneId) == processed.end()) {
					for (auto iso : g5p->isoforms){
						iso->GetTranscript();
						faStream << ">" << iso->isoformIdStr << endl;
						faStream << GetPaddingSeq();
						faStream << iso->GetTranscript();
						faStream << GetPaddingSeq() << endl;
					}
				}

				if (processed.find(g3p->geneId) == processed.end()) {
					for (auto iso : g3p->isoforms){
						iso->GetTranscript();
						faStream << ">" << iso->isoformIdStr << endl;
						faStream << GetPaddingSeq();
						faStream << iso->GetTranscript();
						faStream << GetPaddingSeq() << endl;
					}
				}
			}

			//Build bwa index
			char cmd[25000];
			sprintf(cmd, "%s index %s", c->bwaPath.c_str(), mapFa.c_str());
			if (Utils::ExecuteCommand(cmd, "Main").exit_code != 0){
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "Error building fusion map index. ";
				exit(1);
			}

			//Run alignments in max of n-splits;
			vector<boost::thread *> threads;
			string alignmentFiles;
			try {
				ComputePerTask cpt = 
					ComputePerTask::CalculateComputePerTask(6, 2, 2);
				for (int threadId = 0; threadId < cpt.numSplits; threadId++) {
					threads.push_back(
						new boost::thread(FindFusionGeneMappingReads_worker,
						threadId, cpt, mapFa));
					alignmentFiles += mapFa + "_" + to_string(threadId) + ".bam ";
				}
				for (auto thread : threads) {
					thread->join();
					delete thread;
				}
				//Merge all alignments;
				if (cpt.numSplits == 1) {
					//just rename the file
					sprintf(cmd, "mv %s %s.bam",
						alignmentFiles.c_str(), mapFa.c_str());
					Utils::ExecuteCommand(cmd, "Main");
				}
				else {
					sprintf(cmd, "%s merge -f -1 %s.bam %s",
						c->samtoolsPath.c_str(), mapFa.c_str(), alignmentFiles.c_str());
					if (Utils::ExecuteCommand(cmd, "Main").exit_code != 0) {
						BOOST_LOG_CHANNEL(logger::get(), "Main")
							<< "Error sorting fusion map alignments. (samtools view/sort)";
						exit(1);
					}
				}
			}
			catch (std::exception &e) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occurred while "
					<< "extracting Unaligned Reads. Error: " << e.what();
			}


			//Determine the max number of cores for sorting steps
			ComputePerTask cpt = 
				ComputePerTask::CalculateComputePerTask(c->maxCores, 1, 1);
			sprintf(cmd, "%s sort -n -@ %d %s.bam %s.sorted",
				c->samtoolsPath.c_str(), cpt.numSplits, 
				mapFa.c_str(), mapFa.c_str());
			if (Utils::ExecuteCommand(cmd, "Main").exit_code != 0) {
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "Error sorting fusion map alignments. (samtools view/sort)";
				exit(1);
			}

			if (c->removeTemporaryFiles){
				Utils::DeleteFile(mapFa);
				Utils::DeleteFile(mapFa + ".bam");
				Utils::DeleteFile(mapFa + "_*");
				Utils::DeleteFiles(std::vector<string>{mapFa + ".amb", 
					mapFa + ".ann", mapFa + ".bwt", 
					mapFa + ".pac", mapFa + ".sa"});
			}
		}
		modMappingStatus.IsComplete = true;
		modMappingStatus.SaveStatusFile();

		BamReader reader;
		BamAlignment aln1, aln2;
		reader.Open((mapFa + ".sorted.bam").c_str());
		RefVector headerVect = reader.GetReferenceData();
		unordered_map<string, int> refNameToId;
		for (int rv = 0; rv < headerVect.size(); rv++)
			refNameToId[headerVect[rv].RefName] = rv;

		while (reader.GetNextAlignment(aln1)) {
			reader.GetNextAlignment(aln2);
			while (Read::TrimReadName(aln1.Name) != Read::TrimReadName(aln2.Name))
			{
				aln1 = aln2;
				reader.GetNextAlignment(aln2);
			}
			bool isConcordant = (aln1.RefID == aln2.RefID) ? true : false;
			string readname = Read::TrimReadName(aln1.Name);
			auto aln1_vect =
				DiscordantReadFinder::GetBwaAlignsFromXA(aln1, &refNameToId);
			auto aln2_vect =
				DiscordantReadFinder::GetBwaAlignsFromXA(aln2, &refNameToId);
			
			unordered_map<string, bool> doneLoading;
			for (auto a : aln1_vect) {
				for (auto b : aln2_vect){
					if (isConcordant && a.RefID != b.RefID)
						continue;

					string a_id_str = to_string(a.RefID);
					string b_id_str = to_string(b.RefID);
					string key = a_id_str + "_" + b_id_str;
					if (a.RefID > b.RefID)
						key = b_id_str + "_" + a_id_str;
					if (doneLoading.find(key) != doneLoading.end())
						continue;
					doneLoading[key] = true;

					BamAlnPair p(a, b);
					PairedBamAlignment *pa = new PairedBamAlignment(a, b);
					pa->fExon = FindExonsMappedByRead(a, headerVect[a.RefID].RefName);
					pa->sExon = FindExonsMappedByRead(b, headerVect[b.RefID].RefName);
					pa->fIso = gm->IsoformsMap[headerVect[a.RefID].RefName];
					pa->sIso = gm->IsoformsMap[headerVect[a.RefID].RefName];
					if (pa->fExon == 0 || pa->sExon == 0)
						continue;
					Isoform *isoA = gm->IsoformsMap[headerVect[a.RefID].RefName];
					Isoform *isoB = gm->IsoformsMap[headerVect[b.RefID].RefName];

					string gKey = to_string(isoA->gene->geneId) + "_" +
						to_string(isoB->gene->geneId);
					if (isoA->gene->geneId > isoB->gene->geneId) {
						gKey = to_string(isoB->gene->geneId) + "_" +
							to_string(isoA->gene->geneId);
					}

					if (isoA->gene->geneId == isoB->gene->geneId)
						geneAlnsConcordant[to_string(isoA->gene->geneId)].push_back(pa);
					else
						geneAlnsDiscordant[gKey].push_back(pa);
				}
			}
		}
	}

	void FusionQuant::FindFusionGeneMappingReads_worker(int threadId,
		ComputePerTask cpt, string mapFa)
	{
		Config *c = Config::GetConfig();
		try
		{
			string end1Fifo = Config::GenerateCatCmdsForSplit(1, threadId, 
				cpt.numSplits, true);
			string end2Fifo = Config::GenerateCatCmdsForSplit(2, threadId, 
				cpt.numSplits, true);
			string end1Fq = c->workingDir + "/fastqs/split_" + 
				to_string(threadId) + "_1.fastq";
			string end2Fq = c->workingDir + "/fastqs/split_" + 
				to_string(threadId) + "_2.fastq";

			char cmd[25000];
			//Align all reads to this 'fusion' transcriptome
			sprintf(cmd, "%s \n %s aln -q 15 -R 100 -t %d %s %s > " 
				"%s_%d_aln_1.sai 2> %s_%d_aln_1_output.log ",
				end1Fifo.c_str(), c->bwaPath.c_str(), c->maxCores, mapFa.c_str(), 
				end1Fq.c_str(), mapFa.c_str(), threadId, mapFa.c_str(), threadId);
			Utils::ExecuteCommand(cmd, "Main", true, EXIT_ON_FAIL);

			sprintf(cmd, "%s \n %s aln -q 15 -R 100 -t %d %s %s > "
				"%s_%d_aln_2.sai 2> %s_%d_aln_2_output.log ",
				end2Fifo.c_str(), c->bwaPath.c_str(), c->maxCores, mapFa.c_str(), 
				end2Fq.c_str(), mapFa.c_str(), threadId, mapFa.c_str(), threadId);
			Utils::ExecuteCommand(cmd, "Main", true, EXIT_ON_FAIL);

			sprintf(cmd, "%s \n%s \n%s sampe -A -a 1000 -n 250 -N 250 -c 0.0001" 
				" -P %s %s_%d_aln_1.sai %s_%d_aln_2.sai %s %s 2> "
				" %s_%d_sampe_output.log |  %s view -F 12 -b -S - -o %s_%d.bam ",
				end1Fifo.c_str(), end2Fifo.c_str(), c->bwaPath.c_str(), 
				mapFa.c_str(), mapFa.c_str(), threadId, mapFa.c_str(), threadId, 
				end1Fq.c_str(), end2Fq.c_str(), mapFa.c_str(), threadId, 
				c->samtoolsPath.c_str(), mapFa.c_str(), threadId);
			Utils::ExecuteCommand(cmd, "Main", true, EXIT_ON_FAIL);
		}
		catch (std::exception &e)
		{
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occurred while "
				<< "extracting reads mapping to fusion genes. Exception: " 
				<< e.what();
			exit(1);
		}
	}

	Exon* FusionQuant::FindExonsMappedByRead(BamAlignment aln, string name)
	{
		GeneModel *gm = GeneModel::GetGeneModel();
		int start = aln.Position - PAD_LENGTH;
		int end = aln.Position + aln.Qualities.size() - PAD_LENGTH;

		Isoform *iso = gm->IsoformsMap[name];
		if (iso->GetStrand() == "-") {
			std::swap(start, end);
			start = iso->GetTranscriptLength() - start;
			end = iso->GetTranscriptLength() - end;
		}
		int midPoint = (start + end) / 2;
		vector<Exon *> exons;
		int carat = 0;
		for (vector<Exon *>::size_type i = 0; i < iso->exons.size(); i++) {
			int exonLen = (int)iso->exons[i]->sequence.size();
			if (carat <= midPoint && (carat + exonLen) >= midPoint)
				return iso->exons[i];
			carat += exonLen;
		}
		return 0;
	}
	
	JunctionStats FusionQuant::CompileStatsForJunction(Junction *j)
	{
		JunctionStats stats;
		
		Gene *gA = j->ex5p->gene, *gB = j->ex3p->gene;
		unordered_map<int, bool> upExonsA, downExonsA, upExonsB, downExonsB;

		bool found = false;
		for (auto e : gA->allExons) {
			if ((!found && gA->strand == "+") || (found && gA->strand == "-"))
				upExonsA[e->exonId] = true;
			else
				downExonsA[e->exonId] = true;
			if (e == j->ex5p) {
				upExonsA[e->exonId] = true;
				downExonsA[e->exonId] = true;
				found = true;
			}
		}

		found = false;
		for (auto e : gB->allExons) {
			if ((!found && gB->strand == "+") || (found && gB->strand == "-"))
				upExonsB[e->exonId] = true;
			else
				downExonsB[e->exonId] = true;
			if (e == j->ex3p) {
				upExonsB[e->exonId] = true;
				downExonsB[e->exonId] = true;
				found = true;
			}
		}
		downExonsA.erase(j->ex5p->exonId);
		upExonsB.erase(j->ex3p->exonId);

		unordered_map<int, unordered_map<int, bool> > exonToRead;
		unordered_map<int, bool> concordAA_JunctReads, concordBB_JunctReads;
		unordered_map<int, bool> concordAA_SpanReads, concordBB_SpanReads;
		for (auto pba : geneAlnsConcordant[to_string(gA->geneId)]){
			exonToRead[pba->fExon->exonId][pba->GetReadId()] = true;
			exonToRead[pba->sExon->exonId][pba->GetReadId()] = true;

			if ((upExonsA.find(pba->fExon->exonId) != upExonsA.end() &&
				downExonsA.find(pba->sExon->exonId) != downExonsA.end()) ||
				(upExonsA.find(pba->sExon->exonId) != upExonsA.end() &&
				downExonsA.find(pba->fExon->exonId) != downExonsA.end()))
			{
				concordAA_SpanReads[pba->GetReadId()] = true;
			}
			//Determine if the read maps to the proper junction of the fusion gene
			if (pba->fExon == j->ex5p) {
				if (DoesExonMapToJunction(pba->fIso, j->ex5p, 
					pba->fAln, FRAG_TYPE::FIVE_P))
					concordAA_JunctReads[pba->GetReadId()] = true;
			}
			if (pba->sExon == j->ex5p) {
				if (DoesExonMapToJunction(pba->sIso, j->ex5p, 
					pba->sAln, FRAG_TYPE::FIVE_P))
					concordAA_JunctReads[pba->GetReadId()] = true;
			}
		}
		for (auto pba : geneAlnsConcordant[to_string(gB->geneId)]){
			exonToRead[pba->fExon->exonId][pba->GetReadId()] = true;
			exonToRead[pba->sExon->exonId][pba->GetReadId()] = true;

			if ((upExonsB.find(pba->fExon->exonId) != upExonsB.end() &&
				downExonsB.find(pba->sExon->exonId) != downExonsB.end()) ||
				(upExonsB.find(pba->sExon->exonId) != upExonsB.end() &&
				downExonsB.find(pba->fExon->exonId) != downExonsB.end()))
			{
				concordBB_SpanReads[pba->GetReadId()] = true;
			}

			if (pba->fExon == j->ex3p) {
				if (DoesExonMapToJunction(pba->fIso, j->ex3p, 
					pba->fAln, FRAG_TYPE::THREE_P))
					concordBB_JunctReads[pba->GetReadId()] = true;
			}
			if (pba->sExon == j->ex3p) {
				if (DoesExonMapToJunction(pba->sIso, j->ex3p, 
					pba->sAln, FRAG_TYPE::THREE_P))
					concordBB_JunctReads[pba->GetReadId()] = true;
			}
		}

		string key = to_string(gA->geneId) + "_" + to_string(gB->geneId);
		if ( gA->geneId > gB->geneId)
			key = to_string(gB->geneId) + "_" + to_string(gA->geneId);
		for (auto pba : geneAlnsDiscordant[key]){
			exonToRead[pba->fExon->exonId][pba->GetReadId()] = true;
			exonToRead[pba->sExon->exonId][pba->GetReadId()] = true;
		}

		unordered_map<int, bool> upReadsA, downReadsA, upReadsB, downReadsB;
		for (auto eTR : exonToRead) {
			for (auto read : eTR.second) {
				if (upExonsA.find(eTR.first) != upExonsA.end())
					upReadsA[read.first] = true;
				else if (downExonsA.find(eTR.first) != downExonsA.end())
					downReadsA[read.first] = true;

				if (upExonsB.find(eTR.first) != upExonsB.end())
					upReadsB[read.first] = true;
				else if (downExonsB.find(eTR.first) != downExonsB.end())
					downReadsB[read.first] = true;
			}
		}
		stats.numConcordsA_5p = upReadsA.size();		// paired-read count; 
		stats.numConcordsA_3p = downReadsA.size();
		stats.numConcordsB_5p = upReadsB.size();
		stats.numConcordsB_3p = downReadsB.size();
		stats.numConcordsAA_span = concordAA_SpanReads.size();
		stats.numConcordsBB_span = concordBB_SpanReads.size();
		stats.numConcordsAA_junct = concordAA_JunctReads.size();
		stats.numConcordsBB_junct = concordBB_JunctReads.size();
		
		for (auto a : upReadsA) {
			if (upReadsB.find(a.first) != upReadsB.end()) {
				stats.numDiscords_A5p_B5p++;
				stats.numConcordsA_5p--;
				stats.numConcordsB_5p--;
			}
		}
		for (auto a : upReadsA) {
			if (downReadsB.find(a.first) != downReadsB.end()) {
				stats.numDiscords_A5p_B3p++;
				stats.numConcordsA_5p--;
				stats.numConcordsB_3p--;
			}
		}
		for (auto a : downReadsA) {
			if (upReadsB.find(a.first) != upReadsB.end()) {
				stats.numDiscords_A3p_B5p++;
				stats.numConcordsA_3p--;
				stats.numConcordsB_5p--;
			}
		}
		for (auto a : downReadsA) {
			if (downReadsB.find(a.first) != downReadsB.end()) {
				stats.numDiscords_A3p_B3p++;
				stats.numConcordsA_3p--;
				stats.numConcordsB_3p--;
			}
		}

		int gA_5p_len = gA->GetPartialTranscript(j->ex5p, FRAG_TYPE::FIVE_P).length();
		int gA_3p_len = gA->GetPartialTranscript(j->ex5p, FRAG_TYPE::THREE_P).length();
		int gB_5p_len = gB->GetPartialTranscript(j->ex3p, FRAG_TYPE::FIVE_P).length();
		int gB_3p_len = gB->GetPartialTranscript(j->ex3p, FRAG_TYPE::THREE_P).length();
		int gA_len = gA_5p_len + gA_3p_len;
		int gB_len = gB_5p_len + gB_3p_len;

		Config *c = Config::GetConfig();
		int readcount = c->GetTotalReadcount();
		
		double factor = (1.0 / readcount) * 1e9;
		stats.rpkmA = 
			(( 1.0 * (stats.numConcordsA_5p + stats.numConcordsA_3p)) / gA_len) * factor;
		stats.rpkmB = 
			(( 1.0 * (stats.numConcordsB_5p + stats.numConcordsB_3p)) / gB_len) * factor;
		stats.rpkmA_5p = (1.0 * stats.numConcordsA_5p / gA_5p_len) * factor;
		stats.rpkmA_3p = (1.0 * stats.numConcordsA_3p / gA_3p_len) * factor;
		stats.rpkmB_5p = (1.0 * stats.numConcordsB_5p / gB_5p_len) * factor;
		stats.rpkmB_3p = (1.0 * stats.numConcordsB_3p / gB_3p_len) * factor;
		return stats;
	}

	bool FusionQuant::DoesExonMapToJunction(Isoform *iso, Exon *junctExon,
		BamAlignmentEnd aln, FRAG_TYPE type)
	{
		//check if "aln" maps to exon
		bool found = false;
		int length = 0;
		for (auto e : iso->exons) {
			if (junctExon == e) {
				found = true;
				if (iso->GetStrand() == "+" && type == FRAG_TYPE::FIVE_P)
					length += e->sequence.size();
				if (iso->GetStrand() == "-" && type == FRAG_TYPE::THREE_P)
					length += e->sequence.size();
			}
			if (found)
				break;
			length += e->sequence.size();
		}

		int start = aln.position - PAD_LENGTH;
		int end = aln.position + aln.alignedLength - PAD_LENGTH;

		if (iso->GetStrand() == "-") {
			std::swap(start, end);
			start = iso->GetTranscriptLength() - start;
			end = iso->GetTranscriptLength() - end;
		}
		if (start < length && end > length)
			return true;
		return false;
	}
}