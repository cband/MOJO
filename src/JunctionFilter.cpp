
#include "JunctionFilter.h"

#define MAX_EXOME_SPLITS 12

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = ""));

	//
	// JunctionFilter::Run is the main function of this singleton class.  Takes 
	// discordant clusters with fusion junctions and anchor reads loaded, as input.
	// Filters out anchor reads that map canonically to an alternate region of the 
	// genome and transcriptome (both novel and known).  
	//
	void JunctionFilter::Run(vector<DiscordantCluster *> clusters)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		boost::filesystem::create_directories(c->JunctionFilter_dir);
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Filtering Anchor Reads...";
		AnchorReadMap arMap;
		for (auto cluster : clusters)
			for (auto junct : cluster->junctions)
				for (auto ar : junct->anchorReads)
					arMap[ar->GetSplitReadIdentifier()] = ar;

		Module modStatus = Module::LoadStatusFile(c->statusFile_JunctionFilter);
		
		//First align all anchor reads to the genome 
		string step_fa = c->JunctionFilter_dir + "/step0.fa";
		GenerateFastaForAnchorReads(&arMap, step_fa);
		if (!modStatus.IsComplete) {
			char cmd[3000];
			sprintf(cmd, "%s -p %d -x %s --very-sensitive --no-unal -f %s |"
				" %s view -b -S - -o %s.bam", c->bowtie2Path.c_str(),
				min(12, c->maxCores), c->bowtie2GenomeIndex.c_str(), 
				step_fa.c_str(), c->samtoolsPath.c_str(), step_fa.c_str());
			if (!modStatus.IsComplete &&
				Utils::ExecuteCommand(cmd).exit_code != 0)
			{
				BOOST_LOG_CHANNEL(logger::get(), "Main") 
					<< "Error occurred while " << "executing the command: " << cmd;
				exit(0);
			}
		}
	
		// Iterate through the gneome alignments and flag the anchor read as 'spurious'
		// if the split-end read maps in full to an alternate region of the genome.
		BamTools::BamReader reader;
		BamTools::BamAlignment al;
		reader.Open(step_fa + ".bam");
		while (reader.GetNextAlignment(al)) {
			float num_matches = 0;
			for (int i = 0; i < al.CigarData.size(); i++)
				if (al.CigarData[i].Type == 'M')
					num_matches += al.CigarData[i].Length;
			int mismatches = 0;;
			al.GetTag("NM", mismatches);
			if ((num_matches - mismatches)/al.QueryBases.length() > 0.95) 
				arMap[al.Name]->isSpurious = true;
		}

		// The anchor reads that passed the previous filter are aligned to the exome and
		// genome using blat in iterative steps to maximize sensitivity.  
		FilterStringency stringen[2] = { LOW, HIGH };
		for (int step = 0; step < 2; step++) {
			string faPfx = c->JunctionFilter_dir + "/step" + 
				lexical_cast<string>(step);

			string step_fa = faPfx + "_exome.fa";
			if (!modStatus.IsComplete)
				GenerateFastaForAnchorReads(&arMap, step_fa);
			FilterSplitReadsByExome(step_fa, &arMap, stringen[step], &modStatus);

			step_fa = faPfx + "_genome.fa";
			if (!modStatus.IsComplete)
				GenerateFastaForAnchorReads(&arMap, step_fa);
			FilterSplitReadsByGenome(step_fa, &arMap, stringen[step], &modStatus);

			step_fa = faPfx + "_final.fa";
			if (!modStatus.IsComplete)
				GenerateFastaForAnchorReads(&arMap, step_fa);


			//Filter out fusion junctions if the 80bp fusion junction can be aligned
			//concordantly to anywhere else in the genome.
			if (step == 0)
				FilterSpuriousJunctions(&arMap, &modStatus);
		}
		FilterSpuriousScrambledReads(&arMap, &modStatus);
		if (!modStatus.IsComplete)
			GenerateFastaForAnchorReads(&arMap, c->JunctionFilter_dir +
				"/step_post_scramble.fa");
		FilterJunctionPCRdups(&arMap);
		if (!modStatus.IsComplete)
			GenerateFastaForAnchorReads(&arMap, c->JunctionFilter_dir +
				"/step_post_pcr.fa");

		modStatus.IsComplete = true;
		modStatus.SaveStatusFile();

		int anchorCount = 0;
		for (auto am : arMap)
			if (!am.second->isSpurious)
				anchorCount++;
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "\t" << anchorCount 
			<< " anchor reads passed filtering steps";
	}

	//Outputs all anchor reads that are not designated as spurious to the 'outputFileName'
	void JunctionFilter::GenerateFastaForAnchorReads
		(AnchorReadMap *arMap, string outputFilename)
	{
		ofstream outFile(outputFilename.c_str(), ios::out);
		for (auto ar : (*arMap)) {
			if (ar.second->isSpurious)
				continue;
			outFile << ">" << ar.second->GetSplitReadIdentifier() << endl;
			outFile << ar.second->splitRead.alignedSequence << endl;
		}
		outFile.close();
	}
	
	void JunctionFilter::FilterSplitReadsByExome(string faFile,
		AnchorReadMap *arMap, FilterStringency stringency, Module *modStatus)
	{
		Config *c = Config::GetConfig();
		vector<boost::thread *> threads;
		if (!modStatus->IsComplete) {
			modStatus->NumSplits = 
				ComputePerTask::CalculateComputePerTask(MAX_EXOME_SPLITS, 1, 2).numSplits;
			int numReads = Utils::LineCount(faFile) / 2;

			if (numReads == 0)
				return;

			if (numReads < 1000)
				modStatus->NumSplits = 1;

			int linesPerSplit =
				((numReads + modStatus->NumSplits - 1) / modStatus->NumSplits) * 2;
			char cmd[1000];
			sprintf(cmd, "split -l %d -d %s %s.", linesPerSplit,
				faFile.c_str(), faFile.c_str());
			if (Utils::ExecuteCommand(cmd).exit_code != 0) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "Unable to split file: "
					<< faFile << ". Exiting";
				exit(1);
			}

			for (int i = 0; i < modStatus->NumSplits; i++) {
				string faFileSp = faFile + "." + (boost::format("%02d") % i).str();
				threads.push_back(new boost::thread(FilterSplitReadsByExome_worker,
					faFileSp, stringency));
			}

			for (int i = 0; i < threads.size(); i++) {
				threads[i]->join();
				delete threads[i];
			}
		}

		for (int i = 0; i < MAX_EXOME_SPLITS; i++){
			string pslFileName = faFile + "." + 
				(boost::format("%02d") % i).str() + ".psl";
			
			if (!Utils::FileExists(pslFileName))
				continue;

			PSLParser pslParser(pslFileName);
			AnchorReadPSLLine psl;
			while (pslParser.GetNextLine(psl)) {
				int b_qStart = psl.qStart, b_qEnd = psl.qEnd;
				if (psl.blockCount > 1) {
					int max_block_size = -1, max_block_size_pos = -1;
					b_qStart = 0, b_qEnd = 0;
					for (int d = 0; d < psl.blockCount; d++)
						if (max_block_size < psl.blockSizes[d])
							max_block_size = psl.blockSizes[d], 
							max_block_size_pos = d;

					int max_block_to_include = max_block_size_pos;
					int min_block_to_include = max_block_size_pos;
					// Check for contigs before the "max block size";
					for (int bef = max_block_size_pos - 1; bef >= 0; bef--) {
						int gap = psl.tStarts[bef + 1] - 
							(psl.tStarts[bef] + psl.blockSizes[bef]);
						if (psl.blockSizes[bef] < 8 && psl.blockSizes[bef] < 2 * gap)
							break;
						min_block_to_include = bef;
					}

					for (int aft = max_block_size_pos + 1; 
						aft < psl.blockSizes.size(); aft++) 
					{
						int gap = psl.tStarts[aft] - 
							psl.tStarts[aft - 1] + psl.blockSizes[aft - 1];
						if (psl.blockSizes[aft] < 8 && 
								psl.blockSizes[aft] < 2 * gap)
							break;
						max_block_to_include = aft;
					}
					b_qStart = psl.qStarts[min_block_to_include];
					for (int mm = min_block_to_include; 
						mm <= max_block_to_include; mm++)
						b_qEnd += psl.blockSizes[mm];

					if (psl.strand == "-") {
						int size = abs(b_qStart - b_qEnd);
						b_qStart = psl.qSize - b_qStart - psl.matches;
						b_qEnd = b_qStart + size;
					}
					if ((b_qEnd - b_qStart + 1) < psl.qSize) continue;
				}

				if ((psl.break_at - b_qStart) >= 10 && 
					(b_qEnd - psl.break_at) >= 10)
					(*arMap)[psl.qName]->isSpurious = true;
			}
		}
	}

	void JunctionFilter::FilterSplitReadsByExome_worker(string faFile,
		FilterStringency stringency)
	{
		Config *c = Config::GetConfig();
		string params = JunctionFilter::GetBlatParamsForStringency(stringency);
		char cmd[1000];
		sprintf(cmd, "%s %s %s %s %s.psl", c->blatPath.c_str(), params.c_str(), 
			c->blatReferenceJunctions2Bit.c_str(), faFile.c_str(), faFile.c_str());
		if (Utils::ExecuteCommand(cmd).exit_code != 0){
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error executing command: "
				<< cmd << ", Exiting";
			exit(1);
		}
	}

	void JunctionFilter::FilterSplitReadsByGenome(string faFile,
		AnchorReadMap *arMap, FilterStringency stringency, Module *modStatus)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		//initialize variables for code readability
		vector<string> contigFaFiles, contigPSLFiles;
		for (string contig : c->blatFilterChromsVect)
			contigFaFiles.push_back(faFile + "." + contig + ".fa");

		if (modStatus->IsComplete) {
			for (string contig : c->blatFilterChromsVect) {
				FilterSplitReadsByGenome_post(faFile + "." + contig + ".fa"
					+ ".psl", arMap);
			}
			return;
		}

		ComputePerTask cpt = ComputePerTask::CalculateComputePerTask(12, 1, 2);
		vector<boost::thread *> threads;
		vector<bool> thread_finished;
		int finishedJobs = 0;
		bool queueFull = false;
		for (vector<string>::size_type i = 0; 
			i < c->blatFilterChromsVect.size(); 
			i++) 
		{
			string threadIdStr = lexical_cast<string>(i);
			GenerateFastaForAnchorReads(arMap, contigFaFiles[i]);
			threads.push_back(new boost::thread(FilterSplitReadsByGenome_worker,
				contigFaFiles[i], c->blatFilterChromsVect[i], stringency));
			thread_finished.push_back(false);

			if ((i + 1) == cpt.numSplits)
				queueFull = true;
			for (int j = 0; queueFull && j <= i; j++) {
				if (!thread_finished[j]
					&& threads[j]->try_join_for(boost::chrono::milliseconds(10))) 
				{
					thread_finished[j] = true;
					FilterSplitReadsByGenome_post(contigFaFiles[j] + ".psl", arMap);
					break;
				}
				if (j == i) j = 0;
				boost::this_thread::sleep(boost::posix_time::milliseconds(10));
			}
		}

		for (vector<string>::size_type j = 0; 
			j < c->blatFilterChromsVect.size(); j++) 
		{
			threads[j]->join();
			delete threads[j];
			FilterSplitReadsByGenome_post(contigFaFiles[j] + ".psl", arMap);
		}
	}
	
	void JunctionFilter::FilterSplitReadsByGenome_worker(string faFile, 
		string contig, FilterStringency stringency)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();
		string params = JunctionFilter::GetBlatParamsForStringency(stringency);
		char cmd[1000];
		sprintf(cmd, "%s %s %s/%s.2bit %s %s.psl",
			c->blatPath.c_str(), params.c_str(), c->blatFilterChromsDir.c_str(), 
			contig.c_str(), faFile.c_str(), faFile.c_str(), contig.c_str());
		if (Utils::ExecuteCommand(cmd).exit_code != 0){
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error executing command: "
				<< cmd << ", Exiting";
			exit(1);
		}
	}

	void JunctionFilter::FilterSplitReadsByGenome_post(string pslFile,
		AnchorReadMap *arMap)
	{
		GeneModel *gm = GeneModel::GetGeneModel(); PSLParser pslParser(pslFile);
		AnchorReadPSLLine psl;
		while (pslParser.GetNextLine(psl)) {
			if ((psl.break_at - psl.qStart > 20 && psl.qEnd - psl.break_at > 20)
				|| (psl.qSize - psl.matches) < psl.qSize * 0.15)
			{
				bool gAfound = false, gBfound = false;
				for (int d = 0; d < psl.blockCount; d++) {
					string fetched_gene = GeneModel::GetGeneNameForCoordinates(
						psl.tName, 
						(int)(psl.tStarts[d] + (psl.blockSizes[d] / 2)));

					if (fetched_gene != "-") {
						auto sp = Utils::SplitToVector(psl.qName, "_");
						Exon *exA = gm->ExonsMap[stoi(sp[1])];
						Exon *exB = gm->ExonsMap[stoi(sp[2])];
						if (fetched_gene == exA->gene->name) gAfound = true;
						if (fetched_gene == exB->gene->name) gBfound = true;
					}
				}
				if (!(gAfound && gBfound) && psl.misMatches <= 2)
					(*arMap)[psl.qName]->isSpurious = true;
			}
		}
	}

	void JunctionFilter::GenerateFastaForJunctions(AnchorReadMap *arMap, 
		string outputFa)
	{
		JunctionMap junctionsMap;
		ofstream jf(outputFa.c_str(), ios::out);
		for (auto ar : (*arMap)){
			Junction *j = ar.second->junction;
			if (j->IsSpurious())
				continue;
			if (junctionsMap.find(j->GetJunctionName()) != junctionsMap.end())
				continue;
			
			junctionsMap[j->GetJunctionName()] = j;
			string seq5p(j->c_seq5p), seq3p(j->c_seq3p);
			jf << ">" << j->GetJunctionName() << "_1" << endl;
			jf << seq5p << endl;
			jf << ">" << j->GetJunctionName() << "_2" << endl;
			jf << seq3p << endl;
			if (seq5p.length() < 30 || seq3p.length() < 30)
				continue;
			if (seq5p.length() > 40)
				seq5p = seq5p.substr(seq5p.length() - 40, 40);
			if (seq3p.length() > 40)
				seq3p = seq3p.substr(0, 40);

			jf << ">" << j->GetJunctionName() << "_J" << endl
				<< seq5p << seq3p << endl;
		}
		jf.close();
	}

	void JunctionFilter::FilterSpuriousJunctions(AnchorReadMap *arMap, 
		Module *modStatus)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		string fusJunctsFile = c->JunctionFilter_dir + "/fusions.junctions";
		//initialize variables for code readability
		vector<string> contigFaFiles, contigPSLFiles;
		for (string contig : c->blatFilterChromsVect) 
			contigFaFiles.push_back(fusJunctsFile + "." + contig + ".fa");

		if (modStatus->IsComplete) {
			for (string contig : c->blatFilterChromsVect) {
				FilterSpuriousJunctions_post(fusJunctsFile + "." + contig + ".fa"
					+ ".psl", arMap);
			}
			return;
		}
		ComputePerTask cpt = ComputePerTask::CalculateComputePerTask(12, 1, 2);
		vector<boost::thread *> threads;
		vector<bool> thread_finished;
		int finishedJobs = 0;
		bool queueFull = false;

		for (vector<string>::size_type i = 0; 
			i < c->blatFilterChromsVect.size(); i++) 
		{
			string threadIdStr = lexical_cast<string>(i);
			GenerateFastaForJunctions(arMap, contigFaFiles[i]);
			threads.push_back(new boost::thread(FilterSplitReadsByGenome_worker,
				contigFaFiles[i], c->blatFilterChromsVect[i], 
				FilterStringency::HIGH));
			thread_finished.push_back(false);

			if ((i + 1) == cpt.numSplits)
				queueFull = true;
			for (int j = 0; queueFull && j <= i; j++) {
				if (!thread_finished[j] && 
					threads[j]->try_join_for(boost::chrono::milliseconds(10))) {
					thread_finished[j] = true;
					FilterSpuriousJunctions_post(contigFaFiles[j] + ".psl", arMap);
					break;
				}
				if (j == i) j = 0;
				boost::this_thread::sleep(boost::posix_time::milliseconds(10));
			}
		}
		for (vector<string>::size_type j = 0;
			j < c->blatFilterChromsVect.size(); j++)
		{
			threads[j]->join();
			delete threads[j];
			FilterSpuriousJunctions_post(contigFaFiles[j] + ".psl", arMap);
		}
	}
	
	void JunctionFilter::FilterSpuriousJunctions_post(string pslFile,
		AnchorReadMap *arMap)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel(); 

		JunctionMap junctionsMap;
		for (auto ar : (*arMap)){
			Junction *j = ar.second->junction;
			if (j->IsSpurious())
				continue;
			if (junctionsMap.find(j->GetJunctionName()) != junctionsMap.end())
				continue;
			junctionsMap[j->GetJunctionName()] = j;
		}

		unordered_map<string, vector<string> > forwardJ, reverseJ;
		PSLParser pslParser(pslFile);
		JunctionPSLLine psl;
		while (pslParser.GetNextLine(psl)) {
			string jEnd = psl.qName.substr(psl.qName.length() - 1, 1);
			string jName = psl.qName.substr(0, psl.qName.size() - 2);

			if (junctionsMap.find(jName) == junctionsMap.end())
				continue;

			map<string, bool> same_gene;
			bool gAfound = false, gBfound = false;
			if (((double)psl.misMatches) / psl.matches < 0.1) {
				auto sp = Utils::SplitToVector(psl.qName, "_");
				Exon *exA = gm->ExonsMap[stoi(sp[1])];
				Exon *exB = gm->ExonsMap[stoi(sp[2])];
				//find the middle position of the largest block to use as anchor to search for gene it is in
				int max_size = 0, max_size_pos = 0, prev_end = 0;
				int blockCount = psl.blockCount;
				for (int d = 0; d < psl.blockCount; d++) {
					int bSize = (int)psl.blockSizes[d];
					int tStart = (int)psl.tStarts[d];
					if (abs(prev_end - tStart) <= 5)
						blockCount--;
					prev_end = tStart + bSize;
					if (max_size < bSize) {
						max_size = bSize;
						max_size_pos = (int)(tStart + (bSize / 2));
					}
					string fetched_gene = GeneModel::GetGeneNameForCoordinates(
						psl.tName, (int)(tStart + (bSize / 2)));
					if (fetched_gene != "-") {
						if (fetched_gene == exA->gene->name) gAfound = true;
						if (fetched_gene == exB->gene->name) gBfound = true;
					}
				}

				if (jEnd == "J" && psl.misMatches <= 2) {
					if (gAfound && gBfound) {
						int dAB = Gene::DistanceBetweenGenes(exA->gene, exB->gene);
						if (dAB < c->readThruDist)
							junctionsMap[jName]->SetAsSpurious();
					}
					else if (((double)psl.matches) / psl.qSize > 0.9) {
						junctionsMap[jName]->SetAsSpurious();
					}
				}
				
				//since we are looking at only 40bp pairs, >2 may lead to false negativs
				//allow upto 2bp indels in each block
				if (psl.blockCount > 2)
					continue;

				if (jEnd == "1")
					forwardJ[jName].push_back(
					GeneModel::GetGeneNameForCoordinates(psl.tName, max_size_pos));
				else if (jEnd == "2")
					reverseJ[jName].push_back(
					GeneModel::GetGeneNameForCoordinates(psl.tName, max_size_pos));
			}
		}

		for (auto fJ : forwardJ){
			string jName = fJ.first;
			auto fGenes = fJ.second;
			auto rGenes = reverseJ[jName];
			bool matched = false;
			for (auto fGene : fGenes) {
				for (auto rGene : rGenes) {
					if (fGene != "-" && fGene == rGene)
						matched = true;
				}
			}
			if (matched) 
				junctionsMap[jName]->SetAsSpurious();
		}
	}

	void JunctionFilter::FilterSpuriousScrambledReads(AnchorReadMap *arMap, 
		Module *mod)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		string queryFile = c->JunctionFilter_dir + "scramble.query.fa";
		string targetFile = c->JunctionFilter_dir + "scramble.target.fa";

		if (!mod->IsComplete) {
			ofstream queryFS(queryFile.c_str(), ios::out);
			ofstream targetFS(targetFile.c_str(), ios::out);

			unordered_map<int, bool> scrambledGenes;
			for (auto arIter : (*arMap)){
				string anchorReadName = arIter.first;
				AnchorRead *ar = arIter.second;
				if (ar->isSpurious)
					continue;

				queryFS << ">" << anchorReadName << endl;
				queryFS << ar->splitRead.GetAlignedSeq() << endl;
				Gene *gene5p = ar->junction->ex5p->gene;
				Gene *gene3p = ar->junction->ex3p->gene;

				if (scrambledGenes.find(gene5p->geneId) == scrambledGenes.end()) {
					targetFS << ">" << gene5p->geneId << endl;
					targetFS << gene5p->GetScrambledExonSequence() << endl;
					scrambledGenes[gene5p->geneId] = true;
				}
				if (scrambledGenes.find(gene3p->geneId) == scrambledGenes.end()) {
					targetFS << ">" << gene3p->geneId << endl;
					targetFS << gene3p->GetScrambledExonSequence() << endl;
					scrambledGenes[gene3p->geneId] = true;
				}
			}
		}
		int numTargets = Utils::LineCount(targetFile)/2;
		ComputePerTask cpt = ComputePerTask::CalculateComputePerTask(6, 4, 2);
		//Paralellize this step to save time on index building
		vector<boost::thread *> threads;
	
		if (numTargets == 0)
			return;

		if (numTargets < 100)
			cpt.numSplits = 1;

		if (!mod->IsComplete){
			//Use SplitFileEvenly module to properly implement this;
			int linesPerSplit = ((numTargets + cpt.numSplits - 1) / (cpt.numSplits)) * 2;
			char cmd[1000];
			sprintf(cmd, "split -l %d -d -a 2 %s %s.",
				linesPerSplit, targetFile.c_str(), targetFile.c_str());
			if (Utils::ExecuteCommand(cmd).exit_code != 0) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") 
					<< "Error executing: " << cmd;
				exit(1);
			}
			for (int i = 0; i < cpt.numSplits; i++) {
				string targetSp = 
					targetFile + "." + (boost::format("%02d") % i).str();
				threads.push_back(
					new boost::thread(FilterSpuriousScrambledReads_worker,
					queryFile, targetSp, cpt.numCoresPerSplit));
			}

			for (int i = 0; i < threads.size(); i++) {
				threads[i]->join();
				delete threads[i];
			}
		}
		for (int i = 0; i < cpt.numSplits; i++){
			string bamFile = targetFile + "." +
				(boost::format("%02d") % i).str() + ".alignments.bam";
			BamTools::BamReader reader;
			BamTools::BamAlignment al;
			reader.Open(bamFile.c_str());
			BamTools::RefVector refVect = reader.GetReferenceData();

			while (reader.GetNextAlignment(al)){
				AnchorRead *ar = (*arMap)[al.Name];
				int targetGeneId = stoi(refVect[al.RefID].RefName);

				if (ar->junction->ex5p->gene->geneId == targetGeneId ||
					ar->junction->ex3p->gene->geneId == targetGeneId)
				{
					int targetStart = max<int>(0, ar->splitRead.overhang5p -
						ar->splitRead.GetMinOverhang()) + 1;
					int targetEnd = min<int>(al.Qualities.length(),
						ar->splitRead.overhang5p + ar->splitRead.GetMinOverhang()) - 1;

					int delBefore = 0, delAfter = 0;
					if (al.CigarData[0].Type == 'H' || al.CigarData[0].Type == 'S')
						delBefore = al.CigarData[0].Length;
					if (al.CigarData.size() >= 2 &&
						(al.CigarData[al.CigarData.size() - 1].Type == 'H' ||
						al.CigarData[al.CigarData.size() - 1].Type == 'S')
						)
					{
						delAfter = al.CigarData[al.CigarData.size() - 1].Length;
					}
					int editDistance = 0;
					al.GetTag("NM", editDistance);
					if (editDistance <= (al.QueryBases.length() * 0.05)
						&& delBefore <= targetStart
						&& (al.Qualities.length() - delAfter) >= targetEnd)
					{
						ar->isSpurious = true;
					}
				}
			}
		}
	}

	void JunctionFilter::FilterSpuriousScrambledReads_worker(string query, 
		string target, int cores)
	{
		Config *c = Config::GetConfig();
		char cmd[2000];
		sprintf(cmd, "%s %s %s > /dev/null ", c->bowtie2buildPath.c_str(), 
			target.c_str(), target.c_str());
		Utils::ExecuteCommand(cmd);

		sprintf(cmd, "%s --sensitive-local -a -p %d --no-unal --quiet -x %s "
			"-f %s | %s view -S -b - > %s.alignments.bam 2> /dev/null  ",
			c->bowtie2Path.c_str(), cores, target.c_str(), query.c_str(),
			c->samtoolsPath.c_str(), target.c_str());
		Utils::ExecuteCommand(cmd);
	}

	void JunctionFilter::FilterJunctionPCRdups(AnchorReadMap *arMap)
	{
		unordered_map<string, bool> junctionsProcessed;
		for (auto arIter : (*arMap)){
			Junction *j = arIter.second->junction;
			string junctName = j->GetJunctionName();
			if (junctionsProcessed.find(junctName) != junctionsProcessed.end())
				continue;
			junctionsProcessed[junctName] = true;
			j->MarkPCRduplicates();
		}
	}

	string JunctionFilter::GetBlatParamsForStringency(FilterStringency f)
	{
		switch (f)
		{
			case LOW:
				return "-t=dna -q=dna -noHead -oneOff=1";
				break;
			case MEDIUM:
				return "-t=dna -q=dna -noHead -oneOff=1";
				break;
			case HIGH:
				return "-t=dna -q=dna -noHead -oneOff=1 stepSize=5";
				break;
		}
		return "";
	}

	void JunctionFilter::CleanUp()
	{
		Config *c = Config::GetConfig();
		Utils::DeleteFile(c->JunctionFilter_dir + "step*_genome.fa.*");
		Utils::DeleteFile(c->JunctionFilter_dir + "step*_exome.fa.*");
		Utils::DeleteFile(c->JunctionFilter_dir + "scramble.*");
		Utils::DeleteFile(c->JunctionFilter_dir + "fusions.junctions.chr*");
	}
}