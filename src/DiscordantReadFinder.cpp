
#include "DiscordantReadFinder.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >, 
		(keywords::channel = ""));

	string DiscordantReadFinder::DChannel = "DiscordantReadFinder";
	int DiscordantReadFinder::TotalUnalignedReads = 0;
	int DiscordantReadFinder::TotalDiscordantReads = 0;
	boost::mutex DiscordantReadFinder::statsUpdateMutex1;
	boost::mutex DiscordantReadFinder::statsUpdateMutex2;

	//
	// DiscordantReadFinder::Run is the main function in this singleton class. It 
	// first aligns all fastq reads to all possible isoforms in a given transcriptome.  
	// Unaligned reads are then trimmed and aligned to the unspliced transcriptome to
	// identify discordant reads.
	//
	bool DiscordantReadFinder::Run() 
	{
		Config *c = Config::GetConfig(); 
		GeneModel *gm = GeneModel::GetGeneModel();
		
		Logger::RegisterChannel(c->sampleOutputLogDir +
			"DiscordantReadFinder.log", DChannel);

		Module modExtractUnaligned = 
			Module::LoadStatusFile(c->statusFile_ExtractUnaligned);

		//Given available compute, determines how many steps to parallelize the
		//aligment steps
		ComputePerTask cpt = 
			ComputePerTask::CalculateComputePerTask(6, 2, c->maxBwaMem);
		BOOST_LOG_CHANNEL(logger::get(), "Main")
			<< "Extracting unaligned reads in " << cpt.numSplits << " splits";

		if (!modExtractUnaligned.IsComplete) {
			ExtractUnalignedReads(cpt);
			modExtractUnaligned.NumSplits = cpt.numSplits;
			modExtractUnaligned.IsComplete = true;
			modExtractUnaligned.SaveStatusFile();
		}
		else {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "\tSKIPPING: detected output from previous run";
		}
		ifstream rcFile(c->readcountFile.c_str());
		string s;
		getline(rcFile, s);
		BOOST_LOG_CHANNEL(logger::get(), "Main") << s;
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Total # of unaligned reads: "
			<< TotalUnalignedReads;
		//Parse out the read count;
		int readcount = c->GetTotalReadcount();

		try{
			auto sp = Utils::SplitToVector(c->minSpanFunct, ",");
			int R = stoi(sp[0]), X = stoi(sp[1]), Y = stoi(sp[2]);
			if (Y == 0) {
				c->minSpanCount = stoi(sp[0]);
			}
			else {
				c->minSpanCount =
					std::ceil(R + X * max<double>(0, std::log(((double)readcount) / Y)));
			}
		}
		catch (std::exception e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error parsing --min_span = "
				<< c->minSpanFunct << ". Exiting with error: " << e.what();
			exit(1);
		}

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Finding gene pairs supported by "
			<< c->minSpanCount << " (--min_span) or more discordant reads";

		for (int i = 0; i < modExtractUnaligned.NumSplits; i++) {
			string fq1 = c->workingDir + "fastqs/split_" + lexical_cast<string>(i)
				+"_unaligned_1.fastq";
			string fq2 = c->workingDir + "fastqs/split_" + lexical_cast<string>(i)
				+"_unaligned_2.fastq";
			if (!Utils::FileExists(fq1) || !Utils::FileExists(fq2)) 
				continue;
			c->ExtractUnaligned_fq1.push_back(fq1);
			c->ExtractUnaligned_fq2.push_back(fq2);
		}

		if (c->ExtractUnaligned_fq1.size() != modExtractUnaligned.NumSplits
			|| c->ExtractUnaligned_fq1.size() != c->ExtractUnaligned_fq2.size()) {
			if (modExtractUnaligned.IsComplete)
				BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "Error occurred while trying to resume previously started run. "
				<< "Unaligned reads are not properly accounted for.  Please delete "
				<< "output directory and re-run MOJO";
			else
				BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "Error occurred while extracting aligned reads. Please check "
				<< "input fastqs";
			exit(1);
		}

		Module modDiscordantFinder =
			Module::LoadStatusFile(c->statusFile_DiscordantFinder);

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Identifying discordant reads"
			<< " in " << cpt.numSplits << " splits";
		if (!modDiscordantFinder.IsComplete) {
			FindDiscordantReads(cpt);
			modDiscordantFinder.NumSplits = cpt.numSplits;
			modDiscordantFinder.IsComplete = true;
		}
		else 
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\tSKIPPING: detected output from previous run";

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Total # of discordant reads: "
			<< TotalDiscordantReads;

		if (!Utils::FileExists(c->workingDir + "discordant.flat")) {
			if (modDiscordantFinder.IsComplete)
				BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "Error occurred while trying to resume previously started "
				<< "run. DicordantReadFinder did not finish successfully.  "
				<< "Please delete output directory and re-run MOJO";
			else
				BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "Error occurred while extracting discordant reads ";
			exit(1);
		}

		if (Utils::LineCount(c->workingDir + "discordant.flat") == 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "No discordant reads found. Exiting.";
			exit(1);
		}
		modDiscordantFinder.SaveStatusFile();
		return true;;
	}

	void DiscordantReadFinder::ExtractUnalignedReads(ComputePerTask cpt) 
	{
		Config *c = Config::GetConfig();
		vector<boost::thread *> threads;
		boost::filesystem::create_directories(c->workingDir + "/fastqs/");
		try {
			int residualCores = cpt.numResidualCores;
			for (int threadId = 0; threadId < cpt.numSplits; threadId++) {
				string threadIdStr = lexical_cast<string>(threadId);
				string chan = "ExtractUnaligned." + threadIdStr;
				Logger::RegisterChannel(c->sampleOutputLogDir + chan + ".log", chan);
				threads.push_back(new boost::thread(ExtractUnalignedReads_worker,
					threadId, cpt, chan));
					//threadIdStr, cpt.numCoresPerSplit + (residualCores-- > 0 ? 1 : 0), chan));
			}
			for (vector<boost::thread *>::size_type j = 0; j < threads.size(); j++) {
				threads[j]->join();
				delete threads[j];
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occurred while "
				<< "extracting Unaligned Reads. Error: " << e.what();
		}
	}

	void DiscordantReadFinder::ExtractUnalignedReads_worker(int thread, 
		ComputePerTask cpt, string chan)
	{
		Config *c = Config::GetConfig();
		try
		{
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tStarted split # "
				<< thread << ". Log: ./logs/" << (chan + ".log[.cmds]");
			string fq1ss = Config::GenerateCatCmdsForSplit(1, thread, cpt.numSplits);
			string fq2ss = Config::GenerateCatCmdsForSplit(2, thread, cpt.numSplits);
			
			string pfx = c->workingDir + "/fastqs/split_" + 
				lexical_cast<string>(thread);
			//Align to an index of all canonical isoforms; filter out all reads concordantly 
			//aligning to the transcriptome
			char bt2buf[15000];
			sprintf(bt2buf, "%s\n %s\n %s -p %d -x %s -1 %s_1.fastq -2 %s_2.fastq "
				"--dovetail --no-mixed --no-discordant --un-conc "
				"%s_unaligned_%%.fastq.tmp --score-min L,-2,-0.2 > /dev/null", 
				fq1ss.c_str(), fq2ss.c_str(), c->bowtie2Path.c_str(), 
				cpt.numCoresPerSplit, c->bowtie2AllIsoformIndex.c_str(), 
				pfx.c_str(), pfx.c_str(), pfx.c_str());

			if (Utils::ExecuteCommand(bt2buf, chan).exit_code != 0)
				exit(1);
			if (c->removeTemporaryFiles) {
				Utils::DeleteFiles(std::vector < string > {
					pfx + "_1.fastq", pfx + "_2.fastq"});
			}
			BOOST_LOG_CHANNEL(logger::get(), chan)
				<< "First pass alignment to AllIsoformIndex completed";

			//Perform a more rigorous search to remove paired-end reads with each  
			//end aligning in close proximity (within fragment size) to each other
			sprintf(bt2buf, "%s -p %d -k 4 -x %s -U %s_unaligned_1.fastq.tmp "
				" | cut -f 1,3 | grep -v '*' | grep -P -v \"^@\" "
				"> %s_unaligned_1.align2gene",
				c->bowtie2Path.c_str(), cpt.numCoresPerSplit, 
				c->bowtie2AllIsoformIndex.c_str(), pfx.c_str(), pfx.c_str());

			if (Utils::ExecuteCommand(bt2buf, chan).exit_code != 0) exit(1);
			sprintf(bt2buf, "%s -p %d -k 4 -x %s -U %s_unaligned_2.fastq.tmp "
				" | cut -f 1,3 | grep -v '*' | grep -P -v \"^@\" "
				" > %s_unaligned_2.align2gene",
				c->bowtie2Path.c_str(), cpt.numCoresPerSplit,
				c->bowtie2AllIsoformIndex.c_str(), pfx.c_str(), pfx.c_str());
			if (Utils::ExecuteCommand(bt2buf, chan).exit_code != 0) exit(1);
			BOOST_LOG_CHANNEL(logger::get(), chan) 
				<< "Second pass alignment to AllIsoformIndex completed";

			//scoping this with the expectation that e1Genes and e2Genes 
			//will be released after exiting
			unordered_map<int, bool> concordantReads;
			{
				unordered_map<int, unordered_map<int, bool> > e1Genes;
				ifstream end1ReadToGeneFile((pfx + "_unaligned_1.align2gene").c_str());
				ifstream end2ReadToGeneFile((pfx + "_unaligned_2.align2gene").c_str());
				for (std::string str; getline(end1ReadToGeneFile, str);) {
					vector<string> sp = Utils::SplitToVector(str, "\t");
					int readId = stoi(Read::TrimReadName(sp[0]));
					int geneId = stoi(sp[1].substr(1, sp[1].length() - 1));
					e1Genes[readId][geneId] = 1;
				}
				for (std::string str; getline(end2ReadToGeneFile, str);){
					vector<string> sp = Utils::SplitToVector(str, "\t");
					int readId = stoi(Read::TrimReadName(sp[0]));
					int geneId = stoi(sp[1].substr(1, sp[1].length() - 1));
					if (e1Genes.find(readId) != e1Genes.end())
						if (e1Genes[readId].find(geneId) != e1Genes[readId].end())
							concordantReads[readId] = true;
				}
				end1ReadToGeneFile.close();
				end2ReadToGeneFile.close();
			}

			// Process <>_split_unaligned_1.aligned_gene; both ends of .aligned_gene  
			// are sorted by read name.  Perform intersection of sorted lists.  
			// Previous version implemented this using hashes
			/*	
			vector<string> end1AlignsToGene, end2AlignsToGene, bothAlignsToGene;

			ifstream end1ReadToGeneFile((pfx + "_unaligned_1.align2gene").c_str());
			ifstream end2ReadToGeneFile((pfx + "_unaligned_2.align2gene").c_str());
			for (std::string str; getline(end1ReadToGeneFile, str);) {
				vector<string> sp = Utils::SplitToVector(str, "\t");
				end1AlignsToGene.push_back(Read::TrimReadName(sp[0]) + "_" + sp[1]);
			}
			for (std::string str; getline(end2ReadToGeneFile, str);){
				vector<string> sp = Utils::SplitToVector(str, "\t");
				end2AlignsToGene.push_back(Read::TrimReadName(sp[0]) + "_" + sp[1]);
			}
			std::set_intersection(end1AlignsToGene.begin(), end1AlignsToGene.end(),
				end2AlignsToGene.begin(), end2AlignsToGene.end(),
				std::inserter(bothAlignsToGene, bothAlignsToGene.begin()));

			unordered_map<string, bool> concordantReads;
			for (auto readName : bothAlignsToGene) {
				concordantReads[Read::TrimReadName(
					readName.substr(0, readName.find("||", 0)))] = true;
			}
			*/

			ofstream end1FqOut((pfx + "_unaligned_1.fastq").c_str(), ios::out);
			ofstream end2FqOut((pfx + "_unaligned_2.fastq").c_str(), ios::out);
			//ofstream fqi((pfx + "_unaligned.fastqIndex").c_str(), ios::out);

			FastqParser fpTmp(pfx + "_unaligned_1.fastq.tmp", 
				pfx + "_unaligned_2.fastq.tmp");
			long readId = 1;
			PairedRead pr;
			while (fpTmp.GetNextPairedRead(&pr)) {
				int readId = stoi(pr.FirstRead.GetTrimmedReadName());
				if (concordantReads.find(readId) != concordantReads.end())
					continue;

				//long long end1Pos = end1FqOut.tellp();
				//long long end2Pos = end2FqOut.tellp();
				//fqi << pr.FirstRead.GetTrimmedReadName() << "\t" << readId << "\t"
				//	<< end1Pos << "\t" << end2Pos << endl;
				//pr.FirstRead.ReadName = lexical_cast<string>(readId)+"/1";
				//pr.SecondRead.ReadName = lexical_cast<string>(readId)+"/2";

				char end1buf[500], end2buf[500];
				pr.FirstRead.GetFastqString(end1buf);
				pr.SecondRead.GetFastqString(end2buf);

				end1FqOut.write(end1buf, strlen(end1buf));
				end2FqOut.write(end2buf, strlen(end2buf));
				readId++;
			}
			end1FqOut.close();
			end2FqOut.close();

			//Finally, curate the .readID to only include unaligned reads;
			FastqParser *fp = new FastqParser(pfx + "_unaligned_1.fastq",
				pfx + "_unaligned_2.fastq");
			
			unordered_map<int, bool> unalReadIdMap;
			while (fp->GetNextPairedRead(&pr)) 
				unalReadIdMap[stoi(pr.GetTrimmedReadName())] = true;
			
			ifstream allReadIdFile((pfx + "_1.fastq.readID").c_str(), ios::out); 
			ofstream unalReadIdFile((pfx + "_unaligned.readID").c_str(), ios::out);
			for (string str; getline(allReadIdFile, str);) {
				auto sp = Utils::SplitToVector(str, "\t");
				if (unalReadIdMap.find(stoi(sp[0])) != unalReadIdMap.end())
					unalReadIdFile << str << endl;
			}
			
			if (c->removeTemporaryFiles) {
				Utils::DeleteFiles(std::vector < string > {
					pfx + "_unaligned_1.fastq.tmp", pfx + "_unaligned_2.fastq.tmp"});
				Utils::DeleteFiles(std::vector < string > {pfx + "_unaligned_1.align2gene",
					pfx + "_unaligned_2.align2gene"});
				Utils::DeleteFile(pfx + "_1.fastq.readID");
				Utils::DeleteFiles(std::vector < string > {pfx + "_1.fastq",
					pfx + "_2.fastq"});
			}
			int readCount = (Utils::LineCount(pfx + "_unaligned_1.fastq")) / 4;
			UpdateUnalignedReadCount(readCount);
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tFinished split # " << thread ;

		}
		catch (std::exception &e)
		{
			std::cerr << "Exception occurred in ExtractUnalignedReads-Thread: " << thread
				<< ". Terminating run." << endl << "Exception: " << e.what() << endl;
			exit(1);
		}
	}

	void DiscordantReadFinder::FindDiscordantReads(ComputePerTask cpt)
	{
		Config *c = Config::GetConfig();

		vector<boost::thread *> threads;
		int residualCores = cpt.numResidualCores;
		string alignmentFiles;
		try {
			for (int threadId = 0; threadId < cpt.numSplits; threadId++) {
				string threadIdStr = lexical_cast<string>(threadId);
				string chan = "DiscordantReads." + threadIdStr;
				Logger::RegisterChannel(c->sampleOutputLogDir + chan + ".log", chan);
				threads.push_back(new boost::thread(FindDiscordantReads_worker,
					threadIdStr, cpt.numCoresPerSplit + 
					(residualCores-- > 0 ? 1 : 0), chan));
			}
			for ( auto thread : threads ) {
				thread->join();
				delete thread;
			}

			//Merge all flat files;
			stringstream catSS, rmSS;
			catSS << "cat ";
			for (int i = 0; i < cpt.numSplits; i++) {
				catSS << c->workingDir << "split_" << i << "/split_"
					<< i << ".discordant.flat ";
				rmSS << "rm " << c->workingDir << "split_" << i << "/split_"
					<< i << ".discordant.flat " << endl;
			}
			catSS << " > " << GetFlatFilename();
			if (Utils::ExecuteCommand(catSS.str().c_str()).exit_code != 0) {
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "Error compiling discordant.flat. Exiting.";
				exit(1);
			}
			Utils::ExecuteCommand(rmSS.str().c_str());
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occurred while "
				<< "finding discordant reads. Error: " << e.what();
		}

		////merge alignments
		//char cmd[10000];
		//sprintf(cmd, "%s merge -f -@ %d %s/splits.alignments.bam %s",
		//	c->samtoolsPath.c_str(), min<int>(c->maxCores, 12),
		//	c->workingDir.c_str(), alignmentFiles.c_str());
		//if (Utils::ExecuteCommand(cmd).exit_code != 0) {
		//	BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occured while"
		//		<< " merging alignments files" << endl;
		//	exit(1);
		//}

		////sort merged alignments
		//sprintf(cmd, "%s sort -f -@ %d %s/splits.alignments.bam %s/splits.alignments.sorted.bam",
		//	c->samtoolsPath.c_str(), min<int>(c->maxCores, 12),
		//	c->workingDir.c_str(), c->workingDir.c_str());
		//if (Utils::ExecuteCommand(cmd).exit_code != 0) {
		//	BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occured while"
		//		<< " sorting alignments files" << endl;
		//	exit(1);
		//}

		////index merged alignments
		//sprintf(cmd, "%s index %s/splits.alignments.sorted.bam",
		//	c->samtoolsPath.c_str(), c->workingDir.c_str());
		//if (Utils::ExecuteCommand(cmd).exit_code != 0) {
		//	BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occured while"
		//		<< " indexing sorted alignments file" << endl;
		//	exit(1);
		//}
	}

	// Finds discordant reads in two steps. In Step 1, unaligned reads are trimmed
	// to 36bps and aligned to the transcriptome.  In step 2, unaligned reads in
	// step 1 are trimmed at the 5' and 3' end to skip a potential splice junction
	void DiscordantReadFinder::FindDiscordantReads_worker(string thread, int tCores, 
		string chan)
	{
		Config *c = Config::GetConfig();

		boost::filesystem::create_directories(c->workingDir + "/split_" + thread);
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tStarted split # " << thread
			<< ". Log: ./logs/" << (chan + ".log[.cmds]");
		// Trim to 36bps and peform alignments in two steps;
		string unalFq1 = c->workingDir + "/fastqs/split_" + thread + 
			"_unaligned_1.fastq";
		string unalFq2 = c->workingDir + "/fastqs/split_" + thread + 
			"_unaligned_2.fastq";
		string unalFqIndex = c->workingDir + "/fastqs/split_" + thread + 
			"_unaligned.fastqIndex";
		//
		// Iteration 1
		//
		string splitPrefix = c->workingDir + "/split_" + thread + 
			"/split_" + thread + "_iteration1";
		string iterFq1 = splitPrefix + "_1.fastq", iterFq2 = splitPrefix + "_2.fastq";
		FastqParser p(unalFq1, unalFq2);
		p.CreateTrimmedFiles(iterFq1, iterFq2, 36);
		BOOST_LOG_CHANNEL(logger::get(), DChannel) 
			<< "Created trimmed fastqs for iteration 1 ";
		BOOST_LOG_CHANNEL(logger::get(), DChannel) 
			<< "Starting alignments for iteration 1 ";
		
		DiscordantReadFinder::FindDiscordantReads_ByIteration(splitPrefix, thread, 
			tCores, chan);
		
		BOOST_LOG_CHANNEL(logger::get(), DChannel) 
			<< "Completed alignments for iteration 1 ";
		if (c->removeTemporaryFiles) {
			Utils::DeleteFiles(std::vector < string > { 
				splitPrefix + "_1.fastq", splitPrefix + "_2.fastq"});
		}

		// Extract unaligned reads from Iteration 1
		BamReader reader;
		string iter1Bam = splitPrefix + ".sorted.bam";
		if (!reader.Open(iter1Bam)) {
			BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "Failed to open " << iter1Bam;
			exit(1);
		}

		BamAlignment aln;

		// unmappedReads --> <readname> : <which_end> key pair; <which_end> denotes, 
		//the end that is unaligned and therefore, needs to be slided
		unordered_map<int, int> unmappedReads;	
		unordered_map<int, bool> ignoreReads;

		while (reader.GetNextAlignment(aln)){
			string readName = Read::TrimReadName(aln.Name);
			int readId = stoi(readName);
			if (aln.IsMapped() && aln.IsMateMapped()) {
				//Some alignments have flag 7 set (doesn't make sense) -- exclude those
				if (aln.IsPaired() && !aln.IsMapped() && aln.IsProperPair()) 
					ignoreReads[readId] = true;
				unmappedReads[readId] = 0;	// 0: both ends are aligned;
			}
			else 
				unmappedReads[readId] = aln.IsFirstMate() ? 2 : 1;
		}
		reader.Close();
		FastqParser t(unalFq1, unalFq2);

		//Iteration2 parameters;
		splitPrefix = c->workingDir + "/split_" + thread + "/split_" + 
			thread + "_iteration2";
		iterFq1 = splitPrefix + "_1.fastq", iterFq2 = splitPrefix + "_2.fastq";
		ofstream end1FqOut(iterFq1.c_str(), ios::out);
		ofstream end2FqOut(iterFq2.c_str(), ios::out);
		int iteration1_unaligned = 0;

		PairedRead pr;
		while (t.GetNextPairedRead(&pr)) {
			int readId = stoi(pr.GetTrimmedReadName());
			if (ignoreReads.find(readId) != ignoreReads.end())
				continue;
			auto result = unmappedReads.find(readId);
			char buf1[2000], buf2[2000];
			if (result == unmappedReads.end()) {
				pr.FirstRead.GetFastqSubstring(buf1, 36, 36);
				pr.SecondRead.GetFastqSubstring(buf2, 36, 36);
			}
			else if ((*result).second == 1) {
				pr.FirstRead.GetFastqSubstring(buf1, 36, 36);
				pr.SecondRead.GetFastqSubstring(buf2, 36, 0);
			}
			else if ((*result).second == 2){
				pr.FirstRead.GetFastqSubstring(buf1, 36, 0);
				pr.SecondRead.GetFastqSubstring(buf2, 36, 36);
			}
			else
				continue;

			iteration1_unaligned++;
			end1FqOut.write(buf1, strlen(buf1));
			end2FqOut.write(buf2, strlen(buf2));
		}
		end1FqOut.close(), end2FqOut.close();

		// Iteration 2
		if (iteration1_unaligned > 100)  {
			BOOST_LOG_CHANNEL(logger::get(), chan) 
				<< "Starting alignments for iteration 2 ";
			FindDiscordantReads_ByIteration(splitPrefix, thread, tCores, chan);
			BOOST_LOG_CHANNEL(logger::get(), chan) 
				<< "Completed alignments for iteration 2 ";
		}
		else {
			BOOST_LOG_CHANNEL(logger::get(), chan) 
				<< "Skipping iteration 2 (too few unaligned)";
		}

		// Generate discordant.bam and discordant.flat
		BamAlignment feAl, seAl;
		BamWriter discWriter, alnWriter;
		splitPrefix = c->workingDir + "split_" + thread + "/" + "split_" + thread;
		discWriter.Open((splitPrefix + ".discordant.bam").c_str(), 
			reader.GetHeaderText(), reader.GetReferenceData());
		//alnWriter.Open((splitPrefix + ".alignments.bam").c_str(),
		//	reader.GetHeaderText(), reader.GetReferenceData());
		ofstream discFlat((splitPrefix + ".discordant.flat").c_str(), ios::out);
		FastqParser *unalFp = new FastqParser(unalFq1, unalFq2, unalFqIndex);
		
		for (auto iteration : { "1", "2" }){
			string pairedBam = splitPrefix + "_iteration" + iteration + 
				".paired.sorted.bam";
			if (!Utils::FileExists(pairedBam)) continue;

			BamReader iterAlnReader;
			iterAlnReader.Open(pairedBam.c_str());
			RefVector headerVect = iterAlnReader.GetReferenceData();
			unordered_map<string, int> refNameToId;
			for (int rv = 0; rv < headerVect.size(); rv++)
				refNameToId[headerVect[rv].RefName] = rv;

			bool skip = false;
			while (1) {
				if (!skip && !iterAlnReader.GetNextAlignment(feAl)) break;
				skip = false;
				if (!iterAlnReader.GetNextAlignment(seAl)) break;

				if (feAl.Name != seAl.Name) {
					skip = true;
					feAl = seAl;
					continue;
				}
				string readname = Read::TrimReadName(feAl.Name);
				//int readId = stoi(readname);
				//PairedRead pr = unalFp->GetReadById(readId);
				//string tagValue = thread + "," + readname + "," 
				//	+ pr.FirstRead.Sequence + "," + pr.SecondRead.Sequence ;
				//feAl.AddTag("XJ", "Z", tagValue);
				//seAl.AddTag("XJ", "Z", tagValue);

				vector<BamAlignment> fAligns = GetBwaAlignsFromXA(feAl, &refNameToId);
				vector<BamAlignment> sAligns = GetBwaAlignsFromXA(seAl, &refNameToId);

				if (feAl.RefID == seAl.RefID && 
					(fAligns.size() > 2 || sAligns.size() > 2))
					continue;

				string uniqueFlag = "0";
				if (fAligns.size() == 1 && sAligns.size() == 1)
					uniqueFlag = "1";
				
				int readcount = 0;
				unordered_map<string, bool> genesWritten;  //avoid duplicating the XA tag pairs
				typedef vector<BamAlignment>::size_type ba_type;
				for (ba_type f = 0; f < fAligns.size(); f++) {
					for (ba_type s = 0; s < sAligns.size(); s++) {
						BamAlignment fe = fAligns[f], se = sAligns[s];
						if (fe.RefID == se.RefID) {
							if (f == 0 && s == 0) {
								discWriter.SaveAlignment(fe);
								discWriter.SaveAlignment(se);
							}
						}
						else {
							if (readcount > 0) {
								fAligns[f].Name = readname + "_" + 
									lexical_cast<string>(readcount);
								sAligns[s].Name = readname + "_" + 
									lexical_cast<string>(readcount);
							}
							stringstream key, revKey;
							key << fAligns[f].RefID << "_" << sAligns[s].RefID;
							revKey << sAligns[s].RefID << "_" << fAligns[f].RefID;

							if (genesWritten.find(key.str()) != genesWritten.end() || 
								genesWritten.find(revKey.str()) != genesWritten.end())
								continue;
							genesWritten[key.str()] = "";
							genesWritten[revKey.str()] = "";
							discWriter.SaveAlignment(fAligns[f]);
							discWriter.SaveAlignment(sAligns[s]);

							discFlat << BuildFlatReadString(fAligns[f], sAligns[s], 
								&headerVect);
							discFlat << "\t" << uniqueFlag << endl;
							readcount++;
							if ((fAligns.size() > 3 && sAligns.size() > 3) 
								|| fAligns.size() > 10 || sAligns.size() > 10)
								s = sAligns.size(), f = fAligns.size();
						}
						//alnWriter.SaveAlignment(feAl);
						//alnWriter.SaveAlignment(seAl);
					}
				}
			}
		}
		discWriter.Close();
		alnWriter.Close();
		discFlat.close();

		if (c->removeTemporaryFiles) {
			Utils::DeleteFile(splitPrefix + "_iteration*fastq");
			Utils::DeleteFile(splitPrefix + "_iteration*sorted.bam");
			Utils::DeleteFile(splitPrefix + ".discordant.bam");
		}

		int readCount = (Utils::LineCount(splitPrefix + ".discordant.flat")) / 4;
		UpdateDiscordantReadCount(readCount);
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tFinished split # " << thread ;
	}

	void DiscordantReadFinder::FindDiscordantReads_ByIteration(string prefix, 
		string thread, int tCores, string chan)
	{
		Config *c = Config::GetConfig();

		string end1Fq = prefix + "_1.fastq", end2Fq = prefix + "_2.fastq";
		BOOST_LOG_CHANNEL(logger::get(), DChannel) << "Starting thread " << thread
			<< " of FindDiscordantReads process. Log: " 
			<< (c->sampleOutputLogDir + chan + ".log");

		char cmd[5000];
		sprintf(cmd, "%s aln -q 15 -t %d %s %s > %s_aln_1.sai 2> %s_aln_1_output.log", 
			c->bwaPath.c_str(), tCores, c->bwaTranscriptomeIndex.c_str(), 
			end1Fq.c_str(), prefix.c_str(), prefix.c_str());
		if (Utils::ExecuteCommand(cmd, chan).exit_code != 0) exit(1);
		BOOST_LOG_CHANNEL(logger::get(), DChannel) << "Completed aligning first end ";

		sprintf(cmd, "%s aln -q 15 -t %d %s %s > %s_aln_2.sai 2> %s_aln_2_output.log",
			c->bwaPath.c_str(), tCores, c->bwaTranscriptomeIndex.c_str(), end2Fq.c_str(),
			prefix.c_str(), prefix.c_str());
		if (Utils::ExecuteCommand(cmd, chan).exit_code != 0) exit(1);
		BOOST_LOG_CHANNEL(logger::get(), DChannel) << "Completed aligning second end ";

		sprintf(cmd, "%s sampe -A -a 1000 -N 25 -c 0.0001 -P %s %s_aln_1.sai %s_aln_2.sai"
				" %s %s 2> %s_sampe_output.log | %s view -b -S - -o %s.bam",
				c->bwaPath.c_str(), c->bwaTranscriptomeIndex.c_str(),  prefix.c_str(), 
				prefix.c_str(), end1Fq.c_str(), end2Fq.c_str(), prefix.c_str(), 
				c->samtoolsPath.c_str(), prefix.c_str());
		if (Utils::ExecuteCommand(cmd, chan).exit_code != 0) exit(1);
		BOOST_LOG_CHANNEL(logger::get(), DChannel) << "Completed 'bwa sampe' step.";

		ComputePerTask cpt = ComputePerTask::CalculateComputePerTask(c->maxCores, 1, 1);
		sprintf(cmd, "%s view -F 4 -b %s.bam | %s sort -n -@ %d - %s.sorted",
			c->samtoolsPath.c_str(), prefix.c_str(), c->samtoolsPath.c_str(), 
			cpt.numSplits, prefix.c_str());
		if (Utils::ExecuteCommand(cmd, chan).exit_code != 0) exit(1);

		if (c->removeTemporaryFiles)
			Utils::DeleteFiles(std::vector < string > { prefix + "*sai", 
				prefix + ".bam", prefix + ".bam" });

		BamReader reader;
		BamWriter writer;
		BamAlignment aln;

		reader.Open((prefix + ".sorted.bam").c_str());
		writer.Open((prefix + ".paired.sorted.bam").c_str(), 
			reader.GetHeaderText(), reader.GetReferenceData());
		while (reader.GetNextAlignment(aln)) 
			//For some reason, some alignments have flag 7 set (which doesn't make sense) -- exclude those
			if (!(aln.IsPaired() && !aln.IsMapped() && aln.IsProperPair())) 
				writer.SaveAlignment(aln);
		reader.Close(), writer.Close();
	}

	void DiscordantReadFinder::UpdateUnalignedReadCount(int cnt) 
	{
		boost::mutex::scoped_lock lck(statsUpdateMutex1);
		TotalUnalignedReads += cnt;
	}

	void DiscordantReadFinder::UpdateDiscordantReadCount(int cnt)
	{
		boost::mutex::scoped_lock lck(statsUpdateMutex2);
		TotalDiscordantReads += cnt;
	}


	//BWA stores alternate alignments in the XA tag, this function retrieve them
	vector<BamAlignment> DiscordantReadFinder::GetBwaAlignsFromXA(
		BamAlignment al, unordered_map<string, int> *geneToNameMap) 
	{
		struct XASplit 
		{
			public:
				string gene, cigar;
				int position, nm;
				XASplit(string alignment) {
					vector<string> sp = Utils::SplitToVector(alignment, ",");
					gene = sp[0], cigar = sp[2];
					// offset by -1 to convert to 0-based offset for bam
					position = abs(stoi(sp[1])) - 1;
					nm = stoi(sp[3]);
			}
		};

		vector<BamAlignment> aligns;
		aligns.push_back(al);
		string xa = "";
		int nm = 0;
		al.GetTag("NM", nm);
		if (al.GetTag("XA", xa)) {
			vector<string> splits = Utils::SplitToVector(xa, ";");
			for (vector<string>::size_type i = 0; i < splits.size(); i++) {
				BamAlignment tmpAl = al;
				if (splits[i].length() < 2) continue;
				XASplit s(splits[i]);
				if (nm < s.nm) continue;
				tmpAl.RefID = (*geneToNameMap)[s.gene];
				tmpAl.Position = s.position;
				aligns.push_back(tmpAl);
			}
		}
		return aligns;
	}

	// A .flat file contains all discordant reads identified in the previous steps
	string DiscordantReadFinder::BuildFlatReadString(BamAlignment fAl, 
		BamAlignment sAl, RefVector *refVect) 
	{
		string gA = (*refVect)[fAl.RefID].RefName;
		string gB = (*refVect)[sAl.RefID].RefName;
		string posA = lexical_cast<string>(fAl.Position);
		string posB = lexical_cast<string>(sAl.Position);
		int lenA = 0, lenB = 0;
		for (int i = 0; i < fAl.CigarData.size(); i++)
			if (fAl.CigarData[i].Type == 'M')
				lenA += fAl.CigarData[i].Length;

		for (int i = 0; i < sAl.CigarData.size(); i++)
			if (sAl.CigarData[i].Type == 'M')
				lenB += sAl.CigarData[i].Length;
		string lenAstr = lexical_cast<string>(lenA);
		string lenBstr = lexical_cast<string>(lenB);
		string orientA = "+", orientB = "+";

		if (fAl.IsReverseStrand()) orientA = "-";
		if (sAl.IsReverseStrand()) orientB = "-";

		string flatA = gA + "\t" + posA + "\t" + lenAstr + "\t" + orientA + "\t" + 
			lexical_cast<string> (lexical_cast<int> (fAl.IsSecondMate()) + 1);
		string flatB = gB + "\t" + posB + "\t" + lenBstr + "\t" + orientB + "\t" + 
			lexical_cast<string> (lexical_cast<int> (sAl.IsSecondMate()) + 1);

		if (gA.compare(gB) < 0)
			return flatA + "\t" + flatB + "\t" + fAl.Name;
		else
			return flatB + "\t" + flatA + "\t" + fAl.Name;
	}

	string DiscordantReadFinder::GetFlatFilename()
	{
		Config *c = Config::GetConfig();
		return c->workingDir + "discordant.flat";
	}
};
