
#include "JunctionAligner.h"

#define MAX_JUNCTIONS_PER_SPLIT 2000000

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >, 
		(keywords::channel = ""));
	
	string JunctionAligner::channel = "JunctAlign";

	//
	// JunctionAligner::Run is the main function of this singleton class.  Takes
	// discordatn clusters as input and constructs exon-exon junctions between
	// candidate fusion junctions.  Then, all unaligned reads are mapped to this
	// exon-exon junction library to identify anchor reads
	//
	void JunctionAligner::Run(vector<DiscordantCluster *> clusters)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();
		Module modJunctionAligner = Module::LoadStatusFile(
			c->statusFile_JunctionAligner);

		if (!modJunctionAligner.IsComplete) {
			ComputePerTask cpt = ComputePerTask::CalculateComputePerTask(6, 2, 6);
			modJunctionAligner.NumSplits = cpt.numSplits;	//
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Building candidate fusion junctions for each cluster";
			int junctCount = BuildJunctions(clusters);
			if (junctCount == 0) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") 
					<< "No valid candidate junctions identified. Exiting.";
				exit(0);
			}
			BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "\tConstructed " << junctCount << " junctions";
			AlignToJunctions(&modJunctionAligner);
			GenerateFastqForJunctAlignments(&modJunctionAligner);
			modJunctionAligner.IsComplete = true;
			modJunctionAligner.SaveStatusFile();
		}
		else {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Found junctions from previous MOJO run";
			int junctCount = LoadJunctionsIntoClusters(clusters);
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "\tLoaded " << junctCount << " junctions";
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "Identified " << LoadAnchorReads(clusters) << " anchor reads";
		
		if (c->removeTemporaryFiles) {
			Utils::DeleteFile(c->workingDir + "/junct_aligns/junctions.fa.*");
		}
	}

	int JunctionAligner::BuildJunctions(vector<DiscordantCluster *> clusters) 
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		int neighbor_cnt = 3, junctionCount = 0;
		string junctsFileName = c->workingDir + "/junctions.fa";

		ofstream junctionsFileStream(junctsFileName.c_str(), ios::out);
		ofstream jMapStream((c->workingDir + "/junctions.map.fa").c_str(), ios::out);
		unordered_map<string, char*> junctExonSeq;

		for (auto cluster : clusters) {
			Gene *gA = cluster->gA, *gB = cluster->gB;
			bool isReadThruAB = Gene::IsReadThru(gA, gB, c->readThruDist);
			bool isReadThruBA = Gene::IsReadThru(gB, gA, c->readThruDist);

			unordered_map<int, Exon*> candExonsA, candExonsB;
			for (auto iso : gA->isoforms) {
				for (int i = 0; i < iso->exons.size(); i++) {
					Exon *ex = iso->exons[i];
					if (cluster->exonicReadCount.find(ex->exonId) !=
						cluster->exonicReadCount.end()) {
						for (int k = max(i - neighbor_cnt, 0);
							k < min((int)iso->exons.size(), i + neighbor_cnt + 1);
							k++)
							candExonsA[iso->exons[k]->exonId] = iso->exons[k];
					}
				}
			}
			for (auto iso : gB->isoforms) {
				for (int i = 0; i < iso->exons.size(); i++) {
					Exon *ex = iso->exons[i];
					if (cluster->exonicReadCount.find(ex->exonId) !=
						cluster->exonicReadCount.end()) {
						for (int k = max(i - neighbor_cnt, 0);
							k < min((int)iso->exons.size(), i + neighbor_cnt + 1); 
							k++)
							candExonsB[iso->exons[k]->exonId] = iso->exons[k];
					}
				}
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
								if (seq5p.length() < 30 || seq3p.length() < 30)
									continue;

								string keyA = lexical_cast<string>(exA->exonId) +
									"_5p_" + lexical_cast<string>(idx);
								string keyB = lexical_cast<string>(exB->exonId) +
									"_3p_" + lexical_cast<string>(idx);
								
								char *c_seq5p, *c_seq3p;
								if (junctExonSeq.find(keyA) == junctExonSeq.end()) {
									c_seq5p = new char[seq5p.length() + 1];
									strncpy(c_seq5p, seq5p.c_str(), seq5p.length());
									c_seq5p[seq5p.length()] = '\0';

									jMapStream << ">" << keyA << endl;
									jMapStream << c_seq5p << endl;
									junctExonSeq[keyA] = c_seq5p;
								}
								else
									c_seq5p = junctExonSeq[keyA];

								if (junctExonSeq.find(keyB) == junctExonSeq.end()) {
									c_seq3p = new char[seq3p.length() + 1];
									strncpy(c_seq3p, seq3p.c_str(), seq3p.length());
									c_seq3p[seq3p.length()] = '\0';

									jMapStream << ">" << keyB << endl;
									jMapStream << c_seq3p << endl;
									junctExonSeq[keyB] = c_seq3p;
								}
								else
									c_seq3p = junctExonSeq[keyB];

								Junction *junct = new Junction(cluster->clusterId,
									cluster, exA, exB, "F", idx, c_seq5p, c_seq3p);
								ss << ">" << junct->GetJunctionName() << endl;
								ss << c_seq5p << c_seq3p << endl;
								cluster->junctions.push_back(junct);

								idx++;
							}
						}
					}
					if (!isReadThruBA) {
						vector<string> exB_5p_seq = exB->GetJunctionSeqAs5pExon();
						vector<string> exA_3p_seq = exA->GetJunctionSeqAs3pExon();
						int idx = 1;
						for (auto seq5p : exB_5p_seq) {
							for (auto seq3p : exA_3p_seq)
							{
								if (seq5p.length() < 30 || seq3p.length() < 30)
									continue;

								string keyB = lexical_cast<string>(exB->exonId) +
									"_5p_" + lexical_cast<string>(idx);
								string keyA = lexical_cast<string>(exA->exonId) +
									"_3p_" + lexical_cast<string>(idx);

								char *c_seq5p, *c_seq3p;
								if (junctExonSeq.find(keyB) == junctExonSeq.end()) {
									c_seq5p = new char[seq5p.length() + 1];
									strncpy(c_seq5p, seq5p.c_str(), seq5p.length());
									c_seq5p[seq5p.length()] = '\0';

									jMapStream << ">" << keyB << endl;
									jMapStream << c_seq5p << endl;
									junctExonSeq[keyB] = c_seq5p;
								}
								else
									c_seq5p = junctExonSeq[keyB];

								if (junctExonSeq.find(keyA) == junctExonSeq.end()) {
									c_seq3p = new char[seq3p.length() + 1];
									strncpy(c_seq3p, seq3p.c_str(), seq3p.length());
									c_seq3p[seq3p.length()] = '\0';

									jMapStream << ">" << keyA << endl;
									jMapStream << c_seq3p << endl;
									junctExonSeq[keyA] = c_seq3p;
								}
								else
									c_seq3p = junctExonSeq[keyA];

								Junction *junct = new Junction(cluster->clusterId,
									cluster, exB, exA, "R", idx, c_seq5p, c_seq3p);
								ss << ">" << junct->GetJunctionName() << endl;
								ss << c_seq5p << c_seq3p << endl;
								cluster->junctions.push_back(junct);

								idx++;
							}
						}
					}
					junctionsFileStream << ss.str();
				}
			}
		}
		junctionsFileStream.close();
		jMapStream.close();
		return junctionCount;
	}

	int JunctionAligner::LoadJunctionsIntoClusters(
		vector<DiscordantCluster *> clusters)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		int junctionsCount = 0;
		unordered_map<int, DiscordantCluster*> clustersMap;
		for (auto cluster : clusters) {
			cluster->junctions.clear();
			clustersMap[cluster->clusterId] = cluster;
		}

		//Load junctions.map.fa
		unordered_map<string, char*> junctionSeqs;
		ifstream jMapStream((c->workingDir + "/junctions.map.fa").c_str());
		for (string key, seq; getline(jMapStream, key) && getline(jMapStream, seq);){
			key = key.substr(1, key.length() - 1);
			if (junctionSeqs.find(key) == junctionSeqs.end()) {
				//Note: this will never get deleted; clean up is not necesssary now;
				char *s = new char[seq.length() + 1];
				strncpy(s, seq.c_str(), seq.length());
				s[seq.length()] = '\0';
				junctionSeqs[key] = s;
			}
		}

		ifstream junctionsStream(GetJunctionsFilename().c_str());
		for (string junctName, junctSeq; getline(junctionsStream, junctName) 
			&& getline(junctionsStream, junctSeq);)
		{
			try {
				if ((junctName.c_str())[0] == '>')
					junctName = junctName.substr(1, junctName.length() - 1);
				auto sp = Utils::SplitToVector(junctName, "_");
				Junction *junct = new Junction();
				if (!Junction::ConstructJunctionFromName(junctName, junct)) {
					BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error loading "
						<< "junctions from: " << GetJunctionsFilename()
						<< endl << " at junct: " << junctName << ". Cannot "
						<< "find exons with the ExonIds in the junction";
					exit(1);
				}
				string keyA = lexical_cast<string>(junct->ex5p->exonId) + "_5p_" + 
					lexical_cast<string>(junct->junctionId);
				string keyB = lexical_cast<string>(junct->ex3p->exonId) + "_3p_" +
					lexical_cast<string>(junct->junctionId);
				if (junctionSeqs.find(keyA) == junctionSeqs.end() || 
					junctionSeqs.find(keyA) == junctionSeqs.end()) 
				{
					BOOST_LOG_CHANNEL(logger::get(), "Main") << " Error retrieving: "
						<< junctName << " junction from junctions.map.fa ";
				}
				junct->c_seq5p = junctionSeqs[keyA];
				junct->c_seq3p = junctionSeqs[keyB];
				junctionsCount++;
				junct->clusterPtr = clustersMap[junct->clusterId];
				clustersMap[junct->clusterId]->junctions.push_back(junct);
			}
			catch (std::exception &e) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error loading "
					<< "junctions from: " << GetJunctionsFilename()
					<< endl << " at junct: " << junctName << ". Exiting";
				exit(1);
			}
		}
		return junctionsCount;
	}

	string JunctionAligner::GetJunctionsFilename()
	{
		Config *c = Config::GetConfig();
		return c->workingDir + "junctions.fa";
	}

	string JunctionAligner::GetJunctionReadsFqFilename(int end)
	{
		Config *c = Config::GetConfig();
		return (end == 1 ? (c->workingDir + "junctions.alignments.reads.1.fastq")
			: (c->workingDir + "junctions.alignments.reads.2.fastq"));
	}

	string JunctionAligner::GetJunctionAlignmentsFilename()
	{
		Config *c = Config::GetConfig();
		return c->workingDir + "junctions.alignments";
	}

	string JunctionAligner::GetJunctionAlignmentsBamFilename()
	{
		Config *c = Config::GetConfig();
		return c->workingDir + "junctions.alignments.reads.bam";
	}

	void JunctionAligner::AlignToJunctions(Module *mod)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		int numJuncts = (Utils::LineCount(GetJunctionsFilename())) / 2;
		if (numJuncts == 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Junctions file: "
				<< GetJunctionsFilename() << " is empty";
			exit(1);
		}
		
		//Allow a maximum of 6 splits; with each split containing at least 2 cores;
		if (numJuncts > (6 * MAX_JUNCTIONS_PER_SPLIT))
			mod->NumSplits = (numJuncts / MAX_JUNCTIONS_PER_SPLIT) + 1;
		int linesPerSplit = ((numJuncts + mod->NumSplits) / mod->NumSplits) * 2;

		//split junctions file by linesPerSplit
		boost::filesystem::create_directories( c->workingDir + "junct_aligns/" );
		char buf[2000];
		sprintf(buf, "split -l %d -d -a 2 %s %s/junct_aligns/junctions.fa.", 
			linesPerSplit, GetJunctionsFilename().c_str(), c->workingDir.c_str());
		if (Utils::ExecuteCommand(buf).exit_code != 0){
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error splitting junctions "
				<< "file. Exiting.";
			exit(1);
		}
		
		ComputePerTask cpt = ComputePerTask::CalculateComputePerTask(6, 2, 4);

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Mapping initially unaligned "
			<< "reads in " << mod->NumSplits << " splits to candidate junctions";

		vector<boost::thread *> threads;
		vector<bool> thread_finished;
		int finishedJobs = 0;
		bool queueFull = false;
		for (int i = 0; i < mod->NumSplits; i++) {
			string threadIdStr = lexical_cast<string>(i);
			string chan = channel + ".JunctionAligner.split." + threadIdStr;
			Logger::RegisterChannel(c->sampleOutputLogDir + chan + ".log", chan);

			threads.push_back(new boost::thread(AlignToJunctions_worker, 
				threadIdStr, cpt.numCoresPerSplit, chan));
			thread_finished.push_back(false);

			if ((i + 1) == cpt.numSplits)
				queueFull = true;
			for (int j = 0; queueFull && j <= i; j++) {
				if (!thread_finished[j] 
					&& threads[j]->try_join_for(boost::chrono::seconds(1))) {					
					thread_finished[j] = true;
					break;
				}
				if (j == i) j = 0;				
				boost::this_thread::sleep(boost::posix_time::seconds(1));	
			}
		}

		for (int j = 0; j < mod->NumSplits; j++) {
			threads[j]->join();
			delete threads[j];
		}
	}

	void JunctionAligner::AlignToJunctions_worker(string thread, int cores, 
		string chan)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\tStarted " << "split # " << thread;
		string suffix = (boost::format("%02d") % stoi(thread)).str();
		string junctionsFile = 
			c->workingDir + "junct_aligns/junctions.fa." + suffix;

		char buf[2000];
		sprintf(buf, "%s %s %s", c->bowtie2buildPath.c_str(), 
			junctionsFile.c_str(), junctionsFile.c_str());
		if (Utils::ExecuteCommand(buf, chan).exit_code != 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occured while "
				<< "creating index for junctions file: " << junctionsFile
				<< endl << "Exiting";
			exit(1);
		}

		string fqCsv = "";
		for (auto fq : c->ExtractUnaligned_fq1)
			fqCsv += fq + ",";
		for (auto fq : c->ExtractUnaligned_fq2)
			fqCsv += fq + ",";

		sprintf(buf, "%s -p %d --local --score-min G,25,11 -k 20 --no-unal -x %s "
			"-U %s | %s %s.alignments 10 %f 1 ", c->bowtie2Path.c_str(), cores,
			junctionsFile.c_str(), fqCsv.c_str(), c->filterJunctOutputBinary.c_str(),
			junctionsFile.c_str(), c->maxJunctAlignErrorRate);
		if (Utils::ExecuteCommand(buf, chan).exit_code != 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occured while "
				<< "aligning to junctions. Exiting.";
			exit(1);
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\tFinished split # " << thread;
	}

	void JunctionAligner::GenerateFastqForJunctAlignments(Module *mod)
	{
		Config *c = Config::GetConfig();

		stringstream ss;
		ss << "cat ";
		for (int i = 0; i < mod->NumSplits; i++) {
			ss << c->workingDir << "junct_aligns/junctions.fa."
				<< (boost::format("%02d") % i).str() << ".alignments ";
		}

		ss << " > " << GetJunctionAlignmentsFilename() << ".tmp";
		if (Utils::ExecuteCommand(ss.str().c_str()).exit_code != 0){
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error compiling junction alignments";
			exit(1);
		}

		unordered_map<string, bool> readNames;
		ifstream junctFileStream(
			(GetJunctionAlignmentsFilename() + ".tmp").c_str());
		for (std::string str; getline(junctFileStream, str);) {
			auto sp = Utils::SplitToVector(str, "\t");
			string readName = Read::TrimReadName(sp[7]);
			if (readNames.find(readName) != readNames.end())
				continue;
			readNames[readName] = true;
		}

		unordered_map<string, vector<PairedRead> > pcrCheck;
		for (vector<string>::size_type i = 0; 
			i < c->ExtractUnaligned_fq1.size(); i++) 
		{
			string fq1 = c->ExtractUnaligned_fq1[i];
			string fq2 = c->ExtractUnaligned_fq2[i];

			if (!Utils::FileExists(fq1) || !Utils::FileExists(fq2))
				break;

			FastqParser *unalignedFP = new FastqParser(fq1, fq2);
			PairedRead pr;
			while (unalignedFP->GetNextPairedRead(&pr)) {
				string readName = Read::TrimReadName(pr.FirstRead.ReadName);
				if (readNames.find(readName) != readNames.end()) 
					pcrCheck[pr.FirstRead.Sequence.substr(0, 36) + "_" +
						pr.SecondRead.Sequence.substr(0, 36)].push_back(pr);
			}
		}

		unordered_map<string, bool> uniqueReads;
		ofstream fq1Out(GetJunctionReadsFqFilename(1).c_str(), ios::out);
		ofstream fq2Out(GetJunctionReadsFqFilename(2).c_str(), ios::out);
		for (auto iter : pcrCheck) {
			auto prs = iter.second;
			vector<int> keepReads;
			for (int i = 0; i < prs.size(); i++) {
				bool eliminate = false;
				for (int j = i + 1; j < prs.size() && !eliminate; j++) {
					const char *i1, *i2, *j1, *j2;
					i1 = prs[i].FirstRead.Sequence.c_str();
					i2  = prs[i].SecondRead.Sequence.c_str();
					j1 = prs[j].FirstRead.Sequence.c_str();
					j2 = prs[j].SecondRead.Sequence.c_str();
					if (strncmp(i1, j1, min<size_t>(strlen(i1), strlen(j1))) != 0 ||
						strncmp(i2, j2, min<size_t>(strlen(i2), strlen(j2))) != 0)
					{
						eliminate = true;
					}
				}
				keepReads.push_back(i);
			}
			for (auto keep : keepReads) {
				PairedRead pr = prs[keep];
				char b1[500], b2[500];
				pr.FirstRead.GetFastqString(b1), pr.SecondRead.GetFastqString(b2);
				fq1Out.write(b1, strlen(b1)), fq2Out.write(b2, strlen(b2));
				uniqueReads[Read::TrimReadName(pr.FirstRead.ReadName)] = true;
			}
		}
		fq1Out.close();
		fq2Out.close();

		//now filter out duplicate alignments in the junctionAlignments File
		ofstream junctsDupsRemOut(
			GetJunctionAlignmentsFilename().c_str(), ios::out);
		junctFileStream.clear();
		junctFileStream.seekg(0, ios::beg);
		for (std::string str; getline(junctFileStream, str);){
			auto sp = Utils::SplitToVector(str, "\t");
			if (uniqueReads.find(Read::TrimReadName(sp[7])) != uniqueReads.end())
				junctsDupsRemOut << str << endl;
		}

		// Align junction reads to transcriptome 
		char buf[1000];
		sprintf(buf, "%s -x %s -p %d -k 20 --very-sensitive-local --score-min G,25,11 "
			"-U %s,%s --quiet --no-unal | %s view -S -b - -o %s",
			c->bowtie2Path.c_str(), c->bowtie2AllIsoformIndex.c_str(), 
			min(12, c->maxCores), GetJunctionReadsFqFilename(1).c_str(),
			GetJunctionReadsFqFilename(2).c_str(), c->samtoolsPath.c_str(),
			GetJunctionAlignmentsBamFilename().c_str());
		if (Utils::ExecuteCommand(buf).exit_code != 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error aligning junction reads: "
				<< GetJunctionReadsFqFilename(1) << ","
				<< GetJunctionReadsFqFilename(2) << " to transcriptome. Exiting.";
			exit(1);
		}
	}

	int JunctionAligner::LoadAnchorReads(vector<DiscordantCluster *> clusters)
	{
		Config *c = Config::GetConfig();
		GeneModel *gm = GeneModel::GetGeneModel();

		int anchorCount = 0;
		unordered_map<int, PairedRead > readsFqByName;
		FastqParser *fp = new FastqParser(GetJunctionReadsFqFilename(1),
			GetJunctionReadsFqFilename(2));

		PairedRead pr;
		while (fp->GetNextPairedRead(&pr)) 
			readsFqByName[stoi(pr.GetTrimmedReadName())] = pr;

		unordered_map<string, Junction *> junctionsMap;
		for (auto cluster : clusters)
			for (auto junct : cluster->junctions)
				junctionsMap[junct->GetJunctionName()] = junct;

		//Load alignments
		unordered_map<string, vector<AnchorRead *> > anchorReadsByReadName;
		ifstream alignStream(GetJunctionAlignmentsFilename().c_str());
		for (std::string str; getline(alignStream, str);){
			auto sp = Utils::SplitToVector(str, "\t");
			string readname = sp[7], readname_trimmed = Read::TrimReadName(sp[7]);
			if (readsFqByName.find(stoi(readname_trimmed)) == readsFqByName.end()){
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error occurred while "
					<< "compiling anchor reads. .alignments and .reads.[1/2].fastq "
					<< "are out of sync. Exiting";
				exit(1);
			}
			PairedRead pr = readsFqByName[stoi(Read::TrimReadName(readname))];
			Junction *junct = junctionsMap[sp[0]];
			AnchorRead *ar = new AnchorRead(junct);
			if (pr.FirstRead.ReadName == readname) 
				ar->splitRead = pr.FirstRead, ar->otherRead = pr.SecondRead;
			else if (pr.SecondRead.ReadName == readname)
				ar->splitRead = pr.SecondRead, ar->otherRead = pr.FirstRead;
			else {
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "Anchor read not found in LoadAnchorReads";
				exit(1);
			}

			SplitRead *sr = &ar->splitRead;
			sr->cigar = sp[1], sr->strand = sp[2]; sr->flag = stoi(sp[3]);
			sr->position = stoi(sp[4]), sr->mismatches = stoi(sp[10]);
			sr->delBefore = stoi(sp[11]), sr->delAfter = stoi(sp[12]);
			sr->alignedSequence = sp[8].substr(sr->delBefore, 
				sp[8].length() - sr->delBefore - sr->delAfter);
			sr->alignedLength = sr->alignedSequence.length();
			sr->overhang5p = junct->breakpoint - sr->position;
			sr->overhang3p = sr->alignedSequence.length() - sr->overhang5p;
			sr->mismatchInAnchor = (sp[15] == "1" ? true : false);
			junct->anchorReads.push_back(ar);
			anchorReadsByReadName[readname_trimmed].push_back(ar);
			anchorCount++;
		}

		//Now load the AllIsoform alignments for the junction.reads
		BamTools::BamReader reader;
		BamTools::BamAlignment al;
		reader.Open(GetJunctionAlignmentsBamFilename());
		BamTools::RefVector refvect = reader.GetReferenceData();
		while (reader.GetNextAlignment(al)){
			string readname_trimmed = Read::TrimReadName(al.Name);
			vector<AnchorRead *> *ars = &anchorReadsByReadName[readname_trimmed];
			
			for (auto ar : (*ars)) {
				ReadAlignment readAlign;
				string geneIdStr = refvect[al.RefID].RefName.substr(1, 
					refvect[al.RefID].RefName.length()-1);
				readAlign.gene = gm->GenesMap[stoi(geneIdStr)];
				al.GetTag("NM", readAlign.mismatches);
				readAlign.position = al.Position;
				if (ar->splitRead.ReadName == al.Name)
					ar->splitRead.AddAlignment(readAlign);
				else
					ar->otherRead.AddAlignment(readAlign);
			}
		}
		return anchorCount;
	}
};