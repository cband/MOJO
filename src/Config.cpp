
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "Config.h"

using namespace std;
namespace po = boost::program_options;

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = "Main"));

	Config Config::MOJORunConf;

	Config::Config() { readcount = 0; }

	bool Config::LoadConfiguration(int argc, char *argv[]) 
	{
		po::options_description options_desc("Available options", 100);
		boost::program_options::arg = "";

		options_desc.add_options()
			("config", 
				po::value<string>(&configFile), 
				"[Req] MOJO configuration file")
			("output_dir", 
				po::value<string>(&outputDir),
				"[Req] output directory ")
			("sample_name", 
				po::value<string>(&sampleName),
				"[Req] sample name" )
			("fq1", 
				po::value<string>(&fastqFirstEndcsv),
				"[Req] first-end fastq(.gz) files; comma-separated")
			("fq2", 
				po::value<string>(&fastqSecondEndcsv),
				"[Req] second-end fastq(.gz) files; comma-separated\n")


			("fastq_encoding",
			po::value<string>(&fastqEncodingString),
				"'sanger' or 'illumina' (1.3+) only. If empty, \n"
				"encoding will be inferred. default:[]")
			("cores", 
				po::value<int>(&maxCores)->default_value(0, ""),
				"num cores. If empty, MOJO will infer this. default: [] ")
			("mem", 
				po::value<int>(&maxMem)->default_value(0, ""),
				"memory available (in GB). If empty, MOJO will attempt \n"
				"to infer total system memory. default: [] " )

			("min_span", 
				po::value<string>(&minSpanFunct)->default_value("2,2,80000000", ""),
				"Threshold for minimum number of discordant reads required \n"
				"for each fusion candidate. configured as function of \n"
				"library size: ceil(R + X*max(0,ln(librarSize/Y))), \n"
				"R: constant, X: coefficient and Y a scaling factor \n"
				"for library size. default (R,X,Y): [2,2,80000000]")

			("read_through", 
				po::value<int>(&readThruDist)->default_value(200000, ""),
				"read-through cut-off dist (in bp). default: [200000]")
			
			("junct_mismatch",
				po::value<double>(&maxJunctAlignErrorRate)->default_value(0.03, ""),
				"ignore reads aligning to junctions with mismatch \n"
				"rate higher than this. default: [0.03]\n")
			; 
		
		po::variables_map options_vm, desc_vm;
		po::store(po::parse_command_line(argc, argv, options_desc), options_vm);
		po::notify(options_vm);

		if (argc == 1 || options_vm.count("help")) {
			std::cout << endl 
				<< "Program: MOJO (Minimum Overlap Junction Optimizer)" << endl
				<< "Version: " << MOJO::version << std::endl << std::endl
				<< "Usage: MOJO --config <configfile> [options]" 
				<< std::endl << std::endl << options_desc << std::endl;
			exit(0);
		}
		
		if (!options_vm.count("config")) {
			std::cerr << "Error: --config is options" << endl;
			exit(1);
		}

		options_desc.add_options()
			("mojo_install_dir,m", 
				po::value<string>(&MOJOInstallDir), 
				"path to MOJO installation directory")
			
			("mojo_reference_dir,r", 
				po::value<string>(&MOJOReferenceDir),
				"path to MOJO reference directory")
			
			("mojo_tools_dir,t", 
				po::value<string>(&MOJOToolsDir),
				"path to MOJO tools directory")

			("action_remove_temp", 
				po::value<bool>(&removeTemporaryFiles)->default_value(true), 
				" debug option. 1 to remove temporary files [1]")

			//paths to binaries or supporting tools
			("samtools_binary", 
				po::value<string>(&samtoolsPath), 
				" full path to samtools binary")

			("bwa_binary", 
				po::value<string>(&bwaPath), 
				" full path to bwa binary")
			
			("bowtie2_binary", 
				po::value<string>(&bowtie2Path), 
				" full path to bowtie2 binary")
			
			("bowtie2_build_binary", 
				po::value<string>(&bowtie2buildPath), 
				" full path to bowtie2 binary")
			
			("blat_binary", 
				po::value<string>(&blatPath), 
				" full path to blat binary")

			//index files
			("bwa_transcriptome_index", 
				po::value<string>(&bwaTranscriptomeIndex), 
				" see configuration file for options")

			("bowtie2_all_isoforms_index", 
				po::value<string>(&bowtie2AllIsoformIndex), 
				" see configuration file for options")

			("bowtie2_genome_index", 
				po::value<string>(&bowtie2GenomeIndex),
				" see configuration file for options")
			
			("blat_genome_2bit", 
				po::value<string>(&blatGenome2Bit), 
				" see configuration file for options")
			
			("blat_reference_junctions_2bit", 
				po::value<string>(&blatReferenceJunctions2Bit), 
				" see configuration file for options")
			
			("blat_filter_chroms_dir", 
				po::value<string>(&blatFilterChromsDir), 
				" see configuration file for options")
			
			("blat_filter_chroms_csv", 
				po::value<string>(&blatFilterChromsCsv)->default_value("default"), 
				" see configuration file for options")
			
			("max_bwa_mem", 
				po::value<int>(&maxBwaMem)->default_value(6), 
				"Max memory to allocate (in GB) for each alignment run (6 for humans) [6]")

			//reference files
			("master_gene_file", 
				po::value<string>(&masterGeneFile)->default_value("default"), 
				" ")

			("master_exon_file", 
				po::value<string>(&masterExonFile)->default_value("default"), 
				" ")

			("master_isoform_file", 
				po::value<string>(&masterIsoformFile)->default_value("default"), 
				" ")

			("megablast_output_file", 
				po::value<string>(&megablastOutputFile)->default_value("default"), 
				" ")

			("repeat_masker_file", 
				po::value<string>(&repeatMaskerFile)->default_value("default"), 
				" ")

			("split_fastq_binary", 
				po::value<string>(&splitFastqBinary)->default_value("default"), 
				" ")

			("filter_junct_output_binary", 
				po::value<string>(&filterJunctOutputBinary)->default_value("default"), 
				" ")
			;

		po::store(po::parse_command_line(argc, argv, options_desc), options_vm);
		po::notify(options_vm);

		ifstream ifs(configFile.c_str());
		if (!ifs.is_open()) {
			std::cerr << "Error: cannot open configuration file: " 
				<< configFile << endl;
			exit(1);
		}
		po::store(po::parse_config_file(ifs, options_desc), options_vm);
		po::notify(options_vm);

		if (!options_vm.count("output_dir")) {
			std::cerr << "Missing required parameter: --output_dir " << endl;
			exit(1);
		}

		if (!options_vm.count("sample_name")) {
			std::cerr << "Missing required parameter: --sample_name " << endl;
			exit(1);
		}

		if (!options_vm.count("fq1")) {
			std::cerr << "Missing required parameter: --fq1 " << endl;
			exit(1);
		}

		if (!options_vm.count("fq2")) {
			std::cerr << "Missing required parameter: --fq2 " << endl;
			exit(1);
		}

		if (!options_vm.count("mojo_install_dir")) {
			std::cerr << "Missing required parameter: --mojo_install_dir" << endl;
			exit(1);
		}

		if (!options_vm.count("cores")) {
			std::cerr << "Missing required parameter: --cores " << endl;
			exit(1);
		}

		if (outputDir.back() != '/')	outputDir += "/";
		outputDir += sampleName + "/";
		workingDir = outputDir + "tmp/";
		sampleOutputPrefix = workingDir + "/" + sampleName + ".";
		sampleOutputLogDir = outputDir + "logs/";

		boost::filesystem::create_directories(workingDir);
		boost::filesystem::create_directories(sampleOutputLogDir);

		//Check input for fastqEncoding
		if (fastqEncodingString != "sanger"
			&& fastqEncodingString != "illumina" && fastqEncodingString != "")
		{
			std::cerr << "Invalid option '" << fastqEncodingString
				<< "' for fastq_encoding. Valid options are  'sanger' and 'illumina'. "
				<< " Leave empty to let MOJO infer encoding" << endl, exit(1);
		}

		//Get CPU count and Memory information
		if (maxCores <= 0)  {
			maxCores = SystemInfo::GetSystemCpuCount();
			if (maxCores <= 0) {
				std::cerr << "Error: unable to get number of processor cores. "
					<< "Please specify using --cores" << endl;
				exit(1);
			}
		}
		if (maxMem <= 0) {
			maxMem = (int) ( 0.85 * SystemInfo::GetTotalSystemMemory() );
			if (maxMem <= 0)  {
				std::cerr << "Error: unable to get total available system memory. "
					<< " (min 2GB) " 
					<< "Please specify using --mem" << endl;
				exit(1);
			}
		}

		if (fastqFirstEndcsv.back() == ',')
			fastqFirstEndcsv.pop_back();

		if (fastqSecondEndcsv.back() == ',')
			fastqSecondEndcsv.pop_back();

		//Verify proper MOJO setup
		VerifyFilePathParam(MOJOInstallDir, "mojo_install_directory");
		VerifyFilePathParam(outputDir, "output_dir");

		boost::replace_all(samtoolsPath, "<TOOLS>", MOJOToolsDir);
		boost::replace_all(bwaPath, "<TOOLS>", MOJOToolsDir);
		boost::replace_all(bowtie2Path, "<TOOLS>", MOJOToolsDir);
		boost::replace_all(bowtie2buildPath, "<TOOLS>", MOJOToolsDir);
		boost::replace_all(blatPath, "<TOOLS>", MOJOToolsDir);

		VerifyFilePathParam(samtoolsPath, 
			"samtools_binary");
		VerifyFilePathParam(bwaPath, 
			"bwa_binary");
		VerifyFilePathParam(bowtie2Path, 
			"bowtie2_binary");
		VerifyFilePathParam(bowtie2buildPath, 
			"bowtie2_build_binary");
		VerifyFilePathParam(blatPath, 
			"blat_binary");

		boost::replace_all(bwaTranscriptomeIndex, "<REF>", MOJOReferenceDir);
		boost::replace_all(bowtie2AllIsoformIndex, "<REF>", MOJOReferenceDir);
		boost::replace_all(bowtie2GenomeIndex, "<REF>", MOJOReferenceDir);
		boost::replace_all(blatGenome2Bit, "<REF>", MOJOReferenceDir);
		boost::replace_all(blatReferenceJunctions2Bit, "<REF>", MOJOReferenceDir);
		boost::replace_all(blatFilterChromsDir, "<REF>", MOJOReferenceDir);

		VerifyFilePathParam(bwaTranscriptomeIndex + ".amb", 
			"bwa_transcriptome_index");
		VerifyFilePathParam(bowtie2AllIsoformIndex + ".1.bt2", 
			"bowtie2_all_isoforms_index");
		VerifyFilePathParam(bowtie2GenomeIndex + ".1.bt2", 
			"bowtie2_genome_index");
		VerifyFilePathParam(blatGenome2Bit, 
			"blat_genome_2bit");
		VerifyFilePathParam(blatReferenceJunctions2Bit, 
			"blat_reference_junctions_2bit");
		VerifyFilePathParam(blatFilterChromsDir, 
			"blat_filter_chroms_dir");

		if (masterGeneFile == "default")
			masterGeneFile = MOJOReferenceDir + "/gene_model/Gene.txt";
		if (masterExonFile == "default")  
			masterExonFile = MOJOReferenceDir + "/gene_model/Exon.txt";
		if (masterIsoformFile == "default")  
			masterIsoformFile = MOJOReferenceDir + "/gene_model/Isoform.txt";
		if (megablastOutputFile == "default") 
			megablastOutputFile = MOJOReferenceDir + "/gene_model/gene2gene.megablast.txt";
		if (repeatMaskerFile == "default") 
			repeatMaskerFile = MOJOReferenceDir + "/gene_model/rmsk.regions.txt";
		if (filterJunctOutputBinary == "default") 
			filterJunctOutputBinary = MOJOInstallDir + "/FilterJunctAlignOutput";
		if (splitFastqBinary == "default") 
			splitFastqBinary = MOJOInstallDir + "/StreamNthFastqSplit";

		if (blatFilterChromsCsv == "default" || blatFilterChromsCsv == "" ) {
			blatFilterChromsCsv = "";
			//get a list 
			boost::filesystem::directory_iterator iterator(blatFilterChromsDir);
			for (; iterator != boost::filesystem::directory_iterator(); ++iterator) {
				auto fn = iterator->path().filename();
				if (fn.extension().string() == ".2bit" && 
					fn.string().substr(0, 3) == "chr") 
				{
					string n = fn.string();
					blatFilterChromsCsv += n.substr(0, n.length() - 5)  + ",";
				}
			}
		}

		VerifyFilePathParam(masterGeneFile, "master_gene_file");
		VerifyFilePathParam(masterExonFile, "master_exon_file");
		VerifyFilePathParam(masterIsoformFile, "master_isoform_file");
		VerifyFilePathParam(megablastOutputFile, "megablast_output_file");
		VerifyFilePathParam(repeatMaskerFile, "repeat_masker_file");
		VerifyFilePathParam(splitFastqBinary, "split_fastq_binary");
		VerifyFilePathParam(filterJunctOutputBinary, "filter_junct_output_binary");

		blatFilterChromsVect = Utils::SplitToVector(blatFilterChromsCsv, ",");
	
		vector<string> blatFilterFiles;
		unordered_map<string, string> blatFilterFiles_map;
		for (auto chrom : blatFilterChromsVect) {
			if (chrom == "")
				continue;
			string filename = blatFilterChromsDir + "/" + chrom + ".2bit";
			blatFilterFiles.push_back(filename);
			VerifyFilePathParam(filename, "blat_filter_chroms_dir:" + chrom);
			blatFilterFiles_map[filename] = chrom;
		}
		
		blatFilterFiles = Utils::SortListOfFilesBySize(blatFilterFiles);
		blatFilterChromsVect.clear();
		for (auto f : blatFilterFiles)
			blatFilterChromsVect.push_back(blatFilterFiles_map[f]);

		JunctionFilter_dir = workingDir + "/filters/";
		statusFilesDir = workingDir + "/status/";

		boost::filesystem::create_directories(JunctionFilter_dir);
		boost::filesystem::create_directories(statusFilesDir);
		boost::filesystem::create_directories(workingDir + "/mapstats/");

		statusFile_ExtractUnaligned = statusFilesDir + "ExtractUnaligned.status";
		statusFile_DiscordantFinder = statusFilesDir + "DiscordantFinder.status";
		statusFile_JunctionAligner = statusFilesDir + "JunctionAligner.status";
		statusFile_JunctionFilter = statusFilesDir + "JunctionFilter.status";
		statusFile_FusionStats = statusFilesDir + "FusionStats.status";
		readcountFile = workingDir + "/fastqs/readcount";

		IsConfigurationLoaded = true;
		return true;
	}
	void Config::VerifyFilePathParam(string path, string parameter_name) 
	{
		if (!Utils::FileExists(path)) {
			std:cerr << "file/directory not found: '" << parameter_name
				<< "':'" << path << "' " << endl;
			exit(1);
		}
	}
	
	// Check if the fastq files exist and infer fastq encoding.
	// Also, check if user supplied encoding and inferred matches.
	bool Config::PrerunFastqCheck() 
	{
		vector<FASTQ_ENCODING> encodings;
		auto fastq1s = Utils::SplitToVector(fastqFirstEndcsv, ",");
		auto fastq2s = Utils::SplitToVector(fastqSecondEndcsv, ",");
		for (string::size_type i = 0; i < fastq1s.size(); i++) {
			string end1File = fastq1s[i], end2File = fastq2s[i];
			
			if (!Utils::FileExists(end1File, true))
				exit(1);
			if (!Utils::FileExists(end2File, true))  
				exit(1);

			FastqParser *parser = new FastqParser(end1File, end2File);
			encodings.push_back(parser->GetFastqEncoding());
			firstEndFastqs.push_back(end1File);
			secondEndFastqs.push_back(end2File);
		}

		for (int i = 0; i < encodings.size(); i++) {
			if (encodings[i] == FASTQ_ENCODING::UNSUPPORTED) {
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "WARNING: Unsupported fastq encoding detected. "
					<< "Only 'sanger' and 'illumina' (1.3+) are supported. ";
				
				if (fastqEncodingString == "") {
					BOOST_LOG_CHANNEL(logger::get(), "Main")
						<< "Use --fastq_encoding to override inferred encoding";
					return false;
				}
			}
			if (encodings[i] != encodings[0]) {
				BOOST_LOG_CHANNEL(logger::get(), "Main") << "WARNING: each "
					<< "sequencing lane has different fastq encoding";
				
				if (fastqEncodingString == ""){
					BOOST_LOG_CHANNEL(logger::get(), "Main")
						<< "Use --fastq_encoding to override inferred encoding";
					return false;
				}
			}
		}

		string sampleEncode = FastqParser::GetFastqEncodingString(encodings[0]);
		if (fastqEncodingString != "" && fastqEncodingString != sampleEncode) {
			BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "WARNING: Fastq Encoding inconsistency detected. User supplied: '"
				<< fastqEncodingString << "', MOJO inferred: '" << sampleEncode
				<< "'. Proceeding with '" << fastqEncodingString << "' ";
		}
		if (fastqEncodingString == "" )
			fastqEncodingString = sampleEncode;

		BOOST_LOG_CHANNEL(logger::get(), "Main") << "Detected Fastq Encoding: "
			<< fastqEncodingString;
		return true;
	}

	Config *Config::GetConfig()
	{
		if (Config::MOJORunConf.IsConfigurationLoaded)
			return &Config::MOJORunConf;
		cerr << "Configuration file not loaded. Exiting." << endl;
		exit(1);
	}

	// ./fastqs/readcount file is generated by DiscordantReadFinder after
	// the initial alignment to the spliced transcriptome. This function
	// reads that file and returns the readcount
	int Config::GetTotalReadcount() 
	{
		if (readcount != 0)
			return readcount;

		ifstream rcFile( readcountFile.c_str());
		string s;
		getline(rcFile, s);
		//Parse out the read count;
		boost::match_results<std::string::const_iterator> results;
		boost::regex_match(s, results,
			boost::regex("Total read count:[\\s|\\t]+(\\d+),\\s.*"));

		if (!Utils::FileExists(readcountFile))
			return 0;

		try{
			readcount = stoi(results[1].str());
		}
		catch (std::exception e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "Error getting readcount from: " << readcountFile;
			exit(1);
		}
		return readcount;
	}

	//Splits input fastqs into specified number of splits and outputs to a 
	//named pipe.  This reduces the I/O overhead during the initial preprocessing
	//of the reads during which reads are trimmed (for a certain length or just
	//for quality)
	string Config::GenerateCatCmdsForSplit(int end, int split_id,
		int num_splits, bool silent )
	{
		Config *c = Config::GetConfig();
		vector<string> *fqs;
		if (end == 1)
			fqs = &c->firstEndFastqs;
		else
			fqs = &c->secondEndFastqs;

		int fq_idx = 0;
		stringstream ss, ssFqFifo;
		ssFqFifo << "cat ";
		for (auto fq : (*fqs)){
			string cat = "cat ";
			if (fq.substr(fq.length() - 3, 3) == ".gz")
				cat = "gzip -dc ";

			string fifo = c->workingDir + "/fastqs/input_"
				+ lexical_cast<string>(split_id)+"_"
				+ lexical_cast<string>(fq_idx)+"_"
				+ lexical_cast<string>(end)+".fifo";

			ss << "rm -f " << fifo << endl;
			ss << "sleep 1" << endl;
			ss << "mkfifo " << fifo << endl;
			ss << cat << " " << fq << " > " << fifo << " & " << endl;
			ss << "sleep 1" << endl;
			ssFqFifo << fifo << " ";
			fq_idx++;
		}

		string fqFile = c->workingDir + "/fastqs/split_" + to_string(split_id) +
			"_" + to_string(end) + ".fastq";
		string readCountFile = "-", readMappingFile = "-";
		if (!silent && split_id == 0)
			readCountFile = c->workingDir + "/fastqs/readcount";
		if (!silent && end == 1)
			readMappingFile = fqFile + ".readID";

		ss << "rm -f " << fqFile << endl;
		ss << "sleep 0.2" << endl;
		ss << "mkfifo " << fqFile << endl;
		ss << "sleep 0.2" << endl;

		ss << ssFqFifo.str() << " | " 
			<< c->splitFastqBinary << " " << num_splits << " " << split_id << " " 
			<< end << " 15 " << c->fastqEncodingString << " 1 " << readCountFile 
			<< " " << readMappingFile << " > " << fqFile << " & " << endl;

		return ss.str();
	}

	//Removes all temporary files.  Final fusion output will have the following files
	// .fusions, .fusions.pileup, .log, .discordants.genes and .junctions.alignments.
	// .discordants.genes has the number of unique/non-unique discordant reads mapping
	// between two genes
	void Config::FinalCleanup() 
	{
		char cmd[1000];
		
		sprintf(cmd, "cp %s %s/readcount", readcountFile.c_str(), 
			outputDir.c_str());
		Utils::ExecuteCommand(cmd);

		sprintf(cmd, "gzip -c %s/discordants.genes > %s/discordants.genes.gz", 
			workingDir.c_str(), outputDir.c_str());
		Utils::ExecuteCommand(cmd);

		sprintf(cmd, "gzip -c %s/junctions.alignments > %s/junctions.alignments.gz",
			workingDir.c_str(), outputDir.c_str());
		Utils::ExecuteCommand(cmd);

		if (removeTemporaryFiles) {
			sprintf(cmd, "rm -r %s", workingDir.c_str());
			Utils::ExecuteCommand(cmd);
		}
	}

	//MOJO is run in five discrete steps.  A statusFile is constructed
	//for each module to indicate if the processing of that module has
	//completed.  
	Module Module::LoadStatusFile(string statusFile)
	{
		Module m (statusFile);
		if (!Utils::FileExists(statusFile))
			return m;

		try {
			po::options_description options_status("Status Options", 100);
			boost::program_options::arg = "";
			options_status.add_options()
				("IsComplete", po::value<bool>(&m.IsComplete))
				("NumSplits", po::value<int>(&m.NumSplits))
				;

			ifstream ifs(statusFile.c_str());
			po::variables_map status_vm;
			po::store(po::parse_config_file(ifs, options_status), status_vm);
			po::notify(status_vm);
		}
		catch (std::exception &e) {
			cerr << "Error occurred while trying to load status file: " 
				<< statusFile << " Error: " << e.what() << endl;
			cerr << "Delete the output directory and re-start run." << endl;
			exit(1);
		}
		return m;
	}

	void Module::SaveStatusFile()
	{
		ofstream out(filename.c_str(), ios::out);
		out << "#### DO NOT MODIFY.  MOJO params." << endl;
		out << "IsComplete=" << (IsComplete ? "1" : "0") << endl;
		out << "NumSplits=" << NumSplits << endl;
		out.close();
	}
}