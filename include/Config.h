#ifndef Config_H
#define Config_H

#include <cstring>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "MOJOConfig.h"
#include "Utils.h"
#include "FastqParser.h"

using namespace std;

#define CONTIG_LENGTH 80	// length of the sequence contributed by each end (5'/3') 
							// to the candidate fusion junction
#define ANCHOR_LENGTH 10	// minimum anchor length below which all junction aligning
							// reads are discarded

namespace MOJO 
{
	class Module
	{
		public:
			bool IsComplete;
			int NumSplits;
			string filename;

			Module() : IsComplete(false), NumSplits(0) {}

			Module(string f) : filename(f), IsComplete(false), NumSplits(0) {}
			
			void SaveStatusFile();
			
			static Module LoadStatusFile(string statusFile);
	};

	class Config 
	{
		public:
			static Config MOJORunConf;

			enum Modules {
				EXTRACT_UNALIGNED, 
				FIND_DISCORDANT, 
				JUNCTION_ALINGER, 
				JUNCTION_FILTER
			};

			bool IsConfigurationLoaded;

			//////// Configuration Parameters ////////
			bool runQuantification;

			string configFile;
			string MOJOInstallDir, MOJOReferenceDir, MOJOToolsDir;
			string sampleName, outputDir, workingDir, sampleOutputLogDir;
			string sampleOutputPrefix, fastqFirstEndcsv, fastqSecondEndcsv;

			//paths to binaries of supporting tools
			string samtoolsPath, bwaPath, bowtie2Path, bowtie2buildPath;
			string blatPath;

			//index files
			string bwaTranscriptomeIndex, bowtie2AllIsoformIndex;
			string bowtie2GenomeIndex, blatReferenceJunctions2Bit;
			string blatGenome2Bit, blatFilterChromsDir, blatFilterChromsCsv;

			//reference files
			string masterGeneFile;
			string masterExonFile, masterIsoformFile;
			
			//filtering files
			string megablastOutputFile, repeatMaskerFile;

			//Compute parameters
			int maxCores, maxMem, maxBwaMem;

			//run parameters
			string minSpanFunct;
			int minSpanCount, readThruDist; // , anchorMismatches;
			double maxJunctAlignErrorRate;

			//debug flags
			bool removeTemporaryFiles;

			//supplementary binaries
			string splitFastqBinary, filterJunctOutputBinary;

			//////// MOJO Environment Variables ////////
			vector<string> firstEndFastqs, secondEndFastqs;
			vector<string> blatFilterChromsVect;

			FASTQ_ENCODING fastqEncoding;
			string fastqEncodingString;

			string JunctionFilter_dir;
			vector<string> ExtractUnaligned_fq1, ExtractUnaligned_fq2;

			string statusFilesDir;
			string statusFile_ExtractUnaligned;
			string statusFile_DiscordantFinder;
			string statusFile_JunctionAligner;
			string statusFile_JunctionFilter;
			string statusFile_FusionStats;

			string readcountFile;

			int readcount;

			Config();
			
			bool LoadConfiguration(int argc, char *argv[]);
			
			void VerifyFilePathParam(string path, string parameter_name);
			
			bool PrerunEnvironmentCheck();		
			
			bool PrerunFastqCheck();
			
			static Config *GetConfig();
			
			static string GenerateCatCmdsForSplit(int end, int split_id,
				int num_splits, bool silent=false);

			int GetTotalReadcount();

			void FinalCleanup();
	};
};

#endif

