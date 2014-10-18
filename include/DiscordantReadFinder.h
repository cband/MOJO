#ifndef DISCORDANT_READ_FINDER_H
#define DISCORDANT_READ_FINDER_H

#pragma once

#include <vector>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "Utils.h"
#include "Logger.h"
#include "Config.h"
#include "FastqParser.h"
#include "GeneModel.h"

#include "BamReader.h"
#include "BamWriter.h"
#include "BamAlignment.h"

using namespace BamTools;

namespace MOJO 
{
	//
	// Find Discordant Reads : To efficiently parallelize the steps in this process, 
	// input fastqs are first split. Number of splits is determined by maximum memory
	// consumption of each split in the alignments step.  
	// For human genome: ~6GB per split.
	//
	class DiscordantReadFinder 
	{
		public:
			static bool Run();
			
			static vector<BamAlignment> GetBwaAlignsFromXA(BamAlignment al,
				unordered_map<string, int> *geneToNameMap);

		private:
			static boost::mutex statsUpdateMutex1, statsUpdateMutex2;
			static string DChannel;
			static int TotalUnalignedReads, TotalDiscordantReads;

			DiscordantReadFinder() {}
			
			DiscordantReadFinder(const DiscordantReadFinder&);
			
			DiscordantReadFinder& operator=(const DiscordantReadFinder&);

			static ComputePerTask CalculateComputePerTask(int max_splits, 
				int min_cores_per_split, int max_mem_per_split);
			
			static void ExtractUnalignedReads(ComputePerTask cpt);
			
			static void ExtractUnalignedReads_worker(int thread, ComputePerTask cpt, 
				string chan);
			
			static void FindDiscordantReads(ComputePerTask cpt);
			
			static void FindDiscordantReads_worker(string threadIdStr, int tCores, 
				string chan);
			
			static void FindDiscordantReads_ByIteration(string fqPrefix, 
				string threadIdStr, int tCores, string chan);

			static void UpdateUnalignedReadCount(int cnt);
			
			static void UpdateDiscordantReadCount(int cnt);

			//Utility function
			static string BuildFlatReadString(BamAlignment fAl, BamAlignment sAl, 
				RefVector *refVect);
			
			static string GetFlatFilename();
	};
};
#endif

