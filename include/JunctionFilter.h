#ifndef JUNCTION_FILTER_H
#define JUNCTION_FILTER_H

#pragma once

#include "BamReader.h"
#include "BamAlignment.h"

#include "Utils.h"
#include "Config.h"
#include "GeneModel.h"
#include "FusionJunction.h"
#include "DiscordantClusterFinder.h"
#include "PSLParser.h"

using AnchorReadMap = boost::unordered_map <string, MOJO::AnchorRead*>;
using JunctionMap = boost::unordered_map < string, MOJO::Junction* > ;

namespace MOJO 
{
	enum FilterStringency 
	{
		LOW, 
		MEDIUM, 
		HIGH
	};

	class JunctionFilter
	{
		private:
			static void GenerateFastaForAnchorReads(
				AnchorReadMap *arMap, string outputFile);
	
			static void GenerateFastaForJunctions(AnchorReadMap *arMap,
				string outputFa);

			static void FilterSplitReadsByExome(string faFile,
				AnchorReadMap *arMap, FilterStringency s, Module *m);
			
			static void FilterSplitReadsByExome_worker(string faFile,
				FilterStringency s);
			
			static void FilterSplitReadsByGenome(string faFile,
				AnchorReadMap *arMap, FilterStringency s, Module *m);
			
			static void FilterSplitReadsByGenome_worker(string faFile, 
				string contig, FilterStringency s);
			
			static void FilterSplitReadsByGenome_post(string pslFile,
				AnchorReadMap *arMap);
			
			static void FilterSpuriousJunctions(AnchorReadMap *arMap, Module *m);
			
			static void FilterSpuriousJunctions_post(string pslFile,
				AnchorReadMap *arMap);

			static string GetBlatParamsForStringency(FilterStringency f);

			static void FilterSpuriousScrambledReads(AnchorReadMap *arMap, 
				Module *m);
			
			static void FilterSpuriousScrambledReads_worker(string query,
				string target, int cores);
			
			static void FilterJunctionPCRdups(AnchorReadMap *arMap);

		public:
			static void Run(vector<DiscordantCluster *> clusters);
			
			static void CleanUp();
	};
}

#endif