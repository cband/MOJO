#ifndef GENE_MODEL_H
#define GENE_MODEL_H

#pragma once

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include "Utils.h"
#include "GeneModelObjs.h"
#include "Config.h"

using namespace std;
using namespace boost;

namespace MOJO 
{
	class GeneModel 
	{
		private:
			static map<string, map<int, Gene*> > chromGenesRefSeq;

		public:
			static GeneModel gm;
			unordered_map<int, Gene*> GenesMap;
			unordered_map<int, Exon*> ExonsMap;
			//need to transition uc000abc.2 ids to custom ids 
			unordered_map<string, Isoform*> IsoformsMap;	
			bool IsGeneModelLoaded;

			GeneModel() : IsGeneModelLoaded(false) {};
			
			void LoadGeneModel();
			
			bool ValidateGeneModelsIntegrity();

			Gene* GetGeneFromMap(string geneIdStr);

			static GeneModel *GetGeneModel();

			static string GetGeneNameForCoordinates(string chr, int pos);
	};
}
#endif