#ifndef FUSION_COMPILER_H
#define FUSION_COMPILER_H

#pragma once

#include <cstdio> 
#include <cmath>

#include "Utils.h"
#include "Config.h"
#include "GeneModel.h"
#include "Read.h"
#include "FusionJunction.h"
#include "DiscordantClusterFinder.h"
#include "FusionQuant.h"

using namespace std;

namespace MOJO {
	class FusionCompiler
	{
		private:
			static double GetDinucleotideEntropy(string sequence);

		public:
			static void Run(vector<DiscordantCluster *> clusters);
	};
};

#endif
