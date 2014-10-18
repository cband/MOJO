#ifndef JUNCTIONALIGNER_H
#define JUNCTIONALIGNER_H

#pragma once

#include <boost/format.hpp>

#include "BamReader.h"
#include "BamAlignment.h"

#include "Utils.h"
#include "Config.h"
#include "Logger.h"
#include "GeneModel.h"
#include "Read.h"
#include "DiscordantClusterFinder.h"
#include "FusionJunction.h"

namespace MOJO 
{
	class JunctionAligner
	{
		private:
			static string channel;

			static int BuildJunctions( vector<DiscordantCluster *> clusters);

			static int LoadJunctionsIntoClusters(vector<DiscordantCluster *> clusters);

			static void AlignToJunctions(Module *m);

			static void AlignToJunctions_worker(string thread, int cores, string chan);

			static void GenerateFastqForJunctAlignments(Module *m);

			static void ConstructAnchorReads();

			static int LoadAnchorReads(vector<DiscordantCluster *> clusters);

			static string GetJunctionsFilename();

			static string GetJunctionReadsFqFilename(int end);

			static string GetJunctionAlignmentsFilename();

			static string GetJunctionAlignmentsBamFilename();

		public:
			static void Run(vector<DiscordantCluster *> clusters);
	};
};

#endif