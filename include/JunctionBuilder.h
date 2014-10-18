#ifndef JUNCTIONBUILDER_H
#define JUNCTIONBUILDER_H

#pragma once

#include "BamReader.h"
#include "BamWriter.h"
#include "BamAlignment.h"

#include "Utils.h"
#include "Config.h"
#include "GeneModel.h"
#include "DiscordantClusterFinder.h"
#include "SWalign.h"

namespace MOJO
{
	class JunctionBuilder
	{
		public:
			static int BuildJunctions(vector<DiscordantCluster *> clusters);
	};

};

#endif