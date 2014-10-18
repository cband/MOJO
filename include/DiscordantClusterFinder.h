#ifndef DISCORDANT_CLUSTER_FINDER_H
#define DISCORDANT_CLUSTER_FINDER_H

#pragma once

#include "Utils.h"
#include "Config.h"
#include "GeneModel.h"
#include "FusionJunction.h"

namespace MOJO 
{
	class DiscordantCluster;

	class DiscordantClusterFinder 
	{
		private:
			static string GetFlatFilename();

		public:
			static int NextClusterID;
		
			static vector<DiscordantCluster *> LoadDiscordantClusters();
	};

	class DiscordantClusterAlignment 
	{
		public:
			int startA, startB;
			int endA, endB;
			string strandA, strandB;
			bool isUnique;

			DiscordantClusterAlignment() {};
			
			DiscordantClusterAlignment( string stA, string lenA, string strndA,
				string stB, string lenB, string strndB, string uniq)
			{
				startA = stoi(stA), endA = startA + stoi(lenA);
				startB = stoi(stB), endB = startB + stoi(lenB);
				strandA = strndA, strandB = strndB;
				isUnique = (bool)stoi(uniq);
			}
	};

	class DiscordantCluster {
	public:
		Gene *gA, *gB;
		int clusterId;
		int numDiscordReads, numDiscordReadsUnique;
		//Exons predicted to be involved 
		vector<Exon *> gAexons, gBexons;
		vector<DiscordantClusterAlignment> alignments;
		unordered_map<int, int> exonicReadCount;
		vector<Junction*> junctions;

		DiscordantCluster(Gene *gA_, Gene *gB_) :gA(gA_), gB(gB_)
		{
			clusterId = DiscordantClusterFinder::NextClusterID++;
			numDiscordReads = 0, numDiscordReadsUnique = 0;
		}

		int GetNumDiscordReads(bool uniqueOnly = false);

		bool operator < (const DiscordantCluster& c) const
		{
			return (alignments.size() > c.alignments.size());
		}

		static bool compareClusters(DiscordantCluster* a, DiscordantCluster* b)
		{
			return (a->alignments.size() > b->alignments.size()); 
		}

	};
}
#endif
