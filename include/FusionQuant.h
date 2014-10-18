#ifndef FUSION_QUANT_H
#define FUSION_QUANT_H

#pragma once

#include "Utils.h"
#include "Config.h"
#include "GeneModel.h"
#include "FusionJunction.h"
#include "Read.h"
#include "DiscordantReadFinder.h"

#include "BamReader.h"
#include "BamAlignment.h"

#define PAD_LENGTH 500 

using namespace BamTools;

namespace MOJO
{
	typedef std::pair < BamAlignment, BamAlignment > BamAlnPair;

	class JunctionStats
	{
		public:
			int numConcordsA_5p, numConcordsA_3p;
			int numConcordsB_5p, numConcordsB_3p;
			int numConcordsAA_span, numConcordsBB_span;
			int numConcordsAA_junct, numConcordsBB_junct;
			int numDiscords_A5p_B3p, numDiscords_A3p_B5p;
			int numDiscords_A5p_B5p, numDiscords_A3p_B3p;
			
			double rpkmA, rpkmB;
			double rpkmA_5p, rpkmA_3p, rpkmB_5p, rpkmB_3p;
			JunctionStats() 
			{
				numConcordsA_5p = 0, numConcordsA_3p = 0;
				numConcordsB_5p = 0, numConcordsB_3p = 0;
				numConcordsAA_span = 0, numConcordsBB_span = 0;
				numConcordsAA_junct = 0, numConcordsBB_junct = 0;
				numDiscords_A5p_B3p = 0, numDiscords_A3p_B5p = 0;
				numDiscords_A5p_B5p = 0, numDiscords_A3p_B3p = 0;
				
				rpkmA_5p = 0, rpkmA_3p = 0, rpkmB_5p = 0, rpkmB_3p = 0;
				rpkmA = 0, rpkmB = 0;
			}
	};

	struct BamAlignmentEnd 
	{
		int position;
		int alignedLength;
		int readId;
		int end;
	};

	class PairedBamAlignment
	{
		public:
			BamAlignmentEnd fAln, sAln;
			Exon *fExon, *sExon;
			Isoform *fIso, *sIso;
			
			PairedBamAlignment(BamAlignment a1, BamAlignment a2) ;

			int GetReadId() 
			{
				return fAln.readId;
			}
	};

	class FusionQuant
	{
		private:
			unordered_map<string, vector<PairedBamAlignment *> > geneAlnsConcordant;
			unordered_map<string, vector<PairedBamAlignment *> > geneAlnsDiscordant;
			string paddingSeq;

			bool DoesExonMapToJunction(Isoform *iso, Exon *ex,
				BamAlignmentEnd aln, FRAG_TYPE ftype);

			string GetPaddingSeq();

		public:
			static FusionQuant FQ;

			FusionQuant();

			void QuantifyMappingStats(vector<Junction*> junctions);

			JunctionStats CompileStatsForJunction(Junction *j);

			static void FindFusionGeneMappingReads_worker(int threadId, 
				ComputePerTask cpt, string mapFa);

			Exon* FindExonsMappedByRead(BamAlignment aln, string name);
	};
}

#endif