#ifndef GENE_MODEL_OBJS_H
#define GENE_MODEL_OBJS_H

#pragma once

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include "Utils.h"

using namespace std;
using namespace boost;

namespace MOJO 
{
	class Gene;
	class Isoform;
	class Exon;
	class CdsFragment;

	enum FRAG_TYPE { 
		FIVE_P, 
		THREE_P 
	};

	struct CoordinatePair
	{
		public:
			int startA, endA, startB, endB;

			CoordinatePair() : startA(0), endA(0), startB(0), endB(0) {}

			CoordinatePair(int sA, int eA, int sB, int eB) :
				startA(sA), endA(eA), startB(sB), endB(eB) {}
	};

	class Gene
	{
		private:
			string scrambledExonSequence;

		public:
			string name, chr, strand;
			int geneId;
			int txStart_genomic, txEnd_genomic;
			int txStart_transcriptomic, txEnd_transcriptomic;
			vector<string> isoformsStrVect;
			vector<Isoform *> isoforms;
			vector<Exon *> allExons;

			boost::unordered_map<int, vector<CoordinatePair> > blastHits;
			vector< std::pair<int, int> > repeatRegions;

			// upstream/downstream genes
			int nextGeneDist, prevGeneDist;
			Gene *nextGene, *prevGene;

			Gene() {}

			Gene(string id, string nm, string chr_, string  strand_,
				string isoformsCsv)
			{
				geneId = stoi(id.substr(1, id.length() - 1));
				name = nm;
				chr = chr_, strand = strand_;
				isoformsStrVect = Utils::SplitToVector(isoformsCsv, ",");
				if (isoformsStrVect[isoformsStrVect.size() - 1] == "")
					isoformsStrVect.pop_back();
			}

			//Get the size of exonic region between two coordinates within a gene
			int GetExonicDistanceBetweenCoords(int coordA, int coordB);

			vector<Exon*> FindExonsWithinPoints(int start, int end,
				bool considerIntronic = true);

			bool DoesRegionContainRepeat(int readStart, int readEnd);

			bool operator==(const Gene & g) const
			{
				return this->geneId == g.geneId;
			};

			bool operator < (const Gene& g) const
			{
				if (chr == g.chr)
					return txStart_transcriptomic < g.txStart_transcriptomic;
				return (chr < g.chr);
			}

			string GetScrambledExonSequence();

			string GetTranscribedSequence();

			string GetPartialTranscript(Exon *targetExon, FRAG_TYPE frag);
				
			static bool compareGenes(Gene* a, Gene* b)
			{ 
				return (*a < *b); 
			}
			
			static bool IsReadThru(Gene *gA, Gene *gB, int dist);
			
			static int DistanceBetweenGenes(Gene *gA, Gene *gB);
			
			static bool DoGenesOverlap(Gene *gA, Gene *gB);
	};

	class Isoform
	{
		public:
			Gene *gene;
			string isoformIdStr, transcriptSequence;
			int txStart, txEnd, txStart_genomic, txEnd_genomic;
			int cdsStartDistFromTxStart, cdsEndDistFromTxEnd;
			vector<string> exonsStrVect;
			vector<Exon*> exons;

			string GetChr() 
			{ 
				return ((gene != 0) ? gene->chr : ""); 
			}
			
			string GetStrand() 
			{ 
				return ((gene != 0) ? gene->strand : ""); 
			}

			Isoform(){}
			
			Isoform(Gene *g, string id, string exonCsv)
			{
				gene = g;
				isoformIdStr = id;
				exonsStrVect = Utils::SplitToVector(exonCsv, ",");
				if (exonsStrVect[exonsStrVect.size() - 1] == "")
					exonsStrVect.pop_back();
			}
			
			bool HasExon(Exon *e);
			
			string GetTranscript();
			
			int GetTranscriptLength();
	};

	class Exon
	{
		public:
			Gene *gene;
			vector<Isoform *> exonIsoforms;
			int exonId;
			int exStart, exEnd, exStart_genomic, exEnd_genomic;
			string sequence;
			
			Exon() {}
			
			Exon(Gene *g, string id)
			{
				gene = g;
				exonId = stoi(id.substr(1, id.length() - 1));
				exStart = 0, exEnd = 0;
			}

			int GetExonLength()
			{
				return exEnd - exStart;
			}

			vector<string> GetJunctionSeqAs5pExon(int minLength = 80);
			
			vector<string> GetJunctionSeqAs3pExon(int minLength = 80);

			vector<CdsFragment> GetUpstreamDownstreamSequences(int pos);

			bool operator < (const Exon& ex) const 
			{
				return (exStart < ex.exStart);
			}

			bool operator == (const Exon& ex) const 
			{
				return (exonId == ex.exonId);
			}

			static bool compareExons(Exon* a, Exon* b) 
			{ 
				return (*a < *b); 
			}
			
			static bool compareExonsGenomic(Exon* a, Exon* b) 
			{ 
				return ((*a).exStart_genomic < (*b).exStart_genomic); 
			}
	};

	class CdsFragment
	{
		public:
			enum BreakPos {
				UNKNOWN = 0, 
				UTR5, 
				CDS, 
				UTR3 
			};

			string codingSeq, transcriptSeq;
			Isoform *isoform;
			bool cdsStartIsIncluded, cdsEndIsIncluded;
			BreakPos breakPos;
			
			CdsFragment();

			CdsFragment(string tx_seq, string cds_seq, bool start, bool end,
				Isoform *iso, BreakPos pos);
	};
}
#endif