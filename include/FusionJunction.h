#ifndef FUSION_JUNCTION_H
#define FUSION_JUNCTION_H

#pragma once

#include "Utils.h"
#include "GeneModel.h"
#include "FastqParser.h"
#include "Read.h"

namespace MOJO 
{
	class Junction;
	class FusionFragment;
	class FusionTranscript;

	struct Expression
	{
		string name;
		double value;
		double percentileRank;
	};

	class SplitRead : public Read 
	{
		public:
			string cigar, strand;
			string alignedSequence; // this is independent of the Read.Sequence;
									// only represents the sequence aligned to junction
			int flag, position, alignedLength, sequenceLength;
			int overhang5p, overhang3p;
			int mismatches, delBefore, delAfter;
			bool mapsToScrambledExons, mismatchInAnchor;

			SplitRead() 
			{
				flag = 0, position = -1, alignedLength = -1, sequenceLength = -1;
				overhang5p = -1, overhang3p = -1;
				delAfter = -1, delBefore = -1;
				mapsToScrambledExons = false, mismatchInAnchor = false;
			}

			SplitRead& operator=(const Read& rd) 
			{
				(*this).ReadName = rd.ReadName;
				(*this).Sequence = rd.Sequence;
				(*this).Quality = rd.Quality;
				return *this;
			}

			string GetAlignedSeq()
			{
				return Sequence.substr(delBefore, 
					Sequence.length() - delBefore - delAfter);
			}

			int GetMinOverhang() 
			{
				if (min(overhang5p, overhang3p) == -1) {
					cerr << "Anchor reads not properly loaded"; 
					exit(1);
				}
				return std::min<int>(overhang5p, overhang3p);
			}
	};

	class OtherRead : public Read 
	{
		public:
			OtherRead& operator=(const Read& rd) 
			{
				(*this).ReadName = rd.ReadName;
				(*this).Sequence = rd.Sequence;
				(*this).Quality = rd.Quality;
				return *this;
			}
	};

	class AnchorRead 
	{
		private:
			static int nextAnchorReadId;
		
		public:
			Junction *junction;
			SplitRead splitRead;
			OtherRead otherRead;
			bool isSpurious;
			int anchorReadId;

			AnchorRead(Junction *junct) : isSpurious(false), junction(junct) 
			{
				anchorReadId = GetNextId();
			}

			string GetSplitReadIdentifier();

			void ParseFromFaReadname(string line);

			bool IsHighConfidenceAR();

			bool IsAnchorRead();

			static bool BestAnchorRead(const AnchorRead *a, const AnchorRead *b) 
			{
				return ((a->splitRead.overhang5p + a->splitRead.overhang3p) 
						> (b->splitRead.overhang5p + b->splitRead.overhang3p));
			}

			static bool AreAnchorReadsPCRdups(AnchorRead *a, AnchorRead *b);

			static int GetNextId() 
			{
				return ++nextAnchorReadId;
			}
	};

	class Junction
	{
		private:
			vector<FusionTranscript *> fusionTranscripts;

			vector<FusionFragment *> GetFusionFragments(FRAG_TYPE frag_type);
		
		public:
			int clusterId;
			void *clusterPtr;
			Exon *ex5p, *ex3p;
			string clusterOrientFlag;
			int junctionId, breakpoint;

			Expression gA_TPM, gB_TPM, gAB_TPM;
			vector<Expression> isoA_TPM, isoB_TPM;

			char *c_seq5p, *c_seq3p;
			//string seq5p, seq3p;
			vector<AnchorRead *> anchorReads;
			
			Junction(){
				clusterPtr = 0;
			}

			Junction(int clustId, void *clustPtr, Exon *fiveP, Exon *threeP,
				string flag, int id, char *fivePseq, char *threePseq) : 
				clusterId(clustId), clusterPtr(clustPtr), ex5p(fiveP), 
				ex3p(threeP), clusterOrientFlag(flag), junctionId(id),
				c_seq5p(fivePseq), c_seq3p(threePseq)
			{
				breakpoint = (int) strlen(fivePseq);
			}

			string GetJunctionName();

			bool IsSpurious();

			void SetAsSpurious();

			void MarkARspurious(string trimmedReadName);

			void MarkPCRduplicates();

			int GetTotalSRCount();

			int GetTotalARCount(bool highConfOnly = false, int minAnchorLen = 0);
			
			int GetTotalARCountUniqueBy53(bool fiveP = true);

			static bool ConstructJunctionFromName(string junctName, Junction *j);

			vector<FusionTranscript *> GetFusionTranscripts();

			bool operator > (const Junction& j) const;
			
			bool operator==(const Junction& j) const
			{
				return this->ex5p->exonId == j.ex5p->exonId && 
					this->ex3p->exonId == j.ex3p->exonId;
			};

			static bool compareJunctions(Junction* a, Junction* b);
	};

	//FusionFragment represents the 5' or 3' fusion transcript.  A given fusion between
	//two genes may have multiple 5' and 3' isoforms
	class FusionFragment
	{
		private:
			bool isCdsStartCovered, isCdsEndCovered;
			string codingSeq, transcribedSeq;

			void ConstructSeq();
		
		public:
			enum BreakPos { 
				UNKNOWN = 0, 
				NONCODING, 
				UTR3, 
				UTR5, 
				CDS = 8 
			};

			Isoform *fIsoform;
			Exon *fExon;
			FRAG_TYPE fragType;
			BreakPos breakPosition;

			FusionFragment() {};

			FusionFragment(Exon *ex, Isoform *iso, FRAG_TYPE ftype) : fExon(ex),
				fIsoform(iso), fragType(ftype), breakPosition(UNKNOWN) {};

			bool IsCdsStartCovered();
			
			bool IsCdsEndCovered();
			
			string GetCodingSeq();
			
			string GetTranscribedSeq();

			static string GetStringForBreakPos(BreakPos pos);
	};

	// A FusionTranscript is the full length fusion isoform as a pairing of two 
	// FusionFragments.  A fusion junction does not provide evidence as to which
	// isoform the predicted fusion is expected to generate.  Therefore, all isoforms
	// are considered
	class FusionTranscript 
	{
		private:
			void GetCodonPositions(string sequence, string codon, vector<int> *sites);
		public:
			FusionFragment *fiveP, *threeP;
			Expression expression;

			FusionTranscript() { }

			FusionTranscript(FusionFragment *five, FusionFragment *three) 
			{
				fiveP = five, threeP = three;
			}

			bool IsInFrame();
			
			string GetCodingSequence();
			
			string GetTranscriptSequence();
	};
};
#endif
