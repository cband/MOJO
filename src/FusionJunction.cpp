
#include "FusionJunction.h"
#include "DiscordantClusterFinder.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = ""));

	int AnchorRead::nextAnchorReadId = 1000000;

	string AnchorRead::GetSplitReadIdentifier()
	{
		return (junction->GetJunctionName() + "|" + splitRead.ReadName + "_"
			+ boost::lexical_cast<string> (splitRead.overhang5p));
	}

	void AnchorRead::ParseFromFaReadname(string line)
	{

	}

	//Determines if two anchor reads are PCR duplicates by comparing the 
	//sequence of the first 36bps. If there are 1 or fewer mismatches between 
	//the two split-ends or the two other-ends, the two anchor reads
	//are classified as PCR duplicates
	bool AnchorRead::AreAnchorReadsPCRdups(AnchorRead *a, AnchorRead *b)
	{
		int mmS = 0, mmO = 0;
		string aSplit = a->splitRead.Sequence, aOther = a->otherRead.Sequence;
		string bSplit = b->splitRead.Sequence, bOther = b->otherRead.Sequence;
		
		for (string::size_type i = 0; i < min(35, (int)aSplit.size()); i++) 
			if (aSplit[i] != bSplit[i])
				mmS++;

		if (mmS > 1)
			return false;

		for (string::size_type i = 0; i < min(35, (int)aOther.size()); i++)
			if (aOther[i] != bOther[i])
				mmO++;
		
		if (mmO > 1)
			return false;

		return true;

		//// DEPRECATED --> ////
		// The following code attempts to check by position; 
		//both ends of a/b can have multiple alignments;
		bool split_pcr_dup = false, other_pcr_dup = false;
		//check split reads
		for (auto a_al : a->splitRead.alignments){
			for (auto b_al : b->splitRead.alignments) {
				if ((Gene*)(a_al.gene) != (Gene*)(b_al.gene))
					continue;
				if (a_al.position == b_al.position)
					split_pcr_dup = true;
			}
		}
		for (auto a_al : a->otherRead.alignments){
			for (auto b_al : b->otherRead.alignments) {
				if ((Gene*)(a_al.gene) != (Gene*)(b_al.gene))
					continue;
				if (a_al.position == b_al.position)
					other_pcr_dup = true;
			}
		}
		if (split_pcr_dup && other_pcr_dup)
			return true;
		return false;
	}

	bool AnchorRead::IsAnchorRead()
	{
		//Minimum requirement of anchor read is that a split read maps to the
		//junction and the OtherRead maps to one of the two genes of the fusion
		//The other end should map to one of the two genes of the fusion pair
		Gene* g5p = junction->ex5p->gene, *g3p = junction->ex3p->gene;
		for (auto otherAln : otherRead.alignments) {
			Gene *g = ((Gene *)otherAln.gene);
			if (g == g5p || g == g3p || Gene::DoGenesOverlap(g, g5p) ||
				Gene::DoGenesOverlap(g, g3p))
				return true;
		}
		return false;
	}

	bool AnchorRead::IsHighConfidenceAR()
	{
		if (isSpurious)
			return false;
		
		if (!IsAnchorRead())
			return false;
		
		for (auto aln : otherRead.alignments) {
			//TODO:: Allow this as high confidence in special cases. eg: there could
			//be 3 or few pseudo-genes that share this same sequence.  Allow this under
			//those circumstances; ??
			Gene *gA = junction->ex5p->gene, *gB = junction->ex3p->gene;
			Gene *o = (Gene*)aln.gene;
			if (o != gA && o != gB) {
				// now check if there is an overlap in genomic region between
				// o and A or B
				if (Gene::DoGenesOverlap(o, gA) || Gene::DoGenesOverlap(o, gB))
					continue;
				return false;
			}
		}
		
		/*
		A confidence AnchorRead is required satisfy one of the critera below:
			-- (a) a proper anchor read: the other-read maps to corresponding pair of
			   the gene to which split-read's major overhang maps to
			-- (b) if (a) is not satisified, then consider as 'proper' if 
				minOverhang is > 20bps.
		*/
		if (splitRead.GetMinOverhang() >= 20 &&
			splitRead.mismatches < std::round(0.02 * splitRead.alignedLength) )
		{
			return true;
		}

		//Require the other-end map to oneo fo the two genes involved in the junction
		if (splitRead.overhang5p > splitRead.overhang3p) {
			for (auto al : otherRead.alignments)
				if (junction->ex3p->gene == ((Gene*)al.gene))
					return true;
		}
		else {
			for (auto al : otherRead.alignments)
				if (junction->ex5p->gene == ((Gene*)al.gene))
					return true;
		}
		return false;
	}

	string Junction::GetJunctionName()
	{
		stringstream ss;
		ss << clusterId << "_" << ex5p->exonId 
			<< "_" << ex3p->exonId << "_" << clusterOrientFlag 
			<< "_" << junctionId << "_" << breakpoint;
		return ss.str();
	}

	bool Junction::IsSpurious()
	{
		for (auto ar : anchorReads)
			if (!ar->isSpurious)
				return false;
		return true;
	}

	void Junction::SetAsSpurious()
	{
		for (auto ar : anchorReads)
			ar->isSpurious = true;
	}

	void Junction::MarkARspurious(string trimmedReadName)
	{
		for (auto ar : anchorReads)
			if (ar->splitRead.GetTrimmedReadName() == trimmedReadName)
				ar->isSpurious = true;
	}

	bool Junction::operator > (const Junction& j) const
	{
		return false;
	}

	// In some cases, an anchor read can be aligned to two different junctions. Use
	// the following criteria to compare the two junctions and determine which takes
	// precedence. 
	bool Junction::compareJunctions(Junction* a, Junction* b)
	{
		typedef DiscordantCluster DC;
		DC *ac = (DC*)a->clusterPtr;
		DC *bc = (DC*)b->clusterPtr;

		if (ac->GetNumDiscordReads(true) != bc->GetNumDiscordReads(true))
			return ac->GetNumDiscordReads(true) > bc->GetNumDiscordReads(true);

		if (ac->GetNumDiscordReads() != bc->GetNumDiscordReads())
			return ac->GetNumDiscordReads() > bc->GetNumDiscordReads();

		if (a->GetTotalARCount(true) != b->GetTotalARCount(true))
			return a->GetTotalARCount(true) > b->GetTotalARCount(true);

		if (a->GetTotalARCount() != b->GetTotalARCount())
			return a->GetTotalARCount() > b->GetTotalARCount();

		//Perhaps find additional criteria to distinguish the two...;
		
		//finally resort to lexical sorting;
		if (a->ex5p != b->ex5p)
			return a->ex5p->exonId > b->ex5p->exonId;
		if (a->ex3p != b->ex3p)
			return a->ex3p->exonId > b->ex3p->exonId;

		// Should not get here 
		return false;
	}

	bool Junction::ConstructJunctionFromName(string junctName, Junction *j)
	{
		GeneModel *gm = GeneModel::GetGeneModel();
		auto sp = Utils::SplitToVector(junctName, "_");
		j->clusterId = stoi(sp[0]);
		int ex5p = stoi(sp[1]), ex3p = stoi(sp[2]);
		
		if (gm->ExonsMap.find(ex5p) == gm->ExonsMap.end() ||
			gm->ExonsMap.find(ex3p) == gm->ExonsMap.end())
			return false;

		j->ex5p = gm->ExonsMap[ex5p], j->ex3p = gm->ExonsMap[ex3p];
		j->clusterOrientFlag = sp[3];
		j->junctionId = stoi(sp[4]);
		j->breakpoint = stoi(sp[5]);
		return true;
	}

	void Junction::MarkPCRduplicates()
	{
		//Sort the anchorReads in the descending order of the junction region overlap;
		sort(anchorReads.begin(), anchorReads.end(), AnchorRead::BestAnchorRead);

		typedef vector<AnchorRead *>::size_type size_t_ar;
		//Logic: to mark a PCR duplicate, the alignment positions of split/other 
		//ends of two PE reads should match
		for (size_t_ar i = 0; i < anchorReads.size(); i++)
		{
			if (anchorReads[i]->isSpurious)
				continue;
			for (size_t_ar j = i + 1; j < anchorReads.size(); j++)
				if (AnchorRead::AreAnchorReadsPCRdups(anchorReads[i], anchorReads[j]))
						anchorReads[j]->isSpurious = true;
		}
	}

	int Junction::GetTotalSRCount()
	{
		int count = 0;
		for (auto ar : anchorReads){
			if (ar->isSpurious)
				continue;
			count++;
		}
		return count;
	}

	// Get total number of anchor reads.  If "highConfOnly", then only
	// counts the high confidence anchor reads.  If "mainAnchorLen", only
	// counts the anchor reads with a minimum anchor length
	int Junction::GetTotalARCount(bool highConfOnly, int minAnchorLen )
	{
		int count = 0;
		for (auto ar : anchorReads){
			if (ar->isSpurious || !ar->IsAnchorRead())
				continue;
			if (highConfOnly && !ar->IsHighConfidenceAR())
				continue;
			if (ar->splitRead.GetMinOverhang() >= minAnchorLen)
				count++;
		}
		return count;
	}

	// Get total high confidence anchor reads with one end mapping to the junction 
	// and the other end mapping to the gene indicated by the 'fiveP' flag
	int Junction::GetTotalARCountUniqueBy53(bool fiveP)
	{
		int fivePcount = 0, threePcount = 0;
		for (auto ar : anchorReads) {
			if (!ar->IsHighConfidenceAR() || ar->otherRead.alignments.size() != 1)
				continue;
			if (ex5p->gene == ((Gene*)ar->otherRead.alignments[0].gene))
				fivePcount++;
			if (ex3p->gene == ((Gene*)ar->otherRead.alignments[0].gene))
				threePcount++;
		}
		return (fiveP ? fivePcount : threePcount);
	}

	vector<FusionFragment *> Junction::GetFusionFragments(FRAG_TYPE fragType)
	{
		Exon *e = (fragType == FRAG_TYPE::FIVE_P ? ex5p : ex3p);
		vector<FusionFragment *> fragments;
		for (auto iso : e->gene->isoforms)
			if (iso->HasExon(e))
				fragments.push_back(new FusionFragment(e, iso, fragType));
		return fragments;
	}

	vector<FusionTranscript *> Junction::GetFusionTranscripts()
	{
		if (fusionTranscripts.size() != 0)
			return fusionTranscripts;

		typedef vector<FusionFragment *> FFvect;
		typedef FusionFragment::BreakPos BP;
		FFvect fivePfrags = GetFusionFragments(FRAG_TYPE::FIVE_P);
		FFvect threePfrags = GetFusionFragments(FRAG_TYPE::THREE_P);

		boost::unordered_map<int, FusionTranscript *> dups;
		for (auto frag5 : fivePfrags) {
			for (auto frag3 : threePfrags) {
				int len = frag5->GetTranscribedSeq().length() +
					frag3->GetTranscribedSeq().length();

				auto ft = new FusionTranscript(frag5, frag3);
				if (dups.find(len) != dups.end()) {
					//An exon that is shared by multiple isoforms can be part of
					//translated sequence or a UTR; give preference to CDS
					auto ftPrev = dups[len];
					int scorePrev = (int)ftPrev->fiveP->breakPosition +
						(int)ftPrev->threeP->breakPosition;
					int score = (int)ft->fiveP->breakPosition +
						(int)ft->threeP->breakPosition;

					if (score < scorePrev)
						continue;
				}
				dups[len] = ft;
				fusionTranscripts.push_back(ft);
			}
		}
		return fusionTranscripts;
	}

	bool FusionFragment::IsCdsStartCovered()
	{
		return false;
	}

	bool FusionFragment::IsCdsEndCovered()
	{
		return false;
	}

	string FusionFragment::GetCodingSeq()
	{
		if (codingSeq == "")
			ConstructSeq();
		return codingSeq;
	}

	string FusionFragment::GetTranscribedSeq()
	{
		if (transcribedSeq == "")
			ConstructSeq();
		return transcribedSeq;
	}

	string FusionFragment::GetStringForBreakPos(BreakPos pos)
	{
		switch (pos)
		{
		case UNKNOWN:
			return "Unknown";
		case UTR5:
			return "UTR-5";
		case UTR3:
			return "UTR-3";
		case CDS:
			return "CDS";
		case NONCODING:
			return "ncRNA";
		default:
			return "UNKNOWN";
		}
	}

	void FusionFragment::ConstructSeq()
	{
		string fragment;
		string t; //transcribedSeq;
		string c; //codingSeq;

		bool fExonFound = false;
		for (auto ex : fIsoform->exons){
			string exSeq = (ex->gene->strand == "+" ? ex->sequence : 
				Utils::ReverseComplement(ex->sequence));
			fragment += exSeq;
			if (ex == fExon) {
				if (ex->gene->strand == "+" && fragType == FIVE_P)
					t = fragment;
				else if (ex->gene->strand == "-" && fragType == THREE_P)
					t = Utils::ReverseComplement(fragment);
				fExonFound = true;
				fragment = exSeq;
			}
			if (fExonFound && ex == fIsoform->exons[fIsoform->exons.size() - 1]) {
				if (ex->gene->strand == "+" && fragType == THREE_P)
					t = fragment;
				else if (ex->gene->strand == "-" && fragType == FIVE_P)
					t = Utils::ReverseComplement(fragment);
			}
		}

		transcribedSeq = t;
		int len = fIsoform->GetTranscriptLength();
		int st = fIsoform->cdsStartDistFromTxStart;
		int en = fIsoform->cdsEndDistFromTxEnd;
		
		//check if the gene is non-coding (cdsStart/cdsEnd are -1 in the MasterIsoform)
		if (st < 0 && en < 0) {
			breakPosition = BreakPos::NONCODING;
			codingSeq = "";
			return;
		}

		if (fragType == FIVE_P) {
			if (t.length() <= st) {
				breakPosition = BreakPos::UTR5;
				return;
			}
			else if (t.length() < (len - en)) {
				breakPosition = BreakPos::CDS;
				c = t.substr(st, t.length() - st);
			}
			else {
				breakPosition = BreakPos::UTR3;
				c = t.substr(st, t.length() - st - en);
			}
		}
		else {
			if (t.length() <= en) {
				breakPosition = BreakPos::UTR3;
				return;
			}
			else if (t.length() < (len - st)) {
				breakPosition = BreakPos::CDS;
				c = t.substr(0, t.length() - en);
			}
			else {
				breakPosition = BreakPos::UTR5;
				c = t.substr(st - (len - t.length()), len - st - en);
			}
		}
		codingSeq = c;
	}

	void FusionTranscript::GetCodonPositions(string sequence, string codon,
		vector<int> *sites)
	{
		auto pos = sequence.find(codon);
		while (pos != string::npos && pos >= 0) {
			if (std::find((*sites).begin(), (*sites).end(), pos) == (*sites).end())
				(*sites).push_back(pos);
			pos = sequence.find(codon, pos + 1);
		}
	}

	bool FusionTranscript::IsInFrame()
	{
		/// >>>>>>>>>>>>>>>>>>>>>>>>>
		// Inframe should include an annotated 3'stop codon for proper translation
		// termination
		/// <<<<<<<<<<<<<<<<<<<<<<<<<
		if (fiveP->breakPosition == FusionFragment::BreakPos::UTR3 ||
			threeP->breakPosition == FusionFragment::BreakPos::UTR5)
			return true;	

		if (fiveP->breakPosition == FusionFragment::BreakPos::CDS && 
			threeP->breakPosition == FusionFragment::BreakPos::CDS)
		{
			if ((fiveP->GetCodingSeq().length() + 
				threeP->GetCodingSeq().length()) % 3 == 0)
				return true;
		}

		return false;
	}

	string FusionTranscript::GetCodingSequence()
	{
		if (!IsInFrame())
			return "";
		return fiveP->GetCodingSeq() + threeP->GetCodingSeq();
	}

	string FusionTranscript::GetTranscriptSequence()
	{
		return fiveP->GetTranscribedSeq() + threeP->GetTranscribedSeq();
	}
}