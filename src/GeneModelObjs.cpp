
#include "GeneModelObjs.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = ""));

	// Computes the distance between the two coordinates of a gene that is 
	// spanned only by exons.  Considers all possible isoforms.
	int Gene::GetExonicDistanceBetweenCoords(int coordA, int coordB)
	{
		if (coordA > coordB)
			std::swap(coordA, coordB);

		if (coordA > coordB) {
			cerr << "coordA is greater than coordB, error." << endl;
			exit(1);
		}

		bool coordAexonic = false, coordBexonic = false;
		//allow mapping coordinate to be at most 20bp from the exon end;
		const int PAD = 20; 

		int distance = 0;
		for (auto ex : allExons){
			//if exon between two coords
			if (coordA < ex->exStart && coordB > ex->exEnd) { 
				distance += ex->GetExonLength();
			}
			else if ((ex->exStart - PAD) < coordA && (ex->exEnd + PAD) > coordA) {
				//else if exon has coordA
				distance += abs(ex->exEnd - coordA);
				coordAexonic = true;
			}
			else if ((ex->exStart - PAD) < coordB && (ex->exEnd + PAD) > coordB) {
				//else if exon has coordB
				distance += abs(coordB - ex->exStart);
				coordBexonic = true;
			}
		}
		if (coordAexonic && coordBexonic)
			return distance;
		return 1000000;
	}

	bool Gene::DoesRegionContainRepeat(int readStart, int readEnd)
	{
		for (auto repeat : repeatRegions)
			return (repeat.first < readStart && repeat.second > readEnd);
		return false;
	};

	// For a given gene, generates all possible back-spliced exon junctions to 
	// simulate the circular-RNA or intra-genic rearranements
	string Gene::GetScrambledExonSequence()
	{
		if (scrambledExonSequence.length() > 0)
			return scrambledExonSequence;

		stringstream ss, ss_scram;
		unordered_map<string, bool> seq_added;
		for (vector<Isoform*>::size_type i = 0; i < isoforms.size(); i++)
		{
			for (vector<Exon*>::size_type j = 0;
				j < isoforms[i]->exons.size(); j++)
			{
				Exon *exon1 = isoforms[i]->exons[j];
				string exon1Str = lexical_cast<string>(exon1->exonId);
				if (seq_added.find(exon1Str) == seq_added.end()) {
					string seq = exon1->sequence;
					string seqA = seq, seqB = seq;
					if (seq.length() > 90) {
						seqA = seq.substr(seq.length() - 90, 90);
						seqB = seq.substr(0, 90);
					}
					ss << seqA << seqB;
					seq_added[exon1Str] = true;
				}
				for (vector<Exon*>::size_type k = j + 1;
					k < isoforms[i]->exons.size(); k++)
				{
					Exon *exon2 = isoforms[i]->exons[k];
					string exon2Str = lexical_cast<string>(exon2->exonId);
					if (seq_added.find(exon1Str + exon2Str) != seq_added.end())
						continue;
					string s1, s2;
					if (strand == "+")
						s1 = exon1->sequence, s2 = exon2->sequence;
					else
						s1 = exon2->sequence, s2 = exon1->sequence;

					string s1_scram = s1.length() > 90 ? s1.substr(0, 90) : s1;
					string s2_scram =
						s2.length() > 90 ? s2.substr(s2.length() - 90, 90) : s2;
					string s1_norm =
						s1.length() > 90 ? s1.substr(s1.length() - 90, 90) : s1;
					string s2_norm = s2.length() > 90 ? s2.substr(0, 90) : s2;
					seq_added[exon1Str + exon2Str] = true;
					ss_scram << s2_scram << s1_scram;
					ss << s1_norm << s2_norm;
				}
			}
		}
		scrambledExonSequence = ss.str() + ss_scram.str();
		return scrambledExonSequence;
	}

	string Gene::GetTranscribedSequence()
	{
		std::map<int, char> seqMap;
		for (auto e : allExons){
			const char *s = e->sequence.c_str();
			for (int i = e->exStart; i < e->exEnd; i++)
				seqMap[i] = s[i - e->exStart];
		}
		char *cc = new char[seqMap.size() + 1];
		int i = 0;
		for (auto c : seqMap)
			cc[i++] = c.second; 
		cc[i] = '\0';
		string fullSeq(cc);
		delete[] cc;
		if (strand == "+")
			return fullSeq;
		else
			return Utils::ReverseComplement(fullSeq);
	}

	// Generate the 5' or 3' transcript ending/beginning with the targetExon.
	string Gene::GetPartialTranscript(Exon *targetExon, FRAG_TYPE frag) {
		string seq;
		std::map<int, char> upstreamMap, downstreamMap;
		bool found = false;
		for (auto e : allExons){
			const char *s = e->sequence.c_str();

			for (int i = e->exStart; i < e->exEnd; i++) {
				if (!found)
					upstreamMap[i] = s[i - e->exStart];
				else
					downstreamMap[i] = s[i - e->exStart];

				if (e == targetExon) {
					found = true;
					upstreamMap[i] = s[i - e->exStart];
					downstreamMap[i] = s[i - e->exStart];
				}
			}
		}

		std::map<int, char> *seqMap;
		if ((frag == FRAG_TYPE::FIVE_P && strand == "+") ||
			(frag == FRAG_TYPE::THREE_P && strand == "-"))
			seqMap = &upstreamMap;
		else if ((frag == FRAG_TYPE::FIVE_P && strand == "-") ||
			(frag == FRAG_TYPE::THREE_P && strand == "+"))
			seqMap = &downstreamMap;

		char *cc = new char[(*seqMap).size() + 1];
		int i = 0;
		for (auto c : (*seqMap)) 
			cc[i++] = c.second;
		cc[i] = '\0';
		string fullSeq(cc);
		delete[] cc;
		if (strand == "+")
			return fullSeq;
		else
			return Utils::ReverseComplement(fullSeq);

		return seq;
	}

	//Returns all exons spanned by the coordinates start/end;
	vector<Exon *> Gene::FindExonsWithinPoints(int start, int end,
		bool considerIntronic)
	{
		vector<Exon *> exons;
		for (auto iso : isoforms) {
			Exon *prevExon = 0;
			for (auto ex : iso->exons) {
				if ((ex->exStart <= start && ex->exEnd >= start) ||
					(ex->exStart <= end && ex->exEnd >= end) ||
					(considerIntronic && prevExon != 0 &&
					prevExon->exEnd <= start && ex->exStart >= end))
				{
					exons.push_back(ex);
				}
				prevExon = ex;
			}
		}
		return exons;
	}

	bool Gene::IsReadThru(Gene *left, Gene *right, int dist)
	{
		if (right->chr != left->chr || right->strand != left->strand)
			return false;

		if ((left->strand == "+" && left->txStart_genomic < right->txStart_genomic)
			||
			(left->strand == "-" && left->txStart_genomic > right->txStart_genomic))
		{
			return Gene::DistanceBetweenGenes(left, right) < dist;
		}
		return false;
	}

	int Gene::DistanceBetweenGenes(Gene *gA, Gene *gB) {
		if (gA->chr != gB->chr)
			return -1;

		if (gA->strand == gB->strand) {
			return min(abs(gA->txEnd_genomic - gB->txStart_genomic),
				abs(gA->txStart_genomic - gB->txEnd_genomic));
		}
		else {
			return min(abs(gA->txStart_genomic - gB->txStart_genomic),
				abs(gA->txEnd_genomic - gB->txEnd_genomic));
		}
	};

	bool Gene::DoGenesOverlap(Gene *gA, Gene *gB) {
		if (gA->chr != gB->chr)
			return false;

		int gA_st = gA->txStart_genomic, gA_en = gA->txEnd_genomic;
		int gB_st = gB->txStart_genomic, gB_en = gB->txEnd_genomic;

		if ((gA_st <= gB_st && gA_en >= gB_st) ||
			(gA_st <= gB_en && gA_en >= gB_en) ||
			(gB_st <= gA_st && gB_en >= gA_st) ||
			(gB_st <= gA_en && gB_en >= gA_en))
		{
			return true;
		}
		return false;
	}

	bool Isoform::HasExon(Exon *e)
	{
		for (auto ex : exons)
			if (ex == e)
				return true;
		return false;
	}

	string Isoform::GetTranscript()
	{
		stringstream ss;
		for (auto ex : exons){
			if (GetStrand() == "+")
				ss << ex->sequence;
			else
				ss << Utils::ReverseComplement(ex->sequence);
		}
		if (GetStrand() == "+")
			return ss.str();
		else
			return Utils::ReverseComplement(ss.str());
	}

	int Isoform::GetTranscriptLength() {
		int len = 0;
		for (auto ex : exons)
			len += abs(ex->exStart - ex->exEnd);
		return len;
	}

	//Get exon sequence for junction construction.  If the length of the exon is
	//less minLength, use upstream/downstream exons to generate a sequence of 
	//specified length. If alternate isoforms are involved, multiple sequences
	//are returned
	vector<string> Exon::GetJunctionSeqAs5pExon(int minLength)
	{
		if (sequence.length() >= minLength)
			return vector < string > {
			sequence.substr(sequence.length() - minLength, minLength)};

		unordered_map<string, bool> allPossibleSeqsMap;
		for (auto iso : gene->isoforms) {
			string seq = "";
			bool exonFound = false;
			for (auto exon : iso->exons) {
				seq = seq + exon->sequence;
				if (this == exon) {
					exonFound = true;
					break;
				}
			}
			if (!exonFound)
				continue;
			seq = seq.substr(max((int)seq.length() - minLength, 0), minLength);
			allPossibleSeqsMap[seq] = true;
		}

		vector<string> allPossible;
		for (auto iter : allPossibleSeqsMap)
			allPossible.push_back(iter.first);

		return allPossible;
	};

	//Get exon sequence for junction construction.  If the length of the exon is
	//less minLength, use upstream/downstream exons to generate a sequence of 
	//specified length. If alternate isoforms are involved, multiple sequences
	//are returned
	vector<string> Exon::GetJunctionSeqAs3pExon(int minLength)
	{
		if (sequence.length() >= minLength)
			return vector < string > {sequence.substr(0, minLength)};

		unordered_map<string, bool> allPossibleSeqsMap;
		for (auto iso : gene->isoforms) {
			string seq = "";
			bool exonFound = false;
			for (auto exon : iso->exons) {
				if (this == exon)
					exonFound = true;
				if (exonFound)
					seq = seq + exon->sequence;
			}
			if (!exonFound)
				continue;
			allPossibleSeqsMap[seq.substr(0, minLength)] = true;
		}

		vector<string> allPossible;
		for (auto iter : allPossibleSeqsMap)
			allPossible.push_back(iter.first);

		return allPossible;
	};

	CdsFragment::CdsFragment() {
		cdsStartIsIncluded = false;
		cdsEndIsIncluded = false;
		isoform = 0;
		breakPos = BreakPos::UNKNOWN;
	}

	CdsFragment::CdsFragment(string tx_seq, string cds_seq, bool start, bool end,
		Isoform *iso, BreakPos pos)
	{
		transcriptSeq = tx_seq;
		codingSeq = cds_seq;
		cdsStartIsIncluded = start, cdsEndIsIncluded = end;
		isoform = iso;
		breakPos = pos;
	}
};