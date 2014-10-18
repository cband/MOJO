
#include "Read.h"
#include "GeneModel.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = "Main"));
	
	void Read::GetFastqString(char *buf, int rtrim, int ltrim, 
		bool convertToSanger)
	{
		if (ReadName.length() == 0) {
			sprintf(buf, "\0");
			return;
		}

		auto full_read_len = Sequence.length();
		auto len = full_read_len - (ltrim + rtrim) + 1;
		if (len < 0) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Too much trimming occurred, readname: " << ReadName 
				<< ", sequence: " << Sequence << ", ltrim: " << ltrim
				<< ", rtrim: " << rtrim << ", seqlength: " << Sequence.length();
			exit(1);
		}

		const char *orig_seq = Sequence.c_str(), *orig_qual = Quality.c_str();
		char *seq = new char[len], *qual = new char[len];
		for (string::size_type i = 0; i < full_read_len; i++)
			if (i >= ltrim && i < (full_read_len - rtrim)) {
			seq[i - ltrim] = orig_seq[i];
			if (convertToSanger)
				qual[i - ltrim] = (char)((int)orig_qual[i] - 31);
			else
				qual[i - ltrim] = orig_qual[i];
			}

		seq[len - 1] = '\0';
		qual[len - 1] = '\0';

		sprintf(buf, "@%s\n%s\n+%s\n%s\n\0", ReadName.c_str(), seq,
			ReadName.c_str(), qual);

		delete[] seq;
		delete[] qual;
	}

	void Read::GetFastqSubstring(char *buf, int maxLength, int startAt)
	{
		auto read_length = Sequence.length();
		string::size_type ltrim, rtrim;
		if (maxLength > read_length)
			ltrim = 0, rtrim = 0;
		else if ((startAt + maxLength) > read_length)
			ltrim = read_length - maxLength, rtrim = 0;
		else
			ltrim = startAt, rtrim = read_length - maxLength - startAt;
		GetFastqString(buf, rtrim, ltrim);
	}

	string Read::GetTrimmedReadName()
	{
		return Read::TrimReadName(ReadName);
	}

	//Trim reads names to remove the machine name and the read tail
	string Read::TrimReadName(string read, string findChr)
	{
		auto idx = read.find(findChr);
		if (idx >= 0)
			read = read.substr(idx + 1, read.length() - idx - 1);
		if ((int)read.find("/1") > 0)
			return read.substr(0, read.find("/1", 0));
		if ((int)read.find("/2") > 0)
			return read.substr(0, read.find("/2", 0));
		return read;
	}

	bool Read::AddAlignment(ReadAlignment newAln)
	{
		Gene *newG = ((Gene *)newAln.gene);
		for (auto al : alignments)
			if (((Gene *)al.gene) == newG)
				return false;
		alignments.push_back(newAln);
	}

	void Read::Initialize() 
	{
		ReadName = "";
		Sequence = "";
		Quality = "";
		alignments.clear();
	}

	void PairedRead::Initialize()
	{
		FirstRead.Initialize();
		SecondRead.Initialize();
	}
}