
#include "FastqParser.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >,
		(keywords::channel = "Main"));

	// Gets the next read in the fastq file. The new read is returned
	// in the input argument.  If either the fastqParser is uninitialized or if 
	// the end of the file is reached, returns false;
	bool FastqFile::GetNextRead(Read *read)
	{
		if (fileStream != 0 && EndOfFile()) 
			return false;

		if (fileStream == 0) {
			if (!Utils::FileExists(fileName)) {
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "Error occured in FastqFile::GetNextRead(). File '"
					<< fileName << "' not found";
				exit(1);
			}
			fileStream = new boost::iostreams::filtering_istream();
			if (fileName.substr(fileName.length() - 2, 2) == "gz") {
				std::ifstream *ifs = new std::ifstream(fileName.c_str(),
					std::ios_base::in | std::ios_base::binary);
				fileStream->push(boost::iostreams::gzip_decompressor());
				fileStream->push((*ifs));
			}
			else {
				std::ifstream *ifs = new std::ifstream(fileName.c_str());
				fileStream->push((*ifs));
			}
			fileOpened = true;
		}
		string s;
		if (!std::getline((*fileStream), read->ReadName))
			return false;
		if (!std::getline((*fileStream), read->Sequence))
			return false;
		if (!std::getline((*fileStream), s))
			return false;
		if (!std::getline((*fileStream), read->Quality))
			return false;

		if ((read->ReadName.c_str())[0] == '@')
			read->ReadName = read->ReadName.substr(1, read->ReadName.length() - 1);
		return true;
	}

	//Return true if file is opened and is at eof.
	bool FastqFile::EndOfFile() 
	{
		return !(fileStream==0) && fileStream->eof();
	}

	//Get a paired-end read starting at a specific position in the index file
	bool FastqFile::GetReadStartingAtPos(long long pos, Read *read)
	{
		if (fileStreamByPos == 0) {
			if (!Utils::FileExists(fileName)) {
				BOOST_LOG_CHANNEL(logger::get(), "Main")
					<< "Error occured in FastqFile::GetNextRead(). File '" 
					<< fileName << "' not found";
				exit(1);
			}
			fileStreamByPos = new std::ifstream(fileName.c_str());
		}
		
		(*fileStreamByPos).seekg(pos, ios::beg);
		string s;
		std::getline((*fileStreamByPos), read->ReadName);
		std::getline((*fileStreamByPos), read->Sequence);
		std::getline((*fileStreamByPos), s);
		std::getline((*fileStreamByPos), read->Quality);
		if ((read->ReadName.c_str())[0] == '@')
			read->ReadName = read->ReadName.substr(1, read->ReadName.length() - 1);
		return true;
	}

	FastqParser::FastqParser(string firstFile, string secondFile) 
	{
		end1File = new FastqFile(firstFile);
		end2File = new FastqFile(secondFile);
		numReads = -1;
	}

	// Paired-end fastq parser.  Index file contains read id and the SEEK position 
	// in the index file.  Index file is optional.
	FastqParser::FastqParser(string firstFile, string secondFile, string indexFile)
	{
		if (firstFile.substr(firstFile.length() - 3, 3) == ".gz") {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Cannot use gz files to "
				<< "instantiate a FastqParser with indexing support. "
				<< "Exiting." << endl;
			exit(1); 
		}
		end1File = new FastqFile(firstFile);
		end2File = new FastqFile(secondFile);

		typedef fstream::pos_type fposType;
		ifstream fqIndFile(indexFile.c_str());
		for (string str; getline(fqIndFile, str);){
			auto sp = Utils::SplitToVector(str, "\t");
			long long pos1 = lexical_cast<long long>(sp[2]);
			long long pos2 = lexical_cast<long long>(sp[3]);
			std::pair<long long, long long> p(pos1, pos2);
			index[stoi(sp[1])] = p;
		}
	}

	//Retrieves read by ID from an indexed Fastq file
	PairedRead FastqParser::GetReadById(int id)
	{
		PairedRead pr; 
		typedef fstream::pos_type fposType;
		std::pair<long long, long long> pos = index[id];
		end1File->GetReadStartingAtPos(pos.first, &pr.FirstRead);
		end2File->GetReadStartingAtPos(pos.second, &pr.SecondRead);
		return pr;
	}

	// outputs reads of specific length to output files. 
	// startAt: indicates where to start the trimmed read
	// maxLength: indicates the maximum length of the read.  
	// Both parameters are adjusted based on the length of the actual read.
	void FastqParser::CreateTrimmedFiles(string firstEndOutFile, 
		string secondEndOutFile,  int maxLength, int startAt)
	{
		ofstream firstOut(firstEndOutFile.c_str(), ios::out | ios::app);
		ofstream secondOut(secondEndOutFile.c_str(), ios::out | ios::app);

		PairedRead pr;
		while (GetNextPairedRead(&pr)) {
			char firstBuf[2000], secondBuf[2000];
			auto read_length = pr.FirstRead.Sequence.length();
			string::size_type ltrim, rtrim;

			if (maxLength > read_length)
				ltrim = 0, rtrim = 0;
			else if ((startAt + maxLength) > read_length)
				ltrim = read_length - maxLength, rtrim = 0;
			else
				ltrim = startAt, rtrim = read_length - maxLength - startAt;
			pr.FirstRead.GetFastqString(firstBuf, rtrim, ltrim, false);

			read_length = pr.SecondRead.Sequence.length();
			if (maxLength > read_length)
				ltrim = 0, rtrim = 0;
			else if ((startAt + maxLength) > read_length)
				ltrim = read_length - maxLength, rtrim = 0;
			else
				ltrim = startAt, rtrim = read_length - maxLength - startAt;

			pr.SecondRead.GetFastqString(secondBuf, rtrim, ltrim, false);
			firstOut.write(firstBuf, strlen(firstBuf));
			secondOut.write(secondBuf, strlen(secondBuf));
		}
	}
	
	// This function is deprecated
	PairedRead FastqParser::GetNextPairedRead() 
	{
		//Expectation is that first and second end files exist;
		PairedRead pr;
		end1File->GetNextRead(&pr.FirstRead);
		end2File->GetNextRead(&pr.SecondRead);
		return pr;
	}

	// if the FastqParser is initialized, returns the next paired-end
	// read in the fastq file
	bool FastqParser::GetNextPairedRead(PairedRead *pr) 
	{
		//Expectation is that first and second end files exist;
		pr->Initialize();
		if (!end1File->GetNextRead(&pr->FirstRead))
			return false;
		if (!end2File->GetNextRead(&pr->SecondRead))
			return false;
		return true;
	}

	bool FastqParser::EndOfFile() 
	{
		if (end1File->EndOfFile() || end2File->EndOfFile())
			return true;
		return false;
	}

	FASTQ_ENCODING FastqParser::GetFastqEncoding() 
	{
		Read rd;
		int min = -1, max = -1;
		//SANGER: 33-73,  ILLUMINA-1.3: 64-104
		//SOLEXA: 59-104,	---> These are not supported now
		for (int line_cnt = 0; line_cnt < 5000 && end1File->GetNextRead(&rd); 
			line_cnt++) 
		{
			const char *qualities = rd.Quality.c_str();
			for (size_t i = 0; i < strlen(qualities); i++) {
				auto ascii_code = (int)qualities[i];
				if (min > ascii_code || min == -1) min = ascii_code;
				if (max < ascii_code || max == -1) max = ascii_code;
			}
		}
		if (min <= 58 && max <= 74) return FASTQ_ENCODING::SANGER;
//		if (min >= 64 && max <= 104) return FASTQ_ENCODING::ILLUMINA;

		// have to use 105 to account for a wierd sequencer configuration issue;
		if (min >= 63 && max <= 105) return FASTQ_ENCODING::ILLUMINA;  

		return FASTQ_ENCODING::UNSUPPORTED;
	}

	string FastqParser::GetFastqEncodingString(FASTQ_ENCODING encoding) 
	{
		if (encoding == FASTQ_ENCODING::SANGER) 
			return "sanger";
		else if (encoding == FASTQ_ENCODING::ILLUMINA) 
			return "illumina";
		else 
			return "unsupported";
	}
}
