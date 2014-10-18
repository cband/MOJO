#ifndef FASTQPARSER_H
#define FASTQPARSER_H

#pragma once

#include <cstring>
#include <fstream>
#include <sstream>

#include <boost/unordered_map.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Utils.h"
#include "Read.h"

using namespace std; 
using namespace boost;

namespace MOJO 
{
	enum FASTQ_ENCODING { 
		UNSUPPORTED, 
		SANGER, 
		ILLUMINA 
	};

	class FastqFile 
	{
		private:
			bool fileOpened;
		public:
			string fileName;
			boost::iostreams::filtering_istream *fileStream;
			ifstream *fileStreamByPos;

			FastqFile() : fileStream(0), fileStreamByPos(0){};

			FastqFile(string file) : 
				fileName(file), fileStream(0), fileStreamByPos(0) {};

			bool GetNextRead(Read *read);

			bool GetReadStartingAtPos(long long pos, Read *read);

			bool EndOfFile();
	};

	class FastqParser 
	{
		private:
			FastqFile *end1File, *end2File;
			string firstEndFile, secondEndFile, indexFile;
			int numReads;
			unordered_map<int, std::pair<long long, long long> > index;
		
		public:
			FastqParser(string firstFile, string secondFile="" );

			FastqParser(string firstFile, string secondFile, string indexFile);

			PairedRead GetNextPairedRead();

			bool GetNextPairedRead(PairedRead *pr);

			bool EndOfFile();

			void CreateTrimmedFiles(string firstEndOutFile, string secondEndOutFile, 
				int maxLength, int startAt = 0);

			FASTQ_ENCODING GetFastqEncoding();

			PairedRead GetReadById(int id);

			static string GetFastqEncodingString( FASTQ_ENCODING encoding );
	};
};

#endif

