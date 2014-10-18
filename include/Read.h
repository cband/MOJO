#ifndef READ_H
#define READ_H

#pragma once

#include <cstring>
#include <iostream>

#include "Utils.h"

using namespace std;
using namespace boost;

namespace MOJO 
{
	class ReadAlignment 
	{
		public:
			void *gene;
			int mismatches;
			int position;
	};

	class Read 
	{
		public:
			string ReadName, Sequence, Quality;
			vector<ReadAlignment> alignments;

			Read() {};

			static string TrimReadName(string read, string findChr = "@");
			
			string GetTrimmedReadName();
			
			void GetFastqString(char *buf, int rtrim = 0, int ltrim = 0,
				bool convertToSanger = false);
			
			void GetFastqSubstring(char *buf, int maxLength, int startAt);
			
			void Initialize();
			
			Read& operator=(const Read& rd) 
			{
				(*this).ReadName = rd.ReadName;
				(*this).Sequence = rd.Sequence;
				(*this).Quality = rd.Quality;
				return *this;
			}

			bool AddAlignment(ReadAlignment al);
	};

	class PairedRead 
	{
		public:
			Read FirstRead;
			Read SecondRead;

			string GetTrimmedReadName() 
			{
				return FirstRead.GetTrimmedReadName();
			}

			void Initialize();
	};
};

#endif
