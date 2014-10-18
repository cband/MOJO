#ifndef PSL_PARSER_H
#define PSL_PARSER_H

#pragma once

#include "Utils.h"

namespace MOJO 
{

	class PSLLine
	{
		public:
			int matches, misMatches;
			string strand, qName;
			int qSize, qStart, qEnd;
			string tName;
			int blockCount;
			vector<int> blockSizes, qStarts, tStarts;

			PSLLine();

			virtual void ParseIn(string line);
	};

	class JunctionPSLLine : public PSLLine
	{
		public:
			string end;
			
			virtual void ParseIn (string line) override;
	};

	class AnchorReadPSLLine : public PSLLine
	{
		public:
			int break_at;
			
			virtual void ParseIn (string line) override;
	};

	class PSLParser
	{
		private:
			string filename;
			ifstream *filestream;

		public:
			PSLParser(string file) : filename(file), filestream(0) {};
			
			bool GetNextLine(PSLLine &psl);
			
			bool EndOfFile();
	};
}

#endif