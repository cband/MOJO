
#include "PSLParser.h"

namespace MOJO
{
	PSLLine::PSLLine()
	{
		matches = 0, misMatches = 0, qSize = 0, qStart = 0, qEnd = 0;
		blockCount = 0;
		strand = "", qName = "", tName = "";
		blockSizes.clear();
		qStarts.clear();
		tStarts.clear();
	}

	void PSLLine::ParseIn(string line)
	{
		auto sp = Utils::SplitToVector(line, "\t");
		matches = stoi(sp[0]), misMatches = stoi(sp[1]), strand = sp[8];
		qName = sp[9], qSize = stoi(sp[10]), qStart = stoi(sp[11]);
		qEnd = stoi(sp[12]), tName = sp[13], blockCount = stoi(sp[17]);

		vector<string> str_blockSizes = Utils::SplitToVector(sp[18], ",");
		vector<string> str_qStarts = Utils::SplitToVector(sp[19], ",");
		vector<string> str_tStarts = Utils::SplitToVector(sp[20], ",");

		blockSizes.clear();
		qStarts.clear();
		tStarts.clear();

		for (auto str : str_blockSizes) 
			if (str != "")
				blockSizes.push_back(stoi(str));

		for (auto str : str_qStarts) 
			if (str != "")
				qStarts.push_back(stoi(str));

		for (auto str : str_tStarts) 
			if (str != "")
				tStarts.push_back(stoi(str));
	}

	void JunctionPSLLine::ParseIn(string line)
	{
		PSLLine::ParseIn(line);
		end = qName.substr(qName.length() - 1, 1);
	}

	void AnchorReadPSLLine::ParseIn(string line)
	{
		PSLLine::ParseIn(line);
		int pos = qName.find_last_of('_') + 1;
		break_at = stoi(qName.substr(pos, qName.size() - pos));
	}

	bool PSLParser::GetNextLine(PSLLine &psl) 
	{
		if (filestream == 0) 
			filestream = new ifstream(filename.c_str());

		string str;
		if (!getline((*filestream), str))
			return false;
		psl.ParseIn(str);
		return true;
	}

	bool PSLParser::EndOfFile()
	{
		if (filestream == 0) 
			filestream = new ifstream(filename.c_str());

		if (!filestream->is_open() || filestream->eof())
			return true;

		return false;
	}
};