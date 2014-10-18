
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>
#include "Utils.h"

using namespace std;
using namespace boost;
using namespace MOJO;

int main(int argc, char *argv[])
{
	if (argc != 3){
		std::cout << "Correct Usage: ";
		std::cout << "MegablastOutputCompiler <megablast_output_dir> <output_filename>  "
			<< endl;
		return 0;
	}

	string queryOutputDir(argv[1]), outputFileName(argv[2]);

	ofstream outFile(outputFileName.c_str(), ios::out);

	boost::filesystem::directory_iterator iterator(queryOutputDir);
	unordered_map<string, bool> added;
	for (; iterator != boost::filesystem::directory_iterator(); ++iterator) {
		auto fn = iterator->path().filename();

		
		if (!(fn.extension().string() == ".output" &&
			strncmp(fn.string().c_str(), "query.", 6) == 0))
			continue;
		cout << "Processing: " << queryOutputDir << "/" << fn.filename().string() << endl;
		ifstream queryOut((queryOutputDir + "/" + fn.filename().string()).c_str());
		for (string str; getline(queryOut, str);) {
			auto sp = Utils::SplitToVector(str, "\t");
			double identity = stof(sp[2]);

			if (identity < 96 || sp[0] == sp[1])
				continue;
			string key;
			if (sp[0].compare(sp[1]) > 0) {
				key = "g" + sp[0] + "\tg" + sp[1] + "\t" + sp[6] + "\t" + sp[7]
					+ "\t" + sp[8] + "\t" + sp[9];
			}
			else {
				key = "g" + sp[1] + "\tg" + sp[0] + "\t" + sp[8] + "\t" + sp[9]
					+ "\t" + sp[6] + "\t" + sp[7];
			}
			if (added.find(key) != added.end())
				continue;
			added[key] = true;
			outFile << key << endl;
		}
	}

	return 1;
}
