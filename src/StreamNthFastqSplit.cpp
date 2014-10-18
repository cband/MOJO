
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

using namespace std;

int TrimRead(char* qual, bool toSanger = false);
int getQual(char c);
const int MAX_CHARS_LINE = 1000;
static int quality_offset = 0, quality_threshold = 0;

int main(int argc, char *argv[])
{
	cin.sync_with_stdio(false);
	if (argc != 9) {
		cout << endl << "SplitFastqEvenly - splits fastq (stdin) into specified number of files" << endl << endl;
		cout << "  Usage: cat file.fastq | SplitFastqEvenly <num_splits> <end> <output_file> <trim quality> <to_sanger> <readcount_file>" << endl;
		cout << "      <num_splits>     		  - number of splits to create" << endl;
		cout << "      <Nth split>      		  - Nth split to generate; print every Nth read; for 6 splits, go from 0-5" << endl;
		cout << "      <end>            		  - 1/2 to designate which end of the read pair the input file represents" << endl;
		cout << "      <trim_quality>   		  - trim 3' ends of reads to <trim_quality> (bwa trimming algorithm)" << endl;
		cout << "      <to_sanger>      		  - 1/0  convert to sanger" << endl;
		cout << "      <numeric name>			  - 1/0 1: convert read names to integers" << endl;
		cout << "      <readcount_file> 		  - total number of reads is written to <readcount_file> " << endl;
		cout << "                        		    Use '-' to skip writing the readcount_file" << endl;
		cout << "      <save_numeric_readname>    - if <numeric_name> is 1, then readnames to ID map is written to this file " << endl;
		cout << "                        		    Use '-' to skip writing the readname_to_id mapping" << endl << endl;
		return 0;
	}

	int num_splits = atoi(argv[1]);
	int n_split = atoi(argv[2]); 
	int end = atoi(argv[3]);
	
	quality_threshold = atoi(argv[4]);
	bool toSanger = atoi(argv[5]) == 1 ? true : false;
	quality_offset = toSanger ? 64 : 33;
	bool numericReadName = atoi(argv[6]) == 1 ? true : false;

	string readcount_filename = std::string(argv[7]);
	string readmap_filename = std::string(argv[8]);

	ofstream *readMapFile = 0;
	if (numericReadName && readmap_filename != "-")
		readMapFile = new ofstream(readmap_filename.c_str(), ios::out);

	char tmpline[MAX_CHARS_LINE], readline[MAX_CHARS_LINE];
	char sequence[MAX_CHARS_LINE], quality[MAX_CHARS_LINE];
	int read_count = 0, max_read_length = 0, min_read_length = 0;
	while (fgets(readline, MAX_CHARS_LINE, stdin)) {
		fgets(sequence, MAX_CHARS_LINE, stdin);
		fgets(tmpline, MAX_CHARS_LINE, stdin);
		fgets(quality, MAX_CHARS_LINE, stdin);
		
		if ((num_splits + read_count++) % num_splits != n_split)
			continue;

		char *rd = new char[100];
		int i = 1;
		for (; i < strlen(readline) - 1; i++) {
			if (readline[i] == ' ')
				break;
			if (readline[i] == '/')
				break;
			rd[i - 1] = readline[i];
		}
		rd[i - 1] = '\0';
		int read_length = strlen(sequence);
		sequence[read_length - 1] = '\0';
		quality[read_length - 1] = '\0';

		int trimTo = TrimRead(quality, toSanger);
		sequence[trimTo] = '\0';
		quality[trimTo] = '\0';
		char buf[1000];
		string readname = std::string(rd);
		if (numericReadName) {
			readname = to_string(100000000 + read_count);
			if (readMapFile != 0)
				(*readMapFile) << readname << "\t" << rd << endl;
		}
		sprintf(buf, "@%s/%d\n%s\n+%s/%d\n%s\n",
			readname.c_str(), end, sequence, readname.c_str(), end, quality);

		cout << buf;
		
		delete[] rd;

		if (read_length > max_read_length)
			max_read_length = read_length;
		if (read_length < min_read_length || min_read_length == 0)
			min_read_length = read_length;
	}

	if (readcount_filename != "-") {
		try {
			ofstream o(readcount_filename, ios::out);
			o << "Total read count: " << read_count << ", Min read length: " 
			<< min_read_length << ", Max read length: " << max_read_length << endl;
			o.close();
		}
		catch (std::exception &e) {
			std::cout << "Warning: could not write readcount (" << read_count
				<< ") to '" << readcount_filename << "' with error: "
				<< e.what() << endl;
		}
	}
	return 0;
}

int TrimRead(char *quality, bool toSanger) {
	auto read_length = strlen(quality);

	if (quality_threshold <= 0 ||
		getQual(quality[read_length - 1]) >= quality_threshold)
	{
		if (toSanger)
			for (size_t i = 0; i < read_length; i++)
				quality[i] = (char)((int)quality[i] - 31);
		return read_length;
	}

	int cum_quals[500];
	for (size_t i = 0; i < read_length; i++) {
		int x = i + 1;
		cum_quals[x] = 0;
		for (size_t j = x; j < read_length; j++)
			cum_quals[x] += quality_threshold - getQual(quality[j]);
		if (toSanger)
			quality[i] = (char)((int)quality[i] - 31);
	}

	int trim_till = 0, maxval = 0;
	for (size_t i = 0; i < read_length; i++) {
		if (maxval < cum_quals[i]) {
			trim_till = i;
			maxval = cum_quals[i];
		}
	}

	if (trim_till < 36)
		return 36;
	else
		return trim_till;
}

int getQual(char c) {
	return ((int)c - quality_offset);
}


