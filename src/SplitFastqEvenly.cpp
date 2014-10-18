
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <vector>
#include <iostream>
#include <fstream>

#include <zlib.h>

using namespace std;

int TrimRead(char* qual, bool toSanger = false);
int getQual(char c);
const int MAX_CHARS_LINE = 1000;
static int quality_offset = 0, quality_threshold = 0;

int main(int argc, char *argv[])
{
	cin.sync_with_stdio(false);
	if (argc < 6) {
		cout << endl << "SplitFastqEvenly - splits fastq (stdin) into specified number of files" << endl << endl;
		cout << "  Usage: cat file.fastq | SplitFastqEvenly <num_splits> <end> <output_file> <trim quality> <to_sanger> <readcount_file>" << endl ;
		cout << "      <num_splits>     - number of splits to create" << endl;
		cout << "      <end>            - 1/2 to designate which end of the read pair the input file represents" << endl;
		cout << "      <output_file>    - /outputdir/outputfile_%d.fastq or /outputdir/outputfile_%d.fastq.gz" << endl;
		cout << "                         if suffix is .gz, output will be compressed" << endl;
		cout << "                         Note: It is optimal to split to fastq and then compress using pigz" << endl;
		cout << "      <trim_quality>   - trim 3' ends of reads to <trim_quality> (bwa trimming algorithm)" << endl;
		cout << "      <to_sanger>      - 1/0  convert to sanger" << endl;
		cout << "      <readcount_file> - total number of reads is written to <readcount_file> [optional]" << endl << endl;
		return 0;
	}

	int num_splits = atoi(argv[1]);
	int end = atoi(argv[2]);

	//vector<gzFile> gzFileStreams;
	//for (int i = 0; i < num_splits; i++) {
	//	char filename[200];
	//	sprintf(filename, argv[3], i);
	//	gzFileStreams.push_back(gzopen(filename, "wb"));
	//}
	vector <ostream *> outputStreams;
	for (int i = 0; i < num_splits; i++) {
		char filename[200];
		sprintf(filename, argv[3], i);
		outputStreams.push_back(new ofstream(filename, ios::out));
	}

	quality_threshold = atoi(argv[4]);
	bool toSanger = atoi(argv[5]) == 1 ? true : false;
	quality_offset = toSanger ? 64 : 33;
	string readcount_filename = ((argc > 6) ? std::string(argv[6]) : "");

	char line[MAX_CHARS_LINE], sequence[MAX_CHARS_LINE], quality[MAX_CHARS_LINE];
	int file_to = 0, read_count = 0, max_read_length = 0, min_read_length = 0;

	while (fgets(line, MAX_CHARS_LINE, stdin)) {
		char *rd = new char[100];
		int i = 1;
		bool propersuffix = false;
		for (; i < strlen(line) - 1; i++) {
			if (line[i] == ' ')
				break;
			if (line[i] == '/')
				propersuffix = true;
			rd[i - 1] = line[i];
		}
		rd[i - 1] = '\0';

		fgets(sequence, MAX_CHARS_LINE, stdin);
		int read_length = strlen(sequence);
		sequence[read_length - 1] = '\0';
		fgets(line, MAX_CHARS_LINE, stdin);
		fgets(quality, MAX_CHARS_LINE, stdin);
		quality[read_length - 1] = '\0';

		int trimTo = TrimRead(quality, toSanger);
		sequence[trimTo] = '\0';
		quality[trimTo] = '\0';
		char buf[1000];
		if (propersuffix)
			sprintf(buf, "@%s\n%s\n+%s\n%s\n", rd, sequence, rd, quality);
		else
			sprintf(buf, "@%s/%d\n%s\n+%s/%d\n%s\n", rd, end, sequence, rd, end, quality);

		(*outputStreams[file_to]).write(buf, strlen(buf));
		//gzprintf(gzFileStreams[file_to], buf);

		if (++file_to == num_splits)
			file_to = 0;
		delete [] rd;

		read_count++;
		if (read_length > max_read_length) max_read_length = read_length;
		if (read_length < min_read_length || min_read_length == 0) min_read_length = read_length;
	}

	for (int i = 0; i < num_splits; i++)
		(*outputStreams[i]).~ostream();
		//gzclose(gzFileStreams[i]);
	
	if (readcount_filename != "") {
		try {
			ofstream o(readcount_filename, ios::out);
			o << "Total read count: " << read_count << ", Min read length: " << min_read_length
				<< ", Max read length: " << max_read_length << endl;
			o.close();
		}
		catch (std::exception &e) {
			std::cout << "Warning: could not write readcount (" << read_count 
				<< ") to '" << readcount_filename << "' with error: " << e.what() << endl;
		}

	}

	return 0;
}

int TrimRead(char *quality, bool toSanger) {
	auto read_length = strlen(quality);

	if (quality_threshold <= 0 || getQual(quality[read_length - 1]) >= quality_threshold) {
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


