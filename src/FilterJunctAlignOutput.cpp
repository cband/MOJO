
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <ctype.h>

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include "Utils.h"

using namespace std;
using namespace boost;
using namespace MOJO;

const int MAX_CHARS_LINE = 1000;
struct Cigar 
{
	char type;
	int length;

	Cigar(char ty, int len) 
	{
		type = ty, length = len;
	}
};

int main(int argc, char *argv[])
{
	cin.sync_with_stdio(false);
	if (argc != 5) {
		cout << "Correct usage: cat <bowtie2-output> | FilterBowtie2JunctOutput "
			"<outfile> <anchor> <error-rate> <max-mismatches>" << endl;
		exit(0);
	}

	ofstream bowtie2sam(argv[1], ios::out);
	int anchor_length = atoi(argv[2]);
	double align_error_rate = atof(argv[3]);
	int allowed_mismatches_in_anchor = atoi(argv[4]);

	unordered_map<string, int> break_positions_map_str;

	char line[MAX_CHARS_LINE];
	regex re(".*_\\w+_\\w+_\\w_(\\d+)$");

	fgets(line, MAX_CHARS_LINE, stdin);
	while (fgets(line, MAX_CHARS_LINE, stdin) != NULL) {
		char buf[1000];
		strcpy(buf, line);
		if (line[0] == '@') {
			string tag(strtok(line, "\t"));
			string refname(strtok(NULL, "\t"));
			if (tag != "@SQ")
				break;
			smatch match;
			regex_match(refname, match, re);
			int break_junct_pos = atoi(match[1].str().c_str());
			break_positions_map_str[refname.substr(3, refname.length() - 3)] = 
				break_junct_pos;
		}
	}
	fgets(line, MAX_CHARS_LINE, stdin);
	while (fgets(line, MAX_CHARS_LINE, stdin) != NULL) {
		char buf[2000];
		sprintf(buf, "%s", line);
		vector<string> sp = Utils::SplitToVector(std::string(line), "\t");
		
		strcpy(buf, line);
		string readname(sp[0]);
		int flag = stoi(sp[1]);
		string strand = ((flag & 16) ? "-" : "+");

		string refname(sp[2]);
		int position = stoi(sp[3]);
		const char *cigar = sp[5].c_str();
		string sequence(sp[9]);

		int mismatches = 0, edit_distance = 0;
		vector<int> mismatch_positions;

		for (vector<string>::size_type i = 11; i < sp.size(); i++) {
			string t = sp[i];
			const char* tok = t.c_str();
			if (strlen(tok) > 3 && strncmp(tok, "XM:", 3) == 0) 
				mismatches = stoi(t.substr(5, t.length() - 5));

			if (strlen(tok) > 3 && strncmp(tok, "NM:", 3) == 0) 
				edit_distance = stoi(t.substr(5, t.length() - 5));

			if (strlen(tok) > 3 && strncmp(tok, "MD:", 3) == 0) {
				string mm_str = t.substr(5, t.length() - 5);
				const char *mm = mm_str.c_str();
				int _carat = 0, cur_pos = 0;

				for (size_t j = 0; j < strlen(mm); j++) {
					if (isalpha(mm[j])) {
						cur_pos += stoi(mm_str.substr(_carat, j - _carat));
						mismatch_positions.push_back(cur_pos);
						_carat = j + 1;
					}
					else if (mm[j] == '^') {
						cur_pos += stoi(mm_str.substr(_carat, j - _carat));
						for (size_t k = j + 1; k < strlen(mm); k++) {
							if (!isalpha(mm[k]) || mm[k] == '^') {
								j = k - 1;
								break;
							}
						}
						_carat = j + 1;
					}
				}
			}
		}

		int delbefore = 0, delafter = 0, _carat = 0, cur_pos = 0;
		vector<Cigar> cigar_vect;
		vector<int> indels_vect;
		for (size_t i = 0; i < strlen(cigar); i++) {
			if (isalpha(cigar[i]))  {
				Cigar c(cigar[i], stoi(sp[5].substr(_carat, i - _carat)));
				cigar_vect.push_back(c);
				_carat = i + 1;
				if (c.type == 'H' || c.type == 'S') 
					continue;
				if (c.type == 'I' || c.type == 'D')
					indels_vect.push_back(cur_pos + 1);
				cur_pos += c.length;
			}
		}
		if (cigar_vect[0].type == 'H' || cigar_vect[0].type == 'S')
			delbefore = cigar_vect[0].length;

		if (cigar_vect[cigar_vect.size() - 1].type == 'H' || 
			cigar_vect[cigar_vect.size() - 1].type == 'S')
			delafter = cigar_vect[cigar_vect.size() - 1].length;

		//Ignore if more than 1 indel  --7/3/2013 -- ALLOW ZERO INDELS; indels result in out-of-frame
		if (indels_vect.size() > 0) 
			continue;
		
		bool skipBecauseOfIndelSize = false;
		
		for (vector<Cigar>::size_type m = 0; m < cigar_vect.size(); m++)
			if ((cigar_vect[m].type == 'I' || cigar_vect[m].type == 'D') && 
				(cigar_vect[m].length > 1))
				skipBecauseOfIndelSize = true;
		if (skipBecauseOfIndelSize)
			continue;

		if ((mismatch_positions.size() + indels_vect.size()) > 
			ceil(align_error_rate * sequence.length())) 
			continue;

		int alignedLength = sequence.length() - delbefore - delafter;
		int end = position + alignedLength;
		int break_junct_pos = break_positions_map_str[refname];
		int left_anchor = break_junct_pos - position + 1;
		int right_anchor = position + alignedLength - break_junct_pos - 1;

		//Nov 7, 2013
		//bowtie2 --local trims reads from both ends; make sure that reads are only trimmed from the 3' end.
		if ((strand == "+" && delbefore > 2) || (strand == "-" && delafter > 2))
			if (min(left_anchor, right_anchor) < 20)
				continue;

		// Mismatches/Indels in +/- 15bp anchor region:  mismatches - max of 1, indels - 0
		// 15bp is a hard cut off, 
		bool indelFound = false;
		int mismatches_in_anchor = 0;
		for (int i = 0; i < indels_vect.size() && !indelFound; i++)
			if ((indels_vect[i] >= (left_anchor - anchor_length)) && 
				(indels_vect[i] <= (left_anchor + anchor_length)))
				indelFound = true;
		
		if (indelFound)  
			continue;

		for (int i = 0; i < mismatch_positions.size(); i++)
			if ((mismatch_positions[i] >= (left_anchor - anchor_length)) && 
				(mismatch_positions[i] <= (left_anchor + anchor_length)))
				mismatches_in_anchor++;
		if (mismatches_in_anchor > allowed_mismatches_in_anchor) 
			continue;

		if (!(position < break_junct_pos &&
			min(left_anchor, right_anchor) >= anchor_length))
			continue;

		string indels_csv = "", mismatches_csv = "";
		for (int i = 0; i < mismatch_positions.size(); i++)
			mismatches_csv += lexical_cast<string>(mismatch_positions[i]) + ",";

		stringstream ss;
		//ref_name, cigar, strand, flag, position, readsequence, 
		//sequence, editdistance, mismatches, delbefore, delafter
		ss << refname << "\t" << cigar << "\t" << strand << "\t" << flag
			<< "\t" << position << "\t" << alignedLength << "\t" << sequence.length()
			<< "\t" << readname << "\t" << sequence << "\t" << edit_distance 
			<< "\t" << mismatches << "\t" << delbefore << "\t" << delafter
			<< "\t" << indels_csv << "\t" << mismatches_csv << "\t" 
			<< (mismatches_in_anchor ? 1 : 0 ) << endl;
		bowtie2sam << ss.str();
	}
	return 0;
}
