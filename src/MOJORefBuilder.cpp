
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

#include <boost/unordered_map.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/lexical_cast.hpp>

#include "GeneModelObjs.h"
#include "Utils.h"

using namespace std;
using namespace boost;
using namespace MOJO;

int next_exon_id = 1000000;
const int TX_PAD_LEN = 500;

int main(int argc, char *argv[])
{
	if (argc != 6){
		std::cout << "Correct Usage: ";
		std::cout << "MojoRefBuilder <MOJO.refName.txt> <genome.fa> <rmsk.txt> "
			<< "<outputdir> <reference_prefix>" << endl;
		return 0;
	}
	string fastaFromBed = "/lustre/beagle/cbandlam/libs/bedtools/bin/fastaFromBed";
	string faToTwoBit = "/lustre/beagle/cbandlam/libs/blattools/faToTwoBit";
	string bowtie2build = "/lustre/beagle/cbandlam/MOJO-P/external/bowtie2-2.2.3/bowtie2-build";
	string bwa = "/lustre/beagle/cbandlam/MOJO-P/external/bwa-0.7.10/bwa";
	string blastDir = "/lustre/beagle/cbandlam/libs/blast-2.2.26/bin/";

	bool debug = false;
	string mojoFormattedRefName(argv[1]), genomeFasta(argv[2]);
	string rmskFName(argv[3]), outputdir(argv[4]), refprefix(argv[5]);
	if (outputdir.back() != '/')
		outputdir = outputdir + "/";
	outputdir += "reference." + refprefix + "/";

	string outPrefix = outputdir + "/" + refprefix + ".";

	string blat_filter_dir = outputdir + "/blat_filter_refs/";
	string all_isoforms_dir = outputdir + "/all_isoforms/";
	string genome_dir = outputdir + "/genome/";
	string transcriptome_dir = outputdir + "/transcriptome/";
	string gene_model_dir = outputdir + "/gene_model/";
	
	boost::filesystem::create_directories(blat_filter_dir);
	boost::filesystem::create_directories(all_isoforms_dir);
	boost::filesystem::create_directories(genome_dir);
	boost::filesystem::create_directories(transcriptome_dir);
	boost::filesystem::create_directories(gene_model_dir);

	ofstream txomeBedOutStream((outPrefix + "transcriptome.bed").c_str(), ios::out);
	ofstream geneOutStream((gene_model_dir + "Gene.txt").c_str(), ios::out);
	ofstream isoformOutStream((gene_model_dir + "Isoform.txt").c_str(), ios::out);
	ofstream exonOutStream((gene_model_dir + "Exon.txt").c_str(), ios::out);

	ifstream refFileStream(mojoFormattedRefName.c_str());
	unordered_map<int, vector<string> > geneIsoformStrMap;
	for (string str; getline(refFileStream, str);) {
		auto sp = Utils::SplitToVector(str, "\t");
		if (sp[3].find('_') != string::npos)
			continue;
		geneIsoformStrMap[stoi(sp[0])].push_back(str);
	}

	// the following StreamMap containers are required to exclude any genes
	// for which sequence (from .bed) cannot be extracted. Eg: genes at the
	// ends of the mitochondrial chromosome
	map<int, string> geneOutStreamMap;
	map<int, vector<string> > isoformOutStreamMap;
	map<int, string> exonOutStreamMap;
	map<string, std::pair<int, int> > isoToCDS;

	map<string, int> exonNames;
	map<string, Exon*> exonsMap;
	map<int, vector<Exon*> > geneToExonVectMap;
	vector<Gene *> genesVect;
	map<int, Gene *> genesMap;
	unordered_map<string, bool> isoformAlreadyProcessed;
	for (auto geneIsoVect : geneIsoformStrMap) {
		Gene *gene = new Gene();
		gene->geneId = geneIsoVect.first;
		genesMap[gene->geneId] = gene;
		genesVect.push_back(gene);
		string isoformStrCsv;
		map<string, string> exonIDmap;
		bool firstExon = true;
		for (auto isoformLine : geneIsoVect.second) {
			auto sp = Utils::SplitToVector(isoformLine, "\t");
			string isoformName = sp[2];

			if (isoformAlreadyProcessed.find(isoformName) != 
				isoformAlreadyProcessed.end()) 
			{
				std::cout << "Duplicate isoform detected: " 
					<< isoformName << ". Exiting." << endl;
				exit(1);
			}
			Isoform *iso = new Isoform();
			iso->isoformIdStr = isoformName;
			isoformStrCsv += isoformName + ",";
			iso->gene = gene;
			gene->chr = sp[3], gene->strand = sp[4];
			gene->name = sp[1];
			boost::replace_all(gene->name, " ", "_");

			vector<string> exonStarts, exonEnds;
			if (sp[10].c_str()[sp[10].length() - 1] == ',')
				sp[10].pop_back();
			if (sp[11].c_str()[sp[11].length() - 1] == ',')
				sp[11].pop_back();
			exonStarts = Utils::SplitToVector(sp[10], ",");
			exonEnds = Utils::SplitToVector(sp[11], ",");

			std::pair<int, int> p(stoi(sp[7]), stoi(sp[8]));
			isoToCDS[isoformName] = p;

			//if (sp[7] == sp[8])
			//	upstreamUTR = -1, downstreamUTR = -1;
			iso->cdsStartDistFromTxStart = -1;
			iso->cdsEndDistFromTxEnd = -1;
			if (gene->strand == "+") {
				iso->txStart_genomic = stoi(exonStarts[0]);
				iso->txEnd_genomic = stoi(exonEnds[exonEnds.size() - 1]);
			}
			else {
				iso->txStart_genomic = stoi(exonEnds[exonEnds.size() - 1]);
				iso->txEnd_genomic = stoi(exonStarts[0]);
			}

			if (firstExon) {
				gene->txStart_genomic = stoi(exonStarts[0]);
				gene->txEnd_genomic = stoi(exonEnds[exonEnds.size() - 1]);
			}
			else {
				if (gene->txStart_genomic > stoi(exonStarts[0]))
					gene->txStart_genomic = stoi(exonStarts[0]);
				if (gene->txEnd_genomic < stoi(exonEnds[exonEnds.size() - 1]))
					gene->txEnd_genomic = stoi(exonEnds[exonEnds.size() - 1]);
			}
			firstExon = false;
			for (int ex_idx = 0; ex_idx < exonStarts.size(); ex_idx++) {
				string exon_pos_tag = gene->chr + "|" + gene->strand + "|" +
					exonStarts[ex_idx] + "|" + exonEnds[ex_idx];

				if (exonIDmap.find(exon_pos_tag) == exonIDmap.end())  {
					string idStr = "e" + lexical_cast<string>(next_exon_id++);
					Exon *ex = new Exon(gene, idStr);
					exonsMap[idStr] = ex;
					exonsMap[idStr]->exStart_genomic = stoi(exonStarts[ex_idx]);
					exonsMap[idStr]->exEnd_genomic = stoi(exonEnds[ex_idx]);
					exonIDmap[exon_pos_tag] = idStr;
					geneToExonVectMap[gene->geneId].push_back(ex);
				}
				iso->exons.push_back(exonsMap[exonIDmap[exon_pos_tag]]);
			}
			gene->isoforms.push_back(iso);
		}

		//gene->ref_name = gene->ref_name.substr(1, gene->ref_name.size() - 2);
		gene->txStart_genomic = max(0, gene->txStart_genomic - TX_PAD_LEN);
		gene->txEnd_genomic = gene->txEnd_genomic + TX_PAD_LEN;

		txomeBedOutStream << gene->chr << "\t" << gene->txStart_genomic << "\t"
			<< gene->txEnd_genomic << "\tg" << gene->geneId << "\t0\t"
			<< gene->strand << "\t" << endl;

		gene->txStart_transcriptomic = 0;
		gene->txEnd_transcriptomic = 
			abs(gene->txEnd_genomic - gene->txStart_genomic);

		char buf[5000];
		std::sprintf(buf, "g%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n", gene->geneId,
			gene->name.c_str(), gene->chr.c_str(), gene->strand.c_str(),
			gene->txStart_genomic, gene->txEnd_genomic, gene->txStart_transcriptomic,
			gene->txEnd_transcriptomic, isoformStrCsv.c_str());
		string s(buf);
		geneOutStreamMap[gene->geneId] = s;
	}
	txomeBedOutStream.close();

	//Properly set the downstream and upstream UTRs
	for (auto gene : genesVect) {
		unordered_map<string, bool> exonWritten;
		for (auto iso : gene->isoforms) {
			iso->txStart = abs(gene->txStart_genomic - iso->txStart_genomic);
			iso->txEnd = abs(gene->txStart_genomic - iso->txEnd_genomic);

			string exons_csv = "";
			int startpos = (gene->strand == "+" ?
				gene->txStart_genomic : gene->txEnd_genomic);

			for (int j = 0; j < iso->exons.size(); j++) {
				Exon *exon = iso->exons[j];
				if (gene->strand == "+") {
					exon->exStart = abs(exon->exStart_genomic - startpos);
					exon->exEnd = abs(exon->exEnd_genomic - startpos);
				}
				else {
					exon->exEnd = abs(exon->exStart_genomic - startpos);
					exon->exStart = abs(exon->exEnd_genomic - startpos);
				}
				string exonIdStr = lexical_cast<string>(exon->exonId);
				exons_csv += "e" + exonIdStr + ",";
			}

			//
			// Find Upstream/downstream UTR lengths (cdsStart/End)
			//
			std::pair<int, int> p = isoToCDS[iso->isoformIdStr];
			if (p.first != p.second) {
				std::sort(iso->exons.begin(), iso->exons.end(), 
					Exon::compareExonsGenomic);
				unordered_map<int, bool> upstream, cds, downstream;
				int step = 0; // 0 = upstream, 1= cds, 2 = downstream
				for (auto e : iso->exons) {
					for (int i = e->exStart_genomic; i < e->exEnd_genomic; i++) {
						if (step == 0 && i == p.first)
							step = 1;
						if (step == 1 && i == p.second)
							step = 2;

						if (step == 0) upstream[i] = true;
						if (step == 1) cds[i] = true;
						if (step == 2) downstream[i] = true;
					}
				}
				if (gene->strand == "+") {
					iso->cdsStartDistFromTxStart = upstream.size() + 1;
					iso->cdsEndDistFromTxEnd = downstream.size();
				}
				else {
					iso->cdsStartDistFromTxStart = downstream.size() + 1;
					iso->cdsEndDistFromTxEnd = upstream.size();
				}
			}
			stringstream ss;
			ss << iso->isoformIdStr << "\tg" << iso->gene->geneId
				<< "\t" << iso->txStart << "\t" << iso->txEnd << "\t" << exons_csv;
			ss << "\t" << iso->cdsStartDistFromTxStart << "\t"
				<< iso->cdsEndDistFromTxEnd << endl;
			isoformOutStreamMap[iso->gene->geneId].push_back(ss.str());
		}
	}

	//populate gene.allExons;
	for (auto g : genesVect) {
		unordered_map<int, Exon *> tmpExonMap;
		for (auto i : g->isoforms) 
			for (auto e : i->exons) 
				tmpExonMap[e->exonId] = e;
		for (auto t : tmpExonMap)
			g->allExons.push_back(t.second);
	}


	// use fastaFromBed (bedtools) to create transcriptome.fa
	string txBed = outPrefix + "transcriptome.bed";
	string txFa = outPrefix + "transcriptome.fa";
	char cmd[5000];
	std::sprintf(cmd, "%s -fi %s -bed %s -fo %s -name -s ",
		fastaFromBed.c_str(), genomeFasta.c_str(), txBed.c_str(), txFa.c_str());
	Utils::ExecuteCommand(cmd, "", true, true);

	//Create exons.fa
	ifstream txFaStream(txFa.c_str());
	for (string str; getline(txFaStream, str);){
		int geneId = stoi(str.substr(2, str.length() - 2));
		if (genesMap.find(geneId) == genesMap.end()) {
			std::cout << "Missing gene id: " << geneId << endl;
			exit(1);
		}
		geneOutStream<<geneOutStreamMap[geneId];
		for (auto line : isoformOutStreamMap[geneId]) 
			isoformOutStream << line;

		Gene *gene = genesMap[geneId];
		string fa;
		getline(txFaStream, fa);

		auto exons = geneToExonVectMap[geneId];
		for (auto exon : exons){
			exon->sequence = fa.substr(exon->exStart, exon->exEnd - exon->exStart);
			exonOutStream << "e" << exon->exonId << "\tg" << exon->gene->geneId 
				<< "\t" << exon->exStart << "\t" << exon->exEnd << "\t"
				<< exon->exStart_genomic << "\t" << exon->exEnd_genomic
				<< "\t" << exon->sequence << endl;
		}
	}
	geneOutStream.close();
	exonOutStream.close();
	isoformOutStream.close();

	// Create allisoforms Fa
	ofstream allIsoStream((all_isoforms_dir + "all.isoforms.fa").c_str(), ios::out);
	for (auto gene : genesVect) {
		if (debug)
			continue;
		allIsoStream << ">g" << gene->geneId << endl;
		for (auto iso : gene->isoforms)
			allIsoStream << iso->GetTranscript();
		allIsoStream << endl;
	}
	allIsoStream.close();

	//Construct allReferenceJunction blat files
	string allRefJunctFa = blat_filter_dir + "allReferenceJunctsDB.fa";
	ofstream allRefJunctOutStream(allRefJunctFa.c_str(), ios::out);
	unordered_map<string, bool> seqAdded;
	for (auto gene : genesVect) {
		if (debug)
			continue;
		std::sort(gene->allExons.begin(), gene->allExons.end(), 
			Exon::compareExons);
		allRefJunctOutStream << ">g" << gene->geneId << "_Collapsed" << endl;
		allRefJunctOutStream << gene->GetTranscribedSequence() << endl;
		//for practical reasons, skip genes with > 200 exons;  Only skips ~5 genes;
		if (gene->allExons.size() > 200) 
			continue;
		for (auto iso : gene->isoforms) {
			std::sort(iso->exons.begin(), iso->exons.end(), Exon::compareExons);
			for (int i = 0; i < iso->exons.size(); i++) {
				for (int j = i + 1; j < min( i + 50, (int)iso->exons.size()); j++) {
					string a = iso->exons[i]->sequence;
					string b = iso->exons[j]->sequence;
					if (a.length() > 100)
						a = a.substr(a.length() - 100, 100);
					if (b.length() > 100)
						b = b.substr(0, 100);
					if (seqAdded.find(a + b) != seqAdded.end())
						continue;
					seqAdded[a + b] = true;
					allRefJunctOutStream << ">" << iso->isoformIdStr << "_" << i 
						<< "_" << j << endl << a << b << endl;
				}
			}
		}
	}
	allRefJunctOutStream.close();

	//convert allReferenceJunctions to 2bit
	sprintf(cmd, "%s %sallReferenceJunctsDB.fa %sallReferenceJunctsDB.2bit",
		faToTwoBit.c_str(), blat_filter_dir.c_str(), blat_filter_dir.c_str());
	Utils::ExecuteCommand(cmd, "", true, true);

	sprintf(cmd, "%s %s %s/genome.2bit",
		faToTwoBit.c_str(), genomeFasta.c_str(), blat_filter_dir.c_str());
	Utils::ExecuteCommand(cmd, "", true, true);

	//copy other .2bit files to reference directory
	std::size_t pos = genomeFasta.find_last_of("/");
	string chr2bit_dir = genomeFasta.substr(0, pos);
	sprintf(cmd, "cp %s/chr*2bit %s/", chr2bit_dir.c_str(), blat_filter_dir.c_str());
	Utils::ExecuteCommand(cmd, "", true, true);

	// find repeatMasker regions in genes
	unordered_map<string, unordered_map<int, vector< Gene * > > > geneCoordsIdx;
	int bin = 1000000;
	for (auto g : genesVect) {
		int bucket_start = g->txStart_genomic / bin;
		int bucket_end = g->txEnd_genomic / bin;
		for (int bucket = bucket_start; bucket <= bucket_end; bucket++)
			geneCoordsIdx[g->chr][bucket].push_back(g);
	}
	
	ifstream rmskStream(rmskFName.c_str());
	ofstream rmskOutStream((gene_model_dir + "rmsk.regions.txt").c_str(), ios::out);
	for (string str; getline(rmskStream, str);) {
		if (debug)
			continue;
		auto sp = Utils::SplitToVector(str, "\t");
		string chr = sp[5];
		int start = stoi(sp[6]), end = stoi(sp[7]);
		int mismatches = stoi(sp[2]), del = stoi(sp[3]), ins = stoi(sp[4]);
		int len = abs(start - end);

		if (geneCoordsIdx[chr].find(start / bin) == geneCoordsIdx[chr].end())
			continue;
		vector<Gene *> *genes = &geneCoordsIdx[chr][start / bin];
		for ( auto g : (*genes)) {
			if (g->txStart_genomic <= start && g->txEnd_genomic >= end)
			{
				rmskOutStream << "g" << g->geneId << "\t" << g->name << "\t"
					<< (start - g->txStart_genomic) << "\t"
					<< (start - g->txStart_genomic + len) << "\t"
					<< chr << "\t" << start << "\t" << end << endl;
			}
		}
	}
	rmskOutStream.close();

	//Generate megablast index
	// /lustre/beagle/cbandlam/libs/blast-2.2.26/bin/formatdb -i transcriptome.fa -p F
	// /lustre/beagle/cbandlam/libs/blast-2.2.26/bin/megablast -d transcriptome.fa 
	//							-i query.fa -e 1e-20 -m 8 -a 24 -p 87 -o query.output
	stringstream cmdSS;
	string mbPrefix = outputdir + "/" + refprefix + ".megablast/";
	boost::filesystem::create_directories(mbPrefix);
	//create megablast database of the transcriptome;
	sprintf(cmd, "%s/formatdb -i %s -p F -n %s/transcriptome &", 
		blastDir.c_str(), txFa.c_str(), mbPrefix.c_str());
	cmdSS << cmd << endl << endl;

	//Construct allIsoforms Bowtie2 index
	std::sprintf(cmd, "%s %sall.isoforms.fa %sall.isoforms &",
		bowtie2build.c_str(), all_isoforms_dir.c_str(), all_isoforms_dir.c_str());
	cmdSS << cmd << endl << endl;

	//Construct genome bowtie2 index
	std::sprintf(cmd, "%s %s %sgenome & ",
		bowtie2build.c_str(), genomeFasta.c_str(), genome_dir.c_str());
	cmdSS << cmd << endl << endl;

	//Construct transcriptome BWA index
	std::sprintf(cmd, "%s index %s -p %stranscriptome & ",
		bwa.c_str(), txFa.c_str(), transcriptome_dir.c_str());
	cmdSS << cmd << endl << endl;

	Utils::ExecuteCommand(cmdSS.str().c_str(), "", true, true);

	return 1;
}
