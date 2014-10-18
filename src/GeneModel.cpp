
#include "GeneModel.h"

namespace MOJO 
{
	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >, 
		(keywords::channel = ""));

	GeneModel GeneModel::gm;
	map<string, map<int, Gene*> > GeneModel::chromGenesRefSeq;
	
	//Loads the transcriptome model into memory.
	void GeneModel::LoadGeneModel(){
		Config *c = Config::GetConfig();

		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "Loading transcriptome annotation...";
		//Load Genes;
		try {
			ifstream geneFile(c->masterGeneFile.c_str());
			for (std::string str; getline(geneFile, str);) {
				auto sp = Utils::SplitToVector(str, "\r\t");
				Gene *g = new Gene(sp[0], sp[1], sp[2], sp[3], sp[8]);
				g->txStart_genomic = stoi(sp[4]);
				g->txEnd_genomic = stoi(sp[5]);
				g->txStart_transcriptomic = stoi(sp[6]);
				g->txEnd_transcriptomic = stoi(sp[7]);
				GenesMap[g->geneId] = g;
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error loading Genes from "
				<< c->masterGeneFile << endl << "Error: " << e.what();
			exit(1);
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\t" << GenesMap.size() << " genes";
		//Load Isoforms;
		try  {
			ifstream isoformFile(c->masterIsoformFile.c_str());
			for (std::string str; getline(isoformFile, str);) {
				auto sp = Utils::SplitToVector(str, "\r\t");
				Gene *gene = GetGeneFromMap(sp[1]);
				if (gene == 0) {
					BOOST_LOG_CHANNEL(logger::get(), "Main") 
						<< "Error loading isoforms. Cannot find gene id: "
						<< sp[1] << " for line: " << endl << str 
						<< " in file: " << c->masterIsoformFile;
					exit(1);
				}
				Isoform* iso = new Isoform(gene, sp[0], sp[4]);
				iso->txStart = stoi(sp[2]);
				iso->txEnd = stoi(sp[3]);
				iso->cdsStartDistFromTxStart = stoi(sp[5]) - 1;
				iso->cdsEndDistFromTxEnd = stoi(sp[6]);
				IsoformsMap[iso->isoformIdStr] = iso;
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error loading Isoforms from " << c->masterIsoformFile 
				<< endl << "Error: " << e.what();
			exit(1);
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\t" << IsoformsMap.size() << " isoforms";
		//Load Exons;
		try {
			ifstream exonFile(c->masterExonFile.c_str());
			for (std::string str; getline(exonFile, str);) {
				auto sp = Utils::SplitToVector(str, "\r\t");
				Gene *gene = GetGeneFromMap(sp[1]);
				if (gene == 0) {
					BOOST_LOG_CHANNEL(logger::get(), "Main") 
						<< "Error loading exons. Cannot find gene id: " 
						<< sp[1] << " for line: " << endl << str 
						<< " in file: " << c->masterExonFile;
					exit(1);
				}
				Exon *ex = new Exon(gene, sp[0]);
				ex->exStart = stoi(sp[2]);
				ex->exEnd = stoi(sp[3]);
				ex->exStart_genomic = stoi(sp[4]);
				ex->exEnd_genomic = stoi(sp[5]);
				ex->sequence = sp[6];
				ExonsMap[ex->exonId] = ex;
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error loading Exons from " << c->masterExonFile 
				<< endl << "Error: " << e.what();
			exit(1);
		}
		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\t" << ExonsMap.size() << " exons";

		try {
			//Link nested references in Gene, Isoform and Exon hashes
			for (auto iter = GenesMap.begin(); iter != GenesMap.end(); iter++) {
				Gene *gene = (*iter).second;
				unordered_map<int, Exon*> geneExons;
				for (auto isoformIdStr : gene->isoformsStrVect) {
					Isoform *iso = IsoformsMap[isoformIdStr];
					gene->isoforms.push_back(iso);
					for (auto exonIdStr : iso->exonsStrVect) {
						int exonId = stoi(exonIdStr.substr(1, exonIdStr.length() - 1));
						Exon *ex = ExonsMap[exonId];
						iso->exons.push_back(ex);
						geneExons[ex->exonId] = ex;
						ex->exonIsoforms.push_back(iso);
					}
					std::sort(iso->exons.begin(), iso->exons.end(), Exon::compareExons);
				}
				for (auto iter = geneExons.begin(); iter != geneExons.end(); ++iter)
					gene->allExons.push_back((*iter).second);
				std::sort(gene->allExons.begin(), gene->allExons.end(), 
					Exon::compareExons);
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error building GeneModels "
				<< endl << "Error: " << e.what();
			exit(1);
		}
		//BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tLoading gene proximity file...";
		//try {
		//	//Load proximity file
		//	ifstream proximityFile(c->masterGeneProximityFile.c_str());
		//	for (std::string str; getline(proximityFile, str);) {
		//		auto sp = Utils::SplitToVector(str, "\t\r");
		//		Gene *gene = GetGeneFromMap(sp[0]);
		//		if (sp[1] != "") {
		//			gene->prevGeneDist = stoi(sp[2]);
		//			gene->prevGene = GetGeneFromMap(sp[1]);
		//		}
		//		if (sp[3] != "") {
		//			gene->nextGeneDist = stoi(sp[4]);
		//			gene->nextGene = GetGeneFromMap(sp[3]);
		//		}
		//	}
		//}
		//catch (std::exception &e) {
		//	BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error loading proximity file: "
		//		<< c->masterGeneProximityFile << endl << "Error: " << e.what();
		//	exit(1);
		//}
		
		BOOST_LOG_CHANNEL(logger::get(), "Main") << "\tLoading megablast index...";
		try {
			Gene *gA = 0, *gB = 0;
			string gA_str = "", gB_str = "";
			ifstream megablastFile(c->megablastOutputFile.c_str());
			for (std::string str; getline(megablastFile, str);) {
				auto sp = Utils::SplitToVector(str, "\t\r");
				if (gA_str != sp[0]) {
					gA_str = sp[0];
					gA = GetGeneFromMap(gA_str);
				}
				if (gB_str != sp[1]) {
					gB_str = sp[1];
					gB = GetGeneFromMap(gB_str);
				}
				CoordinatePair cpA(stoi(sp[2]), stoi(sp[3]),
					stoi(sp[4]), stoi(sp[5]));
				gA->blastHits[gB->geneId].push_back(cpA);
				CoordinatePair cpB(stoi(sp[4]), stoi(sp[5]),
					stoi(sp[2]), stoi(sp[3]));
				gB->blastHits[gA->geneId].push_back(cpB);
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error loading megablast index: "
				<< c->megablastOutputFile << endl << "Error: " << e.what();
			exit(1);
		}

		BOOST_LOG_CHANNEL(logger::get(), "Main") 
			<< "\tLoading RepeatMasker track...";
		try {
			ifstream repeatFile(c->repeatMaskerFile.c_str());
			for (std::string str; getline(repeatFile, str);) {
				auto sp = Utils::SplitToVector(str, "\t\r");
				Gene *g = GetGeneFromMap(sp[0]);
				// this happens because of an inconsistency in MOJORefBuilder -- to fix
				if (g == 0) 
					continue;
				pair<int, int> p(stoi(sp[2]) - 18, stoi(sp[3]) + 18);
				g->repeatRegions.push_back(p);
			}
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error loading RepeatMasker file: " << c->repeatMaskerFile 
				<< endl << "Error: " << e.what();
			exit(1);
		}
		ValidateGeneModelsIntegrity();
		IsGeneModelLoaded = true;
	} 

	//This function is mainly for debugging purposes
	bool GeneModel::ValidateGeneModelsIntegrity() {
		Config *c = Config::GetConfig();
		for (auto geneIter = GenesMap.begin(); 
			geneIter != GenesMap.end(); geneIter++) 
		{
			// check if gene.allExons are sorted by start 
			// check if exon.exStart is always less tha exon.exEnd
			int prevExonStart = 0;
			for (auto ex : (*geneIter).second->allExons){
				if (prevExonStart > ex->exStart) {
					BOOST_LOG_CHANNEL(logger::get(), "Main") 
						<< "Error in GeneModels: allExons is not properly "
						<< "sorted. See exonid: " << ex->exonId;
					exit(1);
				}
				if (ex->exStart > ex->exEnd) {
					BOOST_LOG_CHANNEL(logger::get(), "Main") 
						<< "Error in GeneModels: exStart is greater than "
						<< "exEnd. See exonid: " << ex->exonId;
					exit(1);
				}
				prevExonStart = ex->exStart;
			}
		}

		//validate proxmity file
		//int missingGenesInProximityFile = 0;
		//stringstream ss;
		//ifstream proximityFile(c->masterGeneProximityFile.c_str());
		//for (std::string str; getline(proximityFile, str);) {
		//	auto sp = Utils::SplitToVector(str, "\t\r");
		//	if (GetGeneFromMap(sp[0]) == 0 || GetGeneFromMap(sp[1]) == 0 ||
		//		GetGeneFromMap(sp[3]) == 0) {
		//			missingGenesInProximityFile++;
		//			ss << str << endl;
		//		}
		//}
		//if (missingGenesInProximityFile > 3) {
		//	BOOST_LOG_CHANNEL(logger::get(), "Main") << "Error in GeneModels: "
		//		<< "missing geneIds in proximity file: " 
		//		<< ss.str().substr(0, min(1000, (int)ss.str().length())) << endl;
		//	exit(1);
		//}
		return false;
	}

	GeneModel* GeneModel::GetGeneModel()
	{
		if (!GeneModel::gm.IsGeneModelLoaded) {
			cerr << "Attempting to access uninitialized gene model. Exiting." << endl;
			exit(1);
		}
		return &GeneModel::gm;
	}

	Gene* GeneModel::GetGeneFromMap(string geneIdStr)
	{
		try {
			//Allow for geneId to be empty string, a space, a 0 or -1.
			if (geneIdStr.length() <= 2)
				return 0;
			int id = stoi(geneIdStr.substr(1, geneIdStr.length() - 1));
			if (GenesMap.find(id) != GenesMap.end())
				return GenesMap[id];
		}
		catch (std::exception &e) {
			//ignore;
		}
		return 0;
	}

	string GeneModel::GetGeneNameForCoordinates(string chr, int pos)
	{
		if (chromGenesRefSeq.size() == 0) {
			GeneModel *gm = GeneModel::GetGeneModel();
			vector<Gene *> genes;
			for (auto g : gm->GenesMap)
				genes.push_back(g.second);
			std::sort(genes.begin(), genes.end(), Gene::compareGenes);
			for (auto g : genes)
				chromGenesRefSeq[g->chr][g->txStart_genomic] = g;
		}
		// Change this algo-- use hash or binary search to make it fast.
		for (auto iter : chromGenesRefSeq[chr]) {
			Gene *g = iter.second;
			if ((pos > g->txStart_genomic && pos < g->txEnd_genomic) || 
				(pos > g->txEnd_genomic && pos < g->txStart_genomic))
				return g->name;
		}
		return "-";
	}
};