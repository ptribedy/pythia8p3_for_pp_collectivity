#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

// For ROOT 5 / ACLiC: generate dictionaries for the STL types we use
#ifdef __MAKECINT__
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<int>+;
#endif

void log_to_root(TString inFileName = "", TString outFileName = "") {

	// ---------------------------------------------------------------------
	// Resolve input / output names from env if not passed as arguments
	// ---------------------------------------------------------------------
	if (inFileName == "") {
		const char* envIn = gSystem->Getenv("LOGFILE");
		if (!envIn) {
			std::cout << "Error: LOGFILE env var not set and no input argument." << std::endl;
			return;
		}
		inFileName = envIn;
	}

	if (outFileName == "") {
		const char* envOut = gSystem->Getenv("ROOTFILE");
		if (!envOut) {
			std::cout << "Error: ROOTFILE env var not set and no output argument." << std::endl;
			return;
		}
		outFileName = envOut;
	}

	std::cout << "log_to_root: input = " << inFileName
		<< ", output = " << outFileName << std::endl;

	// ---------------------------------------------------------------------
	// Open input log file
	// ---------------------------------------------------------------------
	std::ifstream infile(inFileName.Data());
	if (!infile.is_open()) {
		std::cout << "Error: Cannot open input file: " << inFileName << std::endl;
		return;
	}

	// ---------------------------------------------------------------------
	// Create output ROOT file and TTree
	// ---------------------------------------------------------------------
	TFile* fOut = new TFile(outFileName, "RECREATE");
	if (!fOut || fOut->IsZombie()) {
		std::cout << "Error: Cannot create output ROOT file: " << outFileName << std::endl;
		infile.close();
		return;
	}

	TTree* tree = new TTree("tree", "Pythia Event Tree");

	int  b_eventID = -1;
	int  b_nTracks = 0;
	bool b_isStar  = false;
	bool b_isSPhenix = false;

	std::vector<float> vec_pt, vec_eta, vec_phi, vec_mass;
	std::vector<int>   vec_charge, vec_pid;

	tree->Branch("EventID",   &b_eventID,   "EventID/I");
	tree->Branch("nTracks",   &b_nTracks,   "nTracks/I");
	tree->Branch("isStar",    &b_isStar,    "isStar/O");
	tree->Branch("isSPhenix", &b_isSPhenix, "isSPhenix/O");

	// STL vector branches (use STL dicts via #pragma link above)
	tree->Branch("pT",     "std::vector<float>", &vec_pt,    32000, 0);
	tree->Branch("eta",    "std::vector<float>", &vec_eta,   32000, 0);
	tree->Branch("phi",    "std::vector<float>", &vec_phi,   32000, 0);
	tree->Branch("mass",   "std::vector<float>", &vec_mass,  32000, 0);
	tree->Branch("charge", "std::vector<int>",   &vec_charge,32000, 0);
	tree->Branch("pid",    "std::vector<int>",   &vec_pid,   32000, 0);

	// ---------------------------------------------------------------------
	// Parsing variables
	// ---------------------------------------------------------------------
	int   t_id = -1, t_pid = 0;
	int   t_isStarInt = 0, t_isSPhenixInt = 0;
	float t_pt = 0, t_eta = 0, t_phi = 0, t_mass = 0, t_charge = 0;

	int   currentID = -1;
	long  lineNumber = 0;
	long  nGoodLines = 0;
	long  nEvents    = 0;

	std::string line;

	// ---------------------------------------------------------------------
	// Main read loop: one line at a time, with checks
	// ---------------------------------------------------------------------
	while (std::getline(infile, line)) {
		++lineNumber;

		if (line.empty()) continue;
		if (line[0] == '#') continue;  // in case logs ever get comments

		std::istringstream iss(line);
		if (!(iss >> t_id >> t_isStarInt >> t_isSPhenixInt
					>> t_pt >> t_eta >> t_phi >> t_mass >> t_charge >> t_pid)) {
			std::cerr << "Warning: bad line " << lineNumber
				<< " in " << inFileName << ":\n"
				<< line << std::endl;
			continue;
		}

		++nGoodLines;

		// New event?
		if (t_id != currentID) {
			if (currentID != -1) {
				b_eventID = currentID;
				b_nTracks = static_cast<int>(vec_pt.size());
				tree->Fill();
				++nEvents;
			}

			currentID   = t_id;
			b_isStar    = (t_isStarInt    != 0);
			b_isSPhenix = (t_isSPhenixInt != 0);

			vec_pt.clear();
			vec_eta.clear();
			vec_phi.clear();
			vec_mass.clear();
			vec_charge.clear();
			vec_pid.clear();
		}

		vec_pt.push_back(t_pt);
		vec_eta.push_back(t_eta);
		vec_phi.push_back(t_phi);
		vec_mass.push_back(t_mass);
		vec_charge.push_back(static_cast<int>(t_charge));
		vec_pid.push_back(t_pid);

		// Optional: debug first few lines
		if (nGoodLines <= 5) {
			std::cout << "Debug: line " << lineNumber
				<< " ID=" << t_id
				<< " pt=" << t_pt
				<< " eta=" << t_eta
				<< " pid=" << t_pid << std::endl;
		}
	}

	// Flush last event
	if (currentID != -1) {
		b_eventID = currentID;
		b_nTracks = static_cast<int>(vec_pt.size());
		tree->Fill();
		++nEvents;
	} else {
		std::cerr << "Warning: no valid events found in " << inFileName
			<< " (nGoodLines = " << nGoodLines << ")" << std::endl;
	}

	// ---------------------------------------------------------------------
	// Finalize
	// ---------------------------------------------------------------------
	tree->Write();
	fOut->Close();
	infile.close();

	std::cout << "log_to_root: finished, wrote tree to "
		<< outFileName << " with " << nEvents
		<< " events from " << nGoodLines << " valid lines" << std::endl;
}

