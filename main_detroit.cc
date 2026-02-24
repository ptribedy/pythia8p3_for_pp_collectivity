// main_detroit.cc
// To compile: g++ main_detroit.cc -o main_detroit -O2 `root-config --cflags --glibs` -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"

using namespace Pythia8;
using namespace std;

int main(int argc, char* argv[]) {

	// --- 1. OUTPUT SETUP ---
	int seed = 10;
	if (argc > 1) seed = atoi(argv[1]);

	// Writes to current directory (./output_SEED.log)
	// No "output/" prefix to avoid folder creation issues
	string outName = "output_" + to_string(seed) + ".log";
	cout << "Writing to file: " << outName << endl;

	ofstream myFile;
	myFile.open(outName.c_str());

	// --- 2. PYTHIA CONFIGURATION ---
	Pythia pythia;

	pythia.readString("Beams:idA = 2212");
	pythia.readString("Beams:idB = 2212");
	pythia.readString("Beams:eCM = 200.");
	pythia.readString("SoftQCD:inelastic = on");
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = " + to_string(seed));

	// "Detroit" Tune Settings
	pythia.readString("PDF:pSet = 17");
	pythia.readString("MultipartonInteractions:bProfile = 2");
	pythia.readString("MultipartonInteractions:ecmRef = 200");
	pythia.readString("MultipartonInteractions:pT0Ref          = 1.40");
	pythia.readString("MultipartonInteractions:ecmPow          = 0.135");
	pythia.readString("MultipartonInteractions:coreRadius      = 0.56");
	pythia.readString("MultipartonInteractions:coreFraction    = 0.78");
	pythia.readString("ColourReconnection:range                = 5.4");

	if (!pythia.init()) return 1;

	Hist mult("charged multiplicity", 100, -0.5, 799.5);

	// --- 3. EVENT LOOP ---
	int nEvents = 100000; 

	for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
		if (!pythia.next()) continue;

		bool isStarF = false, isStarB = false, isStar = false;
		bool isSPhenixF = false, isSPhenixB = false, isSPhenix = false;
		int nCharged = 0;

		// Loop 1: Triggers & Mult
		for (int i = 0; i < pythia.event.size(); ++i) {
			Particle& p = pythia.event[i];
			if (!p.isFinal()) continue;
			if (p.isCharged() && fabs(p.eta()) < 1.0 && p.pT() > 0.1) nCharged++;

			if (p.pT() > 0.1) {
				double eta = p.eta();
				if (eta > 3.4 && eta < 5.0)   isStarF = true;
				if (eta > -5.0 && eta < -3.4) isStarB = true;
				if (eta > 3.1 && eta < 3.9)   isSPhenixF = true;
				if (eta > -3.9 && eta < -3.1) isSPhenixB = true;
			}
		}
		if (isStarF && isStarB)       isStar = true;
		if (isSPhenixF && isSPhenixB) isSPhenix = true;

		mult.fill(nCharged);

		// Loop 2: Output
		for (int i = 0; i < pythia.event.size(); ++i) {
			Particle& p = pythia.event[i];
			if (!p.isFinal() || !p.isCharged()) continue;

			double pt  = p.pT();
			double eta = p.eta();
			double phi = p.phi();
			double m   = p.m();
			int id     = p.id();
			double chg = p.charge();

			if (pt > 0.2 && m >= 0 && 
					(fabs(eta) < 1.5 || (fabs(eta) > 2.5 && fabs(eta) < 5.0))) {

				myFile << iEvent << " " << isStar << " " << isSPhenix << " " 
					<< pt << " " << eta << " " << phi << " " << m << " " 
					<< chg << " " << id << endl;
			}
		}
	}

	pythia.stat();
	cout << "Mult: " << mult;
	myFile.close();
	return 0;
}
