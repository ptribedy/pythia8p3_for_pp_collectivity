#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TMath.h"

#include <fstream>
#include "TObjArray.h"
#include "TObjString.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ----------------------------------------------------------------------
// Global binning
// ----------------------------------------------------------------------
const int nCentBins = 9;   // 0 .. 8 (0 = highest mult)
const int nPtBins   = 9;   // 0 .. 8

const double ptMin = 0.2;
const double ptMax = 2.0;

// Sub-event eta gap: we use TPC East eta < -0.65 and West eta > +0.65
// so effective gap between sub-events is 1.3 in eta
const double etaGapSub = 0.65;

// ---------------------
// Centrality: itpcMult in [0, 60]
// 8 = lowest mult, 0 = highest mult
// ---------------------
int GetCentrality(int mult)
{
	if (mult < 0)  return -1;          // invalid
	if (mult <= 6) return 8;           // most peripheral
	if (mult <= 12) return 7;
	if (mult <= 18) return 6;
	if (mult <= 24) return 5;
	if (mult <= 30) return 4;
	if (mult <= 36) return 3;
	if (mult <= 42) return 2;
	if (mult <= 48) return 1;
	if (mult <= 60) return 0;          // most central
	return 0;                          // mult > 60: lump into most central
}

// ---------------------
// pT binning: uniform in [ptMin, ptMax)
// ---------------------
int GetPtBin(double pt)
{
	if (pt < ptMin || pt >= ptMax) return -1;
	double width = (ptMax - ptMin) / nPtBins;
	int bin = static_cast<int>((pt - ptMin) / width);
	if (bin < 0 || bin >= nPtBins) return -1;
	return bin;
}

// ==========================================================================
//                                 MAIN
// ==========================================================================

void flowanalysis(TString inputFiles, TString outputFile)
{
	// ------------------------------------------------------------------
	// Build TChain from either:
	//  1) a .list file (one ROOT path per line), or
	//  2) a space-separated list of ROOT paths
	// ------------------------------------------------------------------
	TChain *chain = new TChain("tree");

	if (inputFiles.EndsWith(".list")) {
		cout << "Reading input file list from: " << inputFiles << endl;
		ifstream fin(inputFiles.Data());
		if (!fin.is_open()) {
			cout << "Error: cannot open list file " << inputFiles << endl;
			return;
		}
		std::string line;
		while (std::getline(fin, line)) {
			TString fname(line.c_str());
			fname = fname.Strip(TString::kBoth, ' ');
			if (fname.Length() == 0) continue;
			cout << "  " << fname << endl;
			chain->Add(fname);
		}
		fin.close();
	} else {
		cout << "Adding input files to chain from space-separated string:" << endl;
		cout << "  [" << inputFiles << "]" << endl;

		TObjArray *tokens = inputFiles.Tokenize(" ");
		int nTokens = tokens ? tokens->GetEntriesFast() : 0;

		for (int i = 0; i < nTokens; ++i) {
			TObjString *os = dynamic_cast<TObjString*>(tokens->At(i));
			if (!os) continue;
			TString fname = os->GetString();
			fname = fname.Strip(TString::kBoth, ' ');
			if (fname.Length() == 0) continue;
			cout << "  " << fname << endl;
			chain->Add(fname);
		}
		delete tokens;
	}

	if (!chain->GetEntries()) {
		cout << "Error: No entries found in any input files!" << endl;
		return;
	}

	cout << "Total events in chain: " << chain->GetEntries() << endl;

	// ------------------------------------------------------------------
	// Branch setup
	// ------------------------------------------------------------------
	vector<float> *v_pt     = 0;
	vector<float> *v_eta    = 0;
	vector<float> *v_phi    = 0;
	vector<int>   *v_charge = 0;
	int nTracks = 0;

	chain->SetBranchAddress("nTracks", &nTracks);
	chain->SetBranchAddress("pT",     &v_pt);
	chain->SetBranchAddress("eta",    &v_eta);
	chain->SetBranchAddress("phi",    &v_phi);
	chain->SetBranchAddress("charge", &v_charge);

	TFile *fOut = new TFile(outputFile, "RECREATE");

	// ================= QA HISTS ===================
	TH1D *hRefMult = new TH1D("hRefMult",
			"RefMult;Multiplicity;",
			100, -0.5, 99.5);

	TH1D *hTofMult = new TH1D("hTofMult",
			"TOF Mult;Multiplicity;",
			100, -0.5, 99.5);

	// ================= Δη–Δφ HISTS (cent x pT bin, data-like) ===============
	TH2D *detadphi   [nCentBins][nPtBins];
	TH2D *detadphi_os[nCentBins][nPtBins];
	TH2D *detadphi_ss[nCentBins][nPtBins];

	for (int c = 0; c < nCentBins; ++c) {
		for (int p = 0; p < nPtBins; ++p) {

			detadphi[c][p] = new TH2D(
					Form("detadphi_%d_%d", c, p),
					Form("detadphi_%d_%d;#Delta#eta;#Delta#phi", c, p),
					60, 0.0, 3.0,
					24, -M_PI/2, 3*M_PI/2
					);

			detadphi_os[c][p] = new TH2D(
					Form("detadphi_os_%d_%d", c, p),
					Form("detadphi_os_%d_%d;#Delta#eta;#Delta#phi", c, p),
					60, 0.0, 3.0,
					24, -M_PI/2, 3*M_PI/2
					);

			detadphi_ss[c][p] = new TH2D(
					Form("detadphi_ss_%d_%d", c, p),
					Form("detadphi_ss_%d_%d;#Delta#eta;#Delta#phi", c, p),
					60, 0.0, 3.0,
					24, -M_PI/2, 3*M_PI/2
					);

			detadphi   [c][p]->Sumw2();
			detadphi_os[c][p]->Sumw2();
			detadphi_ss[c][p]->Sumw2();
		}
	}

	// ================= TRIGGER COUNTERS ===============
	TH2D *Ntrig_detadphi = new TH2D(
			"Ntrig_detadphi",
			"N_{trig};p_{T} bin index;Centrality;",
			nPtBins, -0.5, nPtBins - 0.5,
			nCentBins, -0.5, nCentBins - 0.5
			);

	TH1D *Ntrig_detadphi_eventcount = new TH1D(
			"Ntrig_detadphi_eventcount",
			"Event count;Centrality;",
			nCentBins, -0.5, nCentBins - 0.5
			);

	// ================= v_n2 PROFILES (all, OS, SS) ===============
	TProfile *v12tpc     = new TProfile("v12tpc",
			"v_{1}^{2} (TPC all);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);
	TProfile *v22tpc     = new TProfile("v22tpc",
			"v_{2}^{2} (TPC all);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);
	TProfile *v32tpc     = new TProfile("v32tpc",
			"v_{3}^{2} (TPC all);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);

	TProfile *v12tpc_os  = new TProfile("v12tpc_os",
			"v_{1}^{2} (TPC OS);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);
	TProfile *v22tpc_os  = new TProfile("v22tpc_os",
			"v_{2}^{2} (TPC OS);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);
	TProfile *v32tpc_os  = new TProfile("v32tpc_os",
			"v_{3}^{2} (TPC OS);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);

	TProfile *v12tpc_ss  = new TProfile("v12tpc_ss",
			"v_{1}^{2} (TPC SS);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);
	TProfile *v22tpc_ss  = new TProfile("v22tpc_ss",
			"v_{2}^{2} (TPC SS);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);
	TProfile *v32tpc_ss  = new TProfile("v32tpc_ss",
			"v_{3}^{2} (TPC SS);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);

	v12tpc->Sumw2();    v22tpc->Sumw2();    v32tpc->Sumw2();
	v12tpc_os->Sumw2(); v22tpc_os->Sumw2(); v32tpc_os->Sumw2();
	v12tpc_ss->Sumw2(); v22tpc_ss->Sumw2(); v32tpc_ss->Sumw2();

	// ============== SUB-EVENT v_n2 PROFILES (TPC East–West) ==============
	TProfile *v1tpc_subevent = new TProfile("v1tpc_subevent",
			"v_{1}^{2} (TPC subevent E–W);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);

	TProfile *v2tpc_subevent = new TProfile("v2tpc_subevent",
			"v_{2}^{2} (TPC subevent E–W);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);

	TProfile *v3tpc_subevent = new TProfile("v3tpc_subevent",
			"v_{3}^{2} (TPC subevent E–W);Centrality;",
			nCentBins, -0.5, nCentBins - 0.5);

	v1tpc_subevent->Sumw2();
	v2tpc_subevent->Sumw2();
	v3tpc_subevent->Sumw2();

	// ======================================================================
	//                              EVENT LOOP
	// ======================================================================
	Long64_t nEvents = chain->GetEntries();

	for (Long64_t iev = 0; iev < nEvents; iev++) {
		chain->GetEntry(iev);

		int itpcMult = 0;
		int bTOF     = 0;

		if (!v_pt || v_pt->size() != (unsigned)nTracks) continue;

		// ---------- SUB-EVENT Q-vectors (TPC East/West) ----------
		double Q0tpcE  = 0.0;
		double Q1xtpcE = 0.0, Q1ytpcE = 0.0;
		double Q2xtpcE = 0.0, Q2ytpcE = 0.0;
		double Q3xtpcE = 0.0, Q3ytpcE = 0.0;

		double Q0tpcW  = 0.0;
		double Q1xtpcW = 0.0, Q1ytpcW = 0.0;
		double Q2xtpcW = 0.0, Q2ytpcW = 0.0;
		double Q3xtpcW = 0.0, Q3ytpcW = 0.0;

		// multiplicity for centrality, QA, and build sub-event Q-vectors
		for (int i = 0; i < nTracks; i++) {
			double pt   = v_pt->at(i);
			double eta  = v_eta->at(i);
			double phi  = v_phi->at(i);

			if (pt < 0.2) continue;
			if (fabs(eta) < 1.5) itpcMult++;    // *** max |eta|=1.5 for itpcMult ***
			if (fabs(eta) < 0.9) bTOF++;

			// restrict Q-vectors to same kinematic as correlations
			if (pt < ptMin || pt >= ptMax || fabs(eta) > 1.5) continue;

			// East: eta < -etaGapSub
			if (eta < -etaGapSub) {
				double w = 1.0;
				Q0tpcE  += w;
				Q1xtpcE += cos(1.0 * phi) * w;
				Q1ytpcE += sin(1.0 * phi) * w;
				Q2xtpcE += cos(2.0 * phi) * w;
				Q2ytpcE += sin(2.0 * phi) * w;
				Q3xtpcE += cos(3.0 * phi) * w;
				Q3ytpcE += sin(3.0 * phi) * w;
			}

			// West: eta > +etaGapSub
			if (eta >  etaGapSub) {
				double w = 1.0;
				Q0tpcW  += w;
				Q1xtpcW += cos(1.0 * phi) * w;
				Q1ytpcW += sin(1.0 * phi) * w;
				Q2xtpcW += cos(2.0 * phi) * w;
				Q2ytpcW += sin(2.0 * phi) * w;
				Q3xtpcW += cos(3.0 * phi) * w;
				Q3ytpcW += sin(3.0 * phi) * w;
			}
		}

		hRefMult->Fill(itpcMult);
		hTofMult->Fill(bTOF);

		int centBin = GetCentrality(itpcMult);

		if (iev % 1000 == 0)
			cout << " event " << iev
				<< " nTracks= " << nTracks
				<< " itpcMult= " << itpcMult
				<< " centBin= " << centBin << endl;

		if (centBin < 0 || centBin >= nCentBins) continue;

		Ntrig_detadphi_eventcount->Fill(centBin);

		// ------------- v_n2 accumulators for THIS EVENT (all, OS, SS) ------
		double sumCos1_all = 0, sumCos2_all = 0, sumCos3_all = 0;
		double sumCos1_os  = 0, sumCos2_os  = 0, sumCos3_os  = 0;
		double sumCos1_ss  = 0, sumCos2_ss  = 0, sumCos3_ss  = 0;

		int countLR_all = 0;
		int countLR_os  = 0;
		int countLR_ss  = 0;

		// Trigger–associate loops
		for (int i = 0; i < nTracks; i++) {

			double pti  = v_pt->at(i);
			double etai = v_eta->at(i);

			// *** max |eta| acceptance 1.5 now ***
			if (pti < ptMin || pti >= ptMax || fabs(etai) > 1.5) continue;

			int ptiBin = GetPtBin(pti);
			if (ptiBin < 0 || ptiBin >= nPtBins) continue;

			// count triggers per (pT bin, centrality)
			Ntrig_detadphi->Fill(ptiBin, centBin);

			for (int j = 0; j < nTracks; j++) {
				if (i == j) continue;

				double ptj  = v_pt->at(j);
				double etaj = v_eta->at(j);
				if (ptj < ptMin || ptj >= ptMax || fabs(etaj) > 1.5) continue;

				double dEta = etai - etaj;
				double dPhi = v_phi->at(i) - v_phi->at(j);

				if (dPhi < -M_PI/2)       dPhi += 2*M_PI;
				else if (dPhi > 3*M_PI/2) dPhi -= 2*M_PI;

				detadphi[centBin][ptiBin]->Fill(dEta, dPhi);

				bool isOS = (v_charge->at(i) != v_charge->at(j));
				if (isOS)
					detadphi_os[centBin][ptiBin]->Fill(dEta, dPhi);
				else
					detadphi_ss[centBin][ptiBin]->Fill(dEta, dPhi);

				// Long-range region for flow
				if (fabs(dEta) > 1.3) {
					double c1 = cos(1.0 * dPhi);
					double c2 = cos(2.0 * dPhi);
					double c3 = cos(3.0 * dPhi);

					// all pairs
					sumCos1_all += c1;
					sumCos2_all += c2;
					sumCos3_all += c3;
					countLR_all++;

					// OS / SS split
					if (isOS) {
						sumCos1_os += c1;
						sumCos2_os += c2;
						sumCos3_os += c3;
						countLR_os++;
					} else {
						sumCos1_ss += c1;
						sumCos2_ss += c2;
						sumCos3_ss += c3;
						countLR_ss++;
					}
				}
			}
		}

		// -------- ALL CHARGES --------
		double v12tpc_evt = 0.0, v22tpc_evt = 0.0, v32tpc_evt = 0.0;
		if (countLR_all > 0) {

			v12tpc_evt = sumCos1_all / countLR_all;
			v22tpc_evt = sumCos2_all / countLR_all;
			v32tpc_evt = sumCos3_all / countLR_all;

			v12tpc->Fill(centBin, v12tpc_evt);
			v22tpc->Fill(centBin, v22tpc_evt);
			v32tpc->Fill(centBin, v32tpc_evt);
		}

		// -------- OPPOSITE SIGN --------
		if (countLR_os > 0) {

			double v12tpc_evt_os = sumCos1_os / countLR_os;
			double v22tpc_evt_os = sumCos2_os / countLR_os;
			double v32tpc_evt_os = sumCos3_os / countLR_os;

			v12tpc_os->Fill(centBin, v12tpc_evt_os);
			v22tpc_os->Fill(centBin, v22tpc_evt_os);
			v32tpc_os->Fill(centBin, v32tpc_evt_os);
		}

		// -------- SAME SIGN --------
		if (countLR_ss > 0) {

			double v12tpc_evt_ss = sumCos1_ss / countLR_ss;
			double v22tpc_evt_ss = sumCos2_ss / countLR_ss;
			double v32tpc_evt_ss = sumCos3_ss / countLR_ss;

			v12tpc_ss->Fill(centBin, v12tpc_evt_ss);
			v22tpc_ss->Fill(centBin, v22tpc_evt_ss);
			v32tpc_ss->Fill(centBin, v32tpc_evt_ss);
		}

		// ==================== SUB-EVENT METHOD (TPC E–W, |eta| gap 0.65) ====================
		if (Q0tpcE > 0.0 && Q0tpcW > 0.0) {
			double w_evt = Q0tpcE * Q0tpcW;   // STAR-like event weight

			// n = 1
			double v1_sub_evt = (Q1xtpcE * Q1xtpcW + Q1ytpcE * Q1ytpcW) / (Q0tpcE * Q0tpcW);
			v1tpc_subevent->Fill(centBin, v1_sub_evt, w_evt);

			// n = 2
			double v2_sub_evt = (Q2xtpcE * Q2xtpcW + Q2ytpcE * Q2ytpcW) / (Q0tpcE * Q0tpcW);
			v2tpc_subevent->Fill(centBin, v2_sub_evt, w_evt);

			// n = 3
			double v3_sub_evt = (Q3xtpcE * Q3xtpcW + Q3ytpcE * Q3ytpcW) / (Q0tpcE * Q0tpcW);
			v3tpc_subevent->Fill(centBin, v3_sub_evt, w_evt);
		}

		if (iev % 1000 == 0 && countLR_all > 0)
			cout << " event " << iev
				<< " v12, v22, v32 = "
				<< v12tpc_evt << " "
				<< v22tpc_evt << " "
				<< v32tpc_evt << endl;
	}

	fOut->Write();
	fOut->Close();
}
