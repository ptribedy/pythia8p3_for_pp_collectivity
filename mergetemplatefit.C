#include <iostream>
#include <fstream>
#include <cstring>   // <-- NEW
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

using namespace std;

TH1D * dphi;
TH1D * dphi_nonflow;
TH1D * dphi_nonflow_plot;

TF1 * templatefit;

// ------------------------------
// NEW: map pairTag -> histogram suffix
// pairTag: "all" (or ""), "os", "ss"
// ------------------------------
const char* PairSuffix(const char* pairTag)
{
    if (!pairTag) return "";
    if (strcmp(pairTag, "") == 0)    return "";
    if (strcmp(pairTag, "all") == 0) return "";
    if (strcmp(pairTag, "os")  == 0) return "_os";
    if (strcmp(pairTag, "ss")  == 0) return "_ss";
    // default fallback
    return "";
}

Double_t myfunction(Double_t *x, Double_t *par)
{
    Double_t f = par[4]*par[4]*templatefit->Eval(x[0]) +
                 par[0]*par[0]*(1. + par[1]*2.*cos(4.*x[0]) +
                                par[2]*2.*cos(2.*x[0]) +
                                par[3]*2.*cos(3.*x[0]));
    return f;
}

Double_t myfunction_y2ridge(Double_t *x, Double_t *par)
{
    Float_t xx = 0.;
    TAxis * xaxis = dphi_nonflow->GetXaxis();
    Int_t binx = xaxis->FindBin(xx);
    Double_t f = par[4]*par[4]*dphi_nonflow->GetBinContent(binx) +
                 par[0]*par[0]*(1. + par[1]*2.*cos(4.*x[0]) +
                                par[2]*2.*cos(2.*x[0]) +
                                par[3]*2.*cos(3.*x[0]));
    return f;
}

Double_t myfunction_noridge(Double_t *x, Double_t *par)
{
    Double_t f = par[4]*par[4]*templatefit->Eval(x[0]) +
                 par[0]*par[0]*(1. + par[1]*2.*cos(4.*x[0]) +
                                par[2]*2.*cos(2.*x[0]) +
                                par[3]*2.*cos(3.*x[0]));
    return f;
}

void myfunc()
{
    TF1 *f1 = new TF1("myfunc", myfunction, 0, TMath::Pi(), 5);
    f1->SetParameters(0.9, 0.1, 0.01, 0.01, 0.9);
}

void myfunc_y2ridge()
{
    TF1 *f2 = new TF1("myfunc_y2ridge", myfunction_y2ridge,
                      -TMath::Pi()/2., 3./2.*TMath::Pi(), 5);
    f2->SetParameters(0.9, 0.1, 0.01, 0.01, 0.9);
}

void myfunc_y3ridge()
{
    TF1 *f2 = new TF1("myfunc_y3ridge", myfunction_y2ridge,
                      -TMath::Pi()/2., 3./2.*TMath::Pi(), 5);
    f2->SetParameters(0.9, 0.1, 0.01, 0.01, 0.9);
}

void myfunc_noridge()
{
    TF1 *f2 = new TF1("myfunc_noridge", myfunction_noridge,
                      -TMath::Pi()/2., 3./2.*TMath::Pi(), 5);
    f2->SetParameters(0.9, 0.1, 0.01, 0.01, 0.9);
}

double getNtrig(int cent=6, int ptbinmin=3, int ptbinmax=3, const char* rootFileName="")
{
    TFile * f = new TFile(rootFileName);
    TH2D *Ntrig_detadphi = nullptr;
    char hname[256];

    snprintf(hname, sizeof(hname), "Ntrig_detadphi");
    Ntrig_detadphi = (TH2D*)f->Get(hname);
    if (!Ntrig_detadphi) {
        cerr << "ERROR: cannot find histogram " << hname << " in " << rootFileName << endl;
        f->Close();
        return 0.0;
    }

    double Ntrigger = Ntrig_detadphi->GetBinContent(ptbinmin+1, cent+1);

    for (int i_ = ptbinmin+1; i_ <= ptbinmax; i_++) {
        Ntrigger += Ntrig_detadphi->GetBinContent(i_+1, cent+1);
    }

    f->Close();
    return Ntrigger;
}

// ------------------------------
// NEW: pairTag added
// ------------------------------

TH1D* getfunc_mod(int cent=6, int ptbinmin=3, int ptbinmax=3,
                  const char* rootFileName="",
                  Double_t deltaetaMin = 1.0, Double_t deltaetaMax = 3.0,
                  const char* pairTag="all")
{
    cout << " I am here "
         << " ptbin " << ptbinmin << " ptbinmax " << ptbinmax
         << " cent= " << cent
         << " pairTag= " << (pairTag ? pairTag : "null")
         << endl;

    TFile *f = TFile::Open(rootFileName, "READ");
    if (!f || f->IsZombie()) {
        cerr << "ERROR: cannot open file " << rootFileName << endl;
        return nullptr;
    }

    const char* suf = PairSuffix(pairTag);
    char hname[256];

    // --- Load first histogram and CLONE it so it survives file close ---
    snprintf(hname, sizeof(hname), "detadphi%s_%d_%d", suf, cent, ptbinmin);
    TH2D *h0 = (TH2D*)f->Get(hname);
    if (!h0) {
        cerr << "ERROR: cannot find histogram " << hname << " in " << rootFileName << endl;
        f->Close();
        return nullptr;
    }

    TH2D *detadphi = (TH2D*)h0->Clone(Form("detadphi_work_%s_c%d_p%d", pairTag, cent, ptbinmin));
    detadphi->SetDirectory(nullptr);  // detach from file ownership
    detadphi->Sumw2();

    // --- Add other pT bins (clone/add safely) ---
    for (int i_ = ptbinmin+1; i_ <= ptbinmax; i_++) {
        snprintf(hname, sizeof(hname), "detadphi%s_%d_%d", suf, cent, i_);
        TH2D *ht = (TH2D*)f->Get(hname);
        if (!ht) {
            cerr << "WARNING: missing " << hname << " (skipping)" << endl;
            continue;
        }
        detadphi->Add(ht);
    }

    // Now safe to close file (we cloned what we need)
    f->Close();

    // --- Find Δη bin range ---
    TAxis* xax = detadphi->GetXaxis();
    Int_t binMin = xax->FindBin(deltaetaMin);
    Int_t binMax = xax->FindBin(deltaetaMax);

    std::cout
        << "Projecting Y over #Delta#Eta ["
        << xax->GetBinLowEdge(binMin) << ", "
        << xax->GetBinUpEdge(binMax)  << "] -> bins "
        << binMin << "_" << binMax
        << std::endl;

    // --- Unique projection name to avoid ROOT object collisions ---
    TH1D* dphi_1d = detadphi->ProjectionY(
        Form("dphi_1d_%s_c%d_pt%d_%d_d%g_%g", pairTag, cent, ptbinmin, ptbinmax, deltaetaMin, deltaetaMax),
        binMin, binMax
    );
    dphi_1d->SetDirectory(nullptr);
    dphi_1d->Sumw2();

    // --- Build output histogram (unique name) ---
    TH1D *detadphi_ratio_1d = new TH1D(
        Form("detadphi_ratio_1d_%s_c%d_pt%d_%d_d%g_%g", pairTag, cent, ptbinmin, ptbinmax, deltaetaMin, deltaetaMax),
        "",
        24, -0.5*TMath::Pi(), 1.5*TMath::Pi()
    );
    detadphi_ratio_1d->SetDirectory(nullptr);
    detadphi_ratio_1d->Sumw2();

    // Ntrig (still reads from file, but in a separate open/close)
    Double_t ntrig = getNtrig(cent, ptbinmin, ptbinmax, rootFileName);
    Double_t xbin  = 12./TMath::Pi();

    for (int j=1; j<=detadphi_ratio_1d->GetNbinsX(); j++) {
        Double_t num     = dphi_1d->GetBinContent(j);
        Double_t num_err = dphi_1d->GetBinError(j);

        Double_t ratio     = num;
        Double_t ratio_err = num_err;

        ratio     = (ntrig != 0 ? ratio/ntrig*xbin     : 0.);
        ratio_err = (ntrig != 0 ? ratio_err/ntrig*xbin : 0.);

        detadphi_ratio_1d->SetBinContent(j, ratio);
        detadphi_ratio_1d->SetBinError(j, ratio_err);
    }

    // cleanup intermediate
    delete detadphi;
    delete dphi_1d;

    return detadphi_ratio_1d;
}
// ------------------------------
// NEW: pairTag added
// ------------------------------
void mergetemplatefit(int cent=6, int ptbinmin=0, int ptbinmax=0,
                      const char* rootFileName="",
                      Double_t deltaetaMin = 1.0, Double_t deltaetaMax = 3.0,
                      const char* pairTag="all")
{
    TCanvas *c1 = new TCanvas("c1", "c1", 10, 67, 753, 502);

    char hname[256];

    dphi = getfunc_mod(cent, ptbinmin, ptbinmax, rootFileName, deltaetaMin, deltaetaMax, pairTag);
    dphi_nonflow = getfunc_mod(8, ptbinmin, ptbinmax, rootFileName, deltaetaMin, deltaetaMax, pairTag);
    dphi_nonflow_plot = getfunc_mod(8, ptbinmin, ptbinmax, rootFileName, deltaetaMin, deltaetaMax, pairTag);

    if (!dphi || !dphi_nonflow || !dphi_nonflow_plot) {
        cerr << "ERROR: null histogram returned (check inputs / histogram names)" << endl;
        return;
    }

    dphi->Sumw2();
    dphi_nonflow->Sumw2();
    dphi_nonflow_plot->Sumw2();

    double scale = 1.;

    TF1 *fa = new TF1("fa",
                      "[0] + [1]*[0]*2*cos(x)+ [2]*[0]*2*cos(2*x) + [3]*[0]*2*cos(3*x)",
                      -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fa->SetParameter(0,0.1);
    fa->SetParameter(1,0.01);
    fa->SetParameter(2,0.01);
    fa->SetParameter(3,0.01);

    dphi_nonflow->Fit("fa");
    templatefit = fa;

    myfunc();
    myfunc_y2ridge();
    myfunc_y3ridge();
    myfunc_noridge();

    TF1 *f1 = (TF1 *)gROOT->GetFunction("myfunc");
    f1->SetParameter(0,0.9);
    f1->SetParameter(1,0.01);
    f1->SetParameter(2,0.01);
    f1->SetParameter(3,0.01);
    f1->SetParameter(4,0.9);

    dphi->Fit(f1);

    TF1 *fa_HM = new TF1("fa_HM",
                         "[0] + [1]*[0]*2*cos(x)+ [2]*[0]*2*cos(2*x) + [3]*[0]*2*cos(3*x)",
                         -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fa_HM->SetParameter(0,0.1);
    fa_HM->SetParameter(1,0.01);
    fa_HM->SetParameter(2,0.01);
    fa_HM->SetParameter(3,0.01);
    dphi->Fit(fa_HM,"R");

    dphi->SetLineColor(kBlack);
    dphi->SetMarkerStyle(20);
    dphi->GetXaxis()->SetTitle("#Delta#phi");
    dphi->GetYaxis()->SetTitle("Y(#Delta #phi)");
    dphi->SetStats(0);
    dphi->SetTitle("p+p 200 GeV (Run-24 FastOffline)");
    dphi->Draw("9pze");

    TF1 *fconst = new TF1("fconst","[0]*[0]",-TMath::Pi()/2.,3./2.*TMath::Pi());
    fconst->SetParameter(0, f1->GetParameter(0));
    Double_t scalefact = f1->GetParameter(4)*f1->GetParameter(4);
    dphi_nonflow_plot->Scale(scalefact);
    dphi_nonflow_plot->Add(fconst);
    dphi_nonflow_plot->SetMarkerStyle(4);
    dphi_nonflow_plot->Draw("9pzsame");

    TF1 *y2ridge = (TF1 *)gROOT->GetFunction("myfunc_y2ridge");
    y2ridge->SetParameter(0,f1->GetParameter(0));
    y2ridge->SetParameter(1,0.);
    y2ridge->SetParameter(2,f1->GetParameter(2));
    y2ridge->SetParameter(3,0.);
    y2ridge->SetParameter(4,f1->GetParameter(4));
    y2ridge->SetLineColor(kBlue-3);
    y2ridge->SetLineStyle(4);
    y2ridge->Draw("same");

    TF1 *y3ridge = (TF1 *)gROOT->GetFunction("myfunc_y3ridge");
    y3ridge->SetParameter(0,f1->GetParameter(0));
    y3ridge->SetParameter(1,0.);
    y3ridge->SetParameter(2,0.);
    y3ridge->SetParameter(3,f1->GetParameter(3));
    y3ridge->SetParameter(4,f1->GetParameter(4));
    y3ridge->SetLineColor(kOrange+2);
    y3ridge->SetLineStyle(3);
    y3ridge->Draw("same");

    TF1 *noridge = (TF1 *)gROOT->GetFunction("myfunc_noridge");
    noridge->SetNpx(1000);
    noridge->SetParameter(0,f1->GetParameter(0));
    noridge->SetParameter(1,0.);
    noridge->SetParameter(2,0.);
    noridge->SetParameter(3,0.);
    noridge->SetParameter(4,f1->GetParameter(4));
    noridge->SetLineColor(kGreen+2);
    noridge->SetLineStyle(7);
    noridge->Draw("same");

    auto legend = new TLegend(0.138482,0.7463312,0.517976,0.9475891,NULL,"brNDC");
    legend->AddEntry(dphi, Form("HM (%s)", pairTag), "ep");
    legend->AddEntry(dphi_nonflow_plot, Form("LM (%s)", pairTag), "ep");
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.0524109);
    legend->Draw();

    // ---- output files: include pairTag to avoid overwrite ----
    ofstream datafile;
    snprintf(hname, sizeof(hname),
             "data_dndphi_%s_%d_%d_%d_%g_%g.txt",
             pairTag, cent, ptbinmin, ptbinmax, deltaetaMin, deltaetaMax);
    datafile.open(hname);

    for (int i=1; i<=dphi->GetNbinsX(); i++){
        datafile << dphi->GetBinCenter(i) << " "
                 << dphi->GetBinContent(i)/scale << " "
                 << dphi->GetBinError(i)/scale << " "
                 << dphi_nonflow_plot->GetBinContent(i)/scale << " "
                 << dphi_nonflow_plot->GetBinError(i)/scale << " "
                 << y2ridge->Eval(dphi->GetBinCenter(i))/scale << " "
                 << y3ridge->Eval(dphi->GetBinCenter(i))/scale << " "
                 << noridge->Eval(dphi->GetBinCenter(i))/scale << " "
                 << f1->Eval(dphi->GetBinCenter(i))/scale << " "
                 << fa->Eval(dphi->GetBinCenter(i))/scale
                 << endl;
    }
    datafile.close();

    Double_t v22    = f1->GetParameter(2);
    Double_t v32    = f1->GetParameter(3);
    Double_t v22err = f1->GetParError(2);
    Double_t v32err = f1->GetParError(3);

    Double_t v2 = (v22>0 ? sqrt(v22) : 0.);
    Double_t v3 = (v32>0 ? sqrt(v32) : 0.);
    Double_t v2err = (v22>0 ? 0.5/sqrt(v22)*v22err : 0.);
    Double_t v3err = (v32>0 ? 0.5/sqrt(v32)*v32err : 0.);

    cout << " pairTag=" << pairTag
         << " v22=" << v22 << " ± " << v22err
         << " v32=" << v32 << " ± " << v32err
         << " v2="  << v2  << " ± " << v2err
         << " v3="  << v3  << " ± " << v3err
         << endl;

    Double_t v1_HM     = fa_HM->GetParameter(1);
    Double_t v2_HM     = fa_HM->GetParameter(2);
    Double_t v3_HM     = fa_HM->GetParameter(3);
    Double_t v1_HM_err = fa_HM->GetParError(1);
    Double_t v2_HM_err = fa_HM->GetParError(2);
    Double_t v3_HM_err = fa_HM->GetParError(3);

    ofstream fitfile;
    snprintf(hname, sizeof(hname),
             "fit_dndphi_%s_%d_%d_%d_%g_%g.txt",
             pairTag, cent, ptbinmin, ptbinmax, deltaetaMin, deltaetaMax);
    fitfile.open(hname);

    for(int i=1; i<=dphi->GetNbinsX(); i++){
        Double_t phi_low  = dphi->GetBinLowEdge(i);
        Double_t phi_cent = dphi->GetBinCenter(i);

        auto dumpLine = [&](double phi){
            fitfile << phi
                    << " " << dphi->GetBinContent(i)/scale
                    << " " << dphi->GetBinError(i)/scale
                    << " " << dphi_nonflow_plot->GetBinContent(i)/scale
                    << " " << dphi_nonflow_plot->GetBinError(i)/scale
                    << " " << y2ridge->Eval(phi)/scale
                    << " " << y3ridge->Eval(phi)/scale
                    << " " << noridge->Eval(phi)/scale
                    << " " << f1->Eval(phi)/scale
                    << " " << fa->Eval(phi)/scale
                    << " " << v22
                    << " " << v22err
                    << " " << v32
                    << " " << v32err
                    << " " << v2
                    << " " << v2err
                    << " " << v3
                    << " " << v3err
                    << " " << v1_HM
                    << " " << v1_HM_err
                    << " " << v2_HM
                    << " " << v2_HM_err
                    << " " << v3_HM
                    << " " << v3_HM_err
                    << endl;
        };

        dumpLine(phi_low);
        dumpLine(phi_cent);
    }
    fitfile.close();
}
