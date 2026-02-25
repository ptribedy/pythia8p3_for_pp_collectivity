// plot_closure_LR_vs_SR_pp200.C
// P. Tribedy, Feb 24, 2026
// this one is a closure test to check if nonflow subtraction is working or no
// Produces 5 canvases:
//   Canvas A (3 panels): cn(all,LR) | cn(OS-SS,LR) | cn(OS-SS,SR)
//                        panels 2&3 share global y-axis
//   Canvas B-LR (4 panels): cn(OS,LR) | cn(SS,LR) | cn(CI,LR) | cn(CD,LR)
//                            all 4 panels share global y-axis
//   Canvas B-SR (4 panels): cn(OS,SR) | cn(SS,SR) | cn(CI,SR) | cn(CD,SR)
//                            all 4 panels share global y-axis
//   Canvas C (3 panels): c1(all,LR) | c1(OS-SS,LR) | c1(OS-SS,SR)
//
// Usage:
//   root -l -q 'plot_closure_LR_vs_SR_pp200.C(2,7)'
//   root -l -q 'plot_closure_LR_vs_SR_pp200.C(3,8)'

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <map>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TColor.h"

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
// UTILITIES
// ─────────────────────────────────────────────────────────────────────────────

string binLabel(int b) {
    const char* lab[9]={"49-60","43-48","37-42","31-36","25-30","19-24","13-18","7-12","0-6"};
    return (b>=0&&b<=8) ? string(lab[b]) : "Unknown";
}

TGraphErrors* readVnFile(const string& fname, int col, int mkr, double xs=0.) {
    ifstream in(fname.c_str());
    if (!in.is_open()) { cerr<<"WARNING: cannot open "<<fname<<endl; return nullptr; }
    vector<double> vx,vy,vey; string line;
    while (getline(in,line)) {
        if (line.empty()||line[0]=='#') continue;
        stringstream ss(line); double v; vector<double> row;
        while (ss>>v) row.push_back(v);
        if (row.size()<7) continue;
        vx.push_back(row[6]+xs); vy.push_back(row[4]); vey.push_back(row[5]);
    }
    in.close(); if (vx.empty()) return nullptr;
    auto* g = new TGraphErrors(vx.size());
    for (size_t i=0;i<vx.size();++i) {
        g->SetPoint(i,vx[i],vy[i]); g->SetPointError(i,0.,vey[i]);
    }
    g->SetMarkerStyle(mkr); g->SetMarkerSize(1.5);
    g->SetMarkerColor(col);  g->SetLineColor(col); g->SetLineWidth(2);
    return g;
}

pair<double,double> getRefVal(const string& fname, int refBin) {
    ifstream in(fname.c_str()); pair<double,double> res={0.,0.};
    if (!in.is_open()) return res; string line;
    while (getline(in,line)) {
        if (line.empty()||line[0]=='#') continue;
        stringstream ss(line); double v; vector<double> row;
        while (ss>>v) row.push_back(v);
        if (row.size()<7) continue;
        if (fabs(row[6]-refBin)<0.01) { res={row[4],row[5]}; break; }
    }
    in.close(); return res;
}

// Y = A - (B*D)/C, with full error propagation
TGraphErrors* makeSubtracted(const string& vn_f, const string& v1_f,
                              int refBin, int col, int mkr) {
    auto vn_ref = getRefVal(vn_f, refBin);
    auto v1_ref = getRefVal(v1_f, refBin);
    double C=v1_ref.first, sC=v1_ref.second;
    double D=vn_ref.first, sD=vn_ref.second;
    if (fabs(C)<1e-9) { cerr<<"ERROR: c1(ref)~0 in "<<v1_f<<endl; return nullptr; }

    ifstream invn(vn_f.c_str()), inv1(v1_f.c_str());
    if (!invn.is_open()||!inv1.is_open()) return nullptr;

    map<int,pair<double,double>> v1m; string line;
    while (getline(inv1,line)) {
        if (line.empty()||line[0]=='#') continue;
        stringstream ss(line); vector<double> row;
        while (ss>>row.emplace_back());
        if (row.size()>=7) v1m[(int)row[6]]={row[4],row[5]};
    } inv1.close();

    vector<double> vx,vy,vey;
    while (getline(invn,line)) {
        if (line.empty()||line[0]=='#') continue;
        stringstream ss(line); double v; vector<double> row;
        while (ss>>v) row.push_back(v);
        if (row.size()<7) continue;
        int bin=(int)row[6]; double A=row[4], sA=row[5];
        if (v1m.find(bin)==v1m.end()) continue;
        double B=v1m[bin].first, sB=v1m[bin].second;
        double Y  = A-(B*D)/C;
        double sY = sqrt(sA*sA + pow((D/C)*sB,2) + pow((B/C)*sD,2) + pow((B*D)/(C*C)*sC,2));
        vx.push_back((double)bin); vy.push_back(Y); vey.push_back(sY);
    } invn.close();

    if (vx.empty()) return nullptr;
    auto* g = new TGraphErrors(vx.size());
    for (size_t i=0;i<vx.size();++i) {
        g->SetPoint(i,vx[i],vy[i]); g->SetPointError(i,0.,vey[i]);
    }
    g->SetMarkerStyle(mkr); g->SetMarkerSize(1.5);
    g->SetMarkerColor(col);  g->SetLineColor(col); g->SetLineWidth(2);
    return g;
}

// s1*g1 + s2*g2, matched by bin index
TGraphErrors* combineGraphs(TGraphErrors* g1, TGraphErrors* g2,
                             double s1, double s2, int col, int mkr) {
    if (!g1||!g2) return nullptr;
    map<int,pair<double,double>> m2;
    for (int i=0;i<g2->GetN();++i) {
        double x,y; g2->GetPoint(i,x,y);
        m2[(int)round(x)]={y,g2->GetErrorY(i)};
    }
    vector<double> vx,vy,vey;
    for (int i=0;i<g1->GetN();++i) {
        double x,y1; g1->GetPoint(i,x,y1); double e1=g1->GetErrorY(i);
        int key=(int)round(x);
        if (m2.find(key)==m2.end()) continue;
        double y2=m2[key].first, e2=m2[key].second;
        vx.push_back(x);
        vy.push_back(s1*y1+s2*y2);
        vey.push_back(sqrt(s1*s1*e1*e1+s2*s2*e2*e2));
    }
    if (vx.empty()) return nullptr;
    auto* g = new TGraphErrors(vx.size());
    for (size_t i=0;i<vx.size();++i) {
        g->SetPoint(i,vx[i],vy[i]); g->SetPointError(i,0.,vey[i]);
    }
    g->SetMarkerStyle(mkr); g->SetMarkerSize(1.5);
    g->SetMarkerColor(col);  g->SetLineColor(col); g->SetLineWidth(2);
    return g;
}

// Scan a list of graphs — returns [global_ymin, global_ymax] including error bars
pair<double,double> globalRange(vector<TGraphErrors*> gs) {
    double ymin=1e9, ymax=-1e9;
    for (auto* g : gs) {
        if (!g) continue;
        for (int i=0;i<g->GetN();++i) {
            double x,y; g->GetPoint(i,x,y); double e=g->GetErrorY(i);
            if (y-e < ymin) ymin = y-e;
            if (y+e > ymax) ymax = y+e;
        }
    }
    return {ymin, ymax};
}

TH1F* makeFrame(const char* nm, const char* yt) {
    auto* f = new TH1F(nm,"",9,-0.5,8.5); f->SetStats(0);
    const char* lab[9]={"49-60","43-48","37-42","31-36","25-30","19-24","13-18","7-12","0-6"};
    for (int i=1;i<=9;i++) f->GetXaxis()->SetBinLabel(i,lab[i-1]);
    f->GetXaxis()->SetBit(TAxis::kLabelsVert);
    f->GetXaxis()->CenterTitle(true);
    f->GetXaxis()->SetLabelFont(42); f->GetXaxis()->SetLabelSize(0.050);
    f->GetXaxis()->SetTitleFont(42); f->GetXaxis()->SetTitleSize(0.055);
    f->GetXaxis()->SetTitleOffset(1.4);
    f->GetXaxis()->SetTitle("Activity Class (Tracks |#eta|<1.5)");
    f->GetYaxis()->SetTitle(yt);
    f->GetYaxis()->SetLabelFont(42); f->GetYaxis()->SetLabelSize(0.040);
    f->GetYaxis()->SetTitleFont(42); f->GetYaxis()->SetTitleSize(0.052);
    f->GetYaxis()->SetTitleOffset(1.4);
    f->GetYaxis()->SetMaxDigits(3);
    f->GetYaxis()->SetNdivisions(3000510);
    return f;
}

// autoScale for single panels that keep their own independent range
void autoScale(TH1F* fr, vector<TGraphErrors*> gs, double lo=0.15, double hi=0.45) {
    auto [ymin,ymax] = globalRange(gs);
    double dy=ymax-ymin; if (dy<1e-9) dy=0.001;
    fr->SetMinimum(ymin-lo*dy); fr->SetMaximum(ymax+hi*dy);
}

void setPad(double lm=0.17, double rm=0.06, double bm=0.18, double tm=0.08) {
    gPad->SetLeftMargin(lm);  gPad->SetRightMargin(rm);
    gPad->SetBottomMargin(bm); gPad->SetTopMargin(tm);
}

TLine* zeroLine() {
    auto* z = new TLine(-0.5,0,8.5,0);
    z->SetLineStyle(2); z->SetLineColor(TColor::GetColor("#666666")); z->SetLineWidth(1);
    return z;
}

// ─────────────────────────────────────────────────────────────────────────────
// CANVAS A: cn(all,LR) | cn(OS-SS,LR) [3 approaches] | cn(OS-SS,SR) [3 approaches]
// Panel 1: independent scale
// Panels 2 & 3: shared global y-axis (all OS-SS graphs, both LR and SR)
// ─────────────────────────────────────────────────────────────────────────────
void drawCanvasA(int n, int refBin) {
    string v1_all  = "v12pc_fourierfit_all_LR";
    string v1_os   = "v12pc_fourierfit_os_LR";
    string v1_ss   = "v12pc_fourierfit_ss_LR";
    string v1_osss = "v12pc_fourierfit_OSminusSS_LR";

    string f_all_LR  = Form("v%d2pc_fourierfit_all_LR",        n);
    string f_os_LR   = Form("v%d2pc_fourierfit_os_LR",         n);
    string f_ss_LR   = Form("v%d2pc_fourierfit_ss_LR",         n);
    string f_osss_LR = Form("v%d2pc_fourierfit_OSminusSS_LR",  n);
    string f_os_SR   = Form("v%d2pc_fourierfit_os_SR",         n);
    string f_ss_SR   = Form("v%d2pc_fourierfit_ss_SR",         n);
    string f_osss_SR = Form("v%d2pc_fourierfit_OSminusSS_SR",  n);

    // Panel 1
    auto* gRaw_all      = readVnFile(f_all_LR, kBlack, 24);
    auto* gSubA_all     = makeSubtracted(f_all_LR, v1_all, refBin, kBlack, 20);

    // Panel 2 — OS-SS LR, 3 approaches
    auto* gRaw_osss_LR  = readVnFile(f_osss_LR, kRed+1, 24);
    auto* gSubA_osss_LR = makeSubtracted(f_osss_LR, v1_all,  refBin, kOrange+7, 20);
    auto* gSubC_osss_LR = makeSubtracted(f_osss_LR, v1_osss, refBin, kBlue+1,   22);
    auto* _gOS_LR       = makeSubtracted(f_os_LR,   v1_os,   refBin, kGreen+2,  20);
    auto* _gSS_LR       = makeSubtracted(f_ss_LR,   v1_ss,   refBin, kGreen+2,  20);
    auto* gSubB_osss_LR = combineGraphs(_gOS_LR, _gSS_LR, 1., -1., kGreen+2, 21);

    // Panel 3 — OS-SS SR, 3 approaches
    auto* gRaw_osss_SR  = readVnFile(f_osss_SR, kRed+1, 25);
    auto* gSubA_osss_SR = makeSubtracted(f_osss_SR, v1_all,  refBin, kOrange+7, 20);
    auto* gSubC_osss_SR = makeSubtracted(f_osss_SR, v1_osss, refBin, kBlue+1,   22);
    auto* _gOS_SR       = makeSubtracted(f_os_SR,   v1_os,   refBin, kGreen+2,  20);
    auto* _gSS_SR       = makeSubtracted(f_ss_SR,   v1_ss,   refBin, kGreen+2,  20);
    auto* gSubB_osss_SR = combineGraphs(_gOS_SR, _gSS_SR, 1., -1., kGreen+2, 21);

    // ── Global y-axis for panels 2 & 3: all OS-SS graphs both LR and SR ──
    auto [gYmin, gYmax] = globalRange({
        gRaw_osss_LR, gSubA_osss_LR, gSubB_osss_LR, gSubC_osss_LR,
        gRaw_osss_SR, gSubA_osss_SR, gSubB_osss_SR, gSubC_osss_SR
    });
    double dy23 = gYmax - gYmin; if (dy23 < 1e-9) dy23 = 0.001;
    double yMin_23 = gYmin - 0.15 * dy23;
    double yMax_23 = gYmax + 0.45 * dy23;

    auto* c = new TCanvas(Form("cA_n%d_ref%d",n,refBin), "c2 closure approaches", 1425, 560);
    c->Divide(3,1);

    string subTxt = Form("Periph. sub. w.r.t. Bin %d (%s)", refBin, binLabel(refBin).c_str());
    TLatex lat; lat.SetNDC();

    // ── Panel 1: all LR — independent scale ──
    c->cd(1); setPad(0.148, 0.084, 0.181, 0.060);
    auto* fr1 = makeFrame("frA1", Form("c_{%d} (All pairs, LR)",n));
    autoScale(fr1, {gRaw_all, gSubA_all});
    fr1->Draw(); zeroLine()->Draw();
    gRaw_all->Draw("P same"); gSubA_all->Draw("P same");
    auto* leg1 = new TLegend(0.220, 0.680, 0.959, 0.861, "", "brNDC");
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextFont(42); leg1->SetTextSize(0.038);
    leg1->AddEntry(gRaw_all,  "Raw All pairs", "p");
    leg1->AddEntry(gSubA_all, "Sub-A: c_{1}(all) scaler", "p");
    leg1->Draw();
    lat.SetTextFont(62); lat.SetTextSize(0.044);
    lat.DrawLatex(0.252, 0.919, "PYTHIA 8.316, p+p 200 GeV");
    lat.SetTextFont(42); lat.SetTextSize(0.036);
    lat.DrawLatex(0.232, 0.261, "Long-Range: |#Delta#eta| > 1.3");
    lat.DrawLatex(0.232, 0.215, subTxt.c_str());

    // ── Panel 2: OS-SS LR — global scale ──
    c->cd(2); setPad(0.168, 0.061, 0.181, 0.079);
    auto* fr2 = makeFrame("frA2", Form("c_{%d} (OS#minusSS, LR)",n));
    fr2->SetMinimum(yMin_23);
    fr2->SetMaximum(yMax_23);
    fr2->Draw(); zeroLine()->Draw();
    gRaw_osss_LR->Draw("P same");
    if (gSubA_osss_LR) gSubA_osss_LR->Draw("P same");
    if (gSubB_osss_LR) gSubB_osss_LR->Draw("P same");
    if (gSubC_osss_LR) gSubC_osss_LR->Draw("P same");
    auto* leg2 = new TLegend(0.234, 0.696, 0.975, 0.921, "", "brNDC");
    leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextFont(42); leg2->SetTextSize(0.038);
    leg2->AddEntry(gRaw_osss_LR,  "Raw (OS#minusSS)", "p");
    leg2->AddEntry(gSubA_osss_LR, "Sub-A: c_{1}(all) scaler", "p");
    leg2->AddEntry(gSubB_osss_LR, "Sub-B: [c_{n}(OS)-c_{n}(SS)]_{sub}", "p");
    leg2->AddEntry(gSubC_osss_LR, "Sub-C: c_{1}(OS#minusSS) scaler", "p");
    leg2->Draw();
    lat.SetTextFont(42); lat.SetTextSize(0.036);
    lat.DrawLatex(0.22, 0.24, "Long-Range: |#Delta#eta| > 1.3");

    // ── Panel 3: OS-SS SR — same global scale as panel 2 ──
    c->cd(3); setPad(0.20, 0.03, 0.18, 0.08);
    auto* fr3 = makeFrame("frA3", Form("c_{%d} (OS#minusSS, SR)",n));
    fr3->SetMinimum(yMin_23);
    fr3->SetMaximum(yMax_23);
    fr3->Draw(); zeroLine()->Draw();
    gRaw_osss_SR->Draw("P same");
    if (gSubA_osss_SR) gSubA_osss_SR->Draw("P same");
    if (gSubB_osss_SR) gSubB_osss_SR->Draw("P same");
    if (gSubC_osss_SR) gSubC_osss_SR->Draw("P same");
    auto* leg3 = new TLegend(0.22, 0.56, 0.97, 0.88, "", "brNDC");
    leg3->SetBorderSize(0); leg3->SetFillStyle(0); leg3->SetTextFont(42); leg3->SetTextSize(0.038);
    leg3->AddEntry(gRaw_osss_SR,  "Raw (OS#minusSS)", "p");
    leg3->AddEntry(gSubA_osss_SR, "Sub-A: c_{1}(all) scaler", "p");
    leg3->AddEntry(gSubB_osss_SR, "Sub-B: [c_{n}(OS)-c_{n}(SS)]_{sub}", "p");
    leg3->AddEntry(gSubC_osss_SR, "Sub-C: c_{1}(OS#minusSS) scaler", "p");
    leg3->Draw();
    lat.DrawLatex(0.22, 0.24, "Short-Range: |#Delta#eta| < 1.3");

    c->Update();
    string base = Form("canvasA_n%d_ref%d", n, refBin);
    c->SaveAs((base+".pdf").c_str()); c->SaveAs((base+".png").c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// CANVAS B: cn(OS) | cn(SS) | cn(CI) | cn(CD)  — for LR or SR
// All 4 panels share global y-axis (all 8 graphs: 4 raw + 4 subtracted)
// ─────────────────────────────────────────────────────────────────────────────
void drawCanvasB(int n, int refBin, const string& region) {
    bool isLR = (region == "LR");

    string v1_os = Form("v12pc_fourierfit_os_%s", region.c_str());
    string v1_ss = Form("v12pc_fourierfit_ss_%s", region.c_str());

    string f_os = Form("v%d2pc_fourierfit_os_%s", n, region.c_str());
    string f_ss = Form("v%d2pc_fourierfit_ss_%s", n, region.c_str());
    string f_ci = Form("v%d2pc_fourierfit_CI_%s", n, region.c_str());
    string f_cd = Form("v%d2pc_fourierfit_CD_%s", n, region.c_str());

    auto* gRaw_OS  = readVnFile(f_os, kRed+1,     24);
    auto* gSub_OS  = makeSubtracted(f_os, v1_os, refBin, kRed+1,     20);
    auto* gRaw_SS  = readVnFile(f_ss, kBlue+1,    25);
    auto* gSub_SS  = makeSubtracted(f_ss, v1_ss, refBin, kBlue+1,    21);
    auto* gRaw_CI  = readVnFile(f_ci, kGreen+2,   27);
    auto* gSubB_CI = combineGraphs(gSub_OS, gSub_SS,  0.5,  0.5, kMagenta+1, 20);
    auto* gRaw_CD  = readVnFile(f_cd, kMagenta+1, 27);
    auto* gSubB_CD = combineGraphs(gSub_OS, gSub_SS,  0.5, -0.5, kMagenta+1, 20);

    // ── Global y-range: scan all 8 graphs including error bars ──
    auto [gYmin, gYmax] = globalRange({
        gRaw_OS, gSub_OS,
        gRaw_SS, gSub_SS,
        gRaw_CI, gSubB_CI,
        gRaw_CD, gSubB_CD
    });
    double dy = gYmax - gYmin; if (dy < 1e-9) dy = 0.001;
    double yMin = gYmin - 0.15 * dy;
    double yMax = gYmax + 0.45 * dy;

    const char* etaLabel = isLR ? "Long-Range: |#Delta#eta| > 1.3"
                                : "Short-Range: |#Delta#eta| < 1.3";
    string subTxt = Form("Periph. sub. w.r.t. Bin %d (%s)", refBin, binLabel(refBin).c_str());

    auto* c = new TCanvas(Form("cB_%s_n%d_ref%d", region.c_str(), n, refBin),
                          Form("Individual charges + CI/CD (%s)", region.c_str()),
                          1728, 560);
    c->Divide(4,1);

    TLatex lat; lat.SetNDC();

    double lms[4] = {0.178, 0.162, 0.164, 0.164};
    double rms[4] = {0.051, 0.069, 0.065, 0.065};

    struct PanelDef {
        TGraphErrors* raw;
        TGraphErrors* sub;
        const char*   ytitle;
        const char*   legRaw;
    };
    PanelDef panels[4] = {
        {gRaw_OS,  gSub_OS,  Form("c_{%d} (OS, %s)",  n, region.c_str()), "Raw OS"},
        {gRaw_SS,  gSub_SS,  Form("c_{%d} (SS, %s)",  n, region.c_str()), "Raw SS"},
        {gRaw_CI,  gSubB_CI, Form("c_{%d} (CI, %s)",  n, region.c_str()), "Raw CI = (OS+SS)/2"},
        {gRaw_CD,  gSubB_CD, Form("c_{%d} (CD, %s)",  n, region.c_str()), "Raw CD = (OS#minusSS)/2"},
    };

    for (int p=0; p<4; ++p) {
        c->cd(p+1); setPad(lms[p], rms[p], 0.181, 0.081);
        auto* fr = makeFrame(Form("frB_%s_%d", region.c_str(), p+1), panels[p].ytitle);

        // Same global y-range for all 4 panels
        fr->SetMinimum(yMin);
        fr->SetMaximum(yMax);

        fr->Draw(); zeroLine()->Draw();
        if (panels[p].raw) panels[p].raw->Draw("P same");
        if (panels[p].sub) panels[p].sub->Draw("P same");

        auto* leg = new TLegend(0.214, 0.762, 0.953, 0.921, "", "brNDC");
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.038);
        if (panels[p].raw) leg->AddEntry(panels[p].raw, panels[p].legRaw, "p");
        if (panels[p].sub) leg->AddEntry(panels[p].sub, "Subtracted (own c_{1})", "p");
        leg->Draw();

        lat.SetTextFont(42); lat.SetTextSize(0.038);
        lat.DrawLatex(0.22, 0.24, etaLabel);

        if (p==0) {
            lat.SetTextFont(62); lat.SetTextSize(0.044);
            lat.DrawLatex(0.307, 0.930, "PYTHIA 8.316, p+p 200 GeV");
            lat.SetTextFont(42); lat.SetTextSize(0.038);
            lat.DrawLatex(0.281, 0.222, subTxt.c_str());
        }
    }

    c->Update();
    string base = Form("canvasB_%s_n%d_ref%d", region.c_str(), n, refBin);
    c->SaveAs((base+".pdf").c_str()); c->SaveAs((base+".png").c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// CANVAS C: c1(all,LR) | c1(OS-SS,LR) [3 approaches] | c1(OS-SS,SR) [3 approaches]
// Each panel keeps its own independent scale
// ─────────────────────────────────────────────────────────────────────────────
void drawCanvasC(int refBin) {
    string v1_all  = "v12pc_fourierfit_all_LR";
    string v1_os   = "v12pc_fourierfit_os_LR";
    string v1_ss   = "v12pc_fourierfit_ss_LR";
    string v1_osss = "v12pc_fourierfit_OSminusSS_LR";

    string f_all_LR  = "v12pc_fourierfit_all_LR";
    string f_os_LR   = "v12pc_fourierfit_os_LR";
    string f_ss_LR   = "v12pc_fourierfit_ss_LR";
    string f_osss_LR = "v12pc_fourierfit_OSminusSS_LR";
    string f_os_SR   = "v12pc_fourierfit_os_SR";
    string f_ss_SR   = "v12pc_fourierfit_ss_SR";
    string f_osss_SR = "v12pc_fourierfit_OSminusSS_SR";

    // Panel 1
    auto* gRaw_all      = readVnFile(f_all_LR, kBlack, 24);
    auto* gSubA_all     = makeSubtracted(f_all_LR, v1_all, refBin, kBlack, 20);

    // Panel 2 — c1 OS-SS LR
    auto* gRaw_osss_LR  = readVnFile(f_osss_LR, kRed+1,    24);
    auto* gSubA_osss_LR = makeSubtracted(f_osss_LR, v1_all,  refBin, kOrange+7, 20);
    auto* gSubC_osss_LR = makeSubtracted(f_osss_LR, v1_osss, refBin, kBlue+1,   22);
    auto* _gOS_LR       = makeSubtracted(f_os_LR,   v1_os,   refBin, kGreen+2,  20);
    auto* _gSS_LR       = makeSubtracted(f_ss_LR,   v1_ss,   refBin, kGreen+2,  20);
    auto* gSubB_osss_LR = combineGraphs(_gOS_LR, _gSS_LR, 1., -1., kGreen+2, 21);

    // Panel 3 — c1 OS-SS SR
    auto* gRaw_osss_SR  = readVnFile(f_osss_SR, kRed+1,    25);
    auto* gSubA_osss_SR = makeSubtracted(f_osss_SR, v1_all,  refBin, kOrange+7, 20);
    auto* gSubC_osss_SR = makeSubtracted(f_osss_SR, v1_osss, refBin, kBlue+1,   22);
    auto* _gOS_SR       = makeSubtracted(f_os_SR,   v1_os,   refBin, kGreen+2,  20);
    auto* _gSS_SR       = makeSubtracted(f_ss_SR,   v1_ss,   refBin, kGreen+2,  20);
    auto* gSubB_osss_SR = combineGraphs(_gOS_SR, _gSS_SR, 1., -1., kGreen+2, 21);

    auto* c = new TCanvas(Form("cC_ref%d",refBin), "c1 closure approaches", 1425, 560);
    c->Divide(3,1);

    string subTxt = Form("Periph. sub. w.r.t. Bin %d (%s)", refBin, binLabel(refBin).c_str());
    TLatex lat; lat.SetNDC();

    // ── Panel 1: c1 all LR (tautology) ──
    c->cd(1); setPad(0.169, 0.062, 0.161, 0.060);
    auto* fr1 = makeFrame("frC1", "c_{1} (All pairs, LR)");
    autoScale(fr1, {gRaw_all, gSubA_all});
    fr1->Draw(); zeroLine()->Draw();
    gRaw_all->Draw("P same");
    if (gSubA_all) gSubA_all->Draw("P same");
    auto* leg1 = new TLegend(0.214, 0.760, 0.953, 0.881, "", "brNDC");
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextFont(42); leg1->SetTextSize(0.036);
    leg1->AddEntry(gRaw_all,  "Raw c_{1}(all)", "p");
    leg1->AddEntry(gSubA_all, "Sub-A: c_{1}(all)/c_{1}(all) [tautology #equiv 0]", "p");
    leg1->Draw();
    lat.SetTextFont(62); lat.SetTextSize(0.044);
    lat.DrawLatex(0.292, 0.928, "PYTHIA 8.316, p+p 200 GeV");
    lat.SetTextFont(42); lat.SetTextSize(0.036);
    lat.DrawLatex(0.206, 0.255, "Long-Range: |#Delta#eta| > 1.3");
    lat.DrawLatex(0.206, 0.205, subTxt.c_str());

    // ── Panel 2: c1 OS-SS LR — legend at BOTTOM ──
    c->cd(2); setPad(0.150, 0.079, 0.177, 0.083);
    auto* fr2 = makeFrame("frC2", "c_{1} (OS#minusSS, LR)");
    autoScale(fr2, {gRaw_osss_LR, gSubA_osss_LR, gSubB_osss_LR, gSubC_osss_LR});
    fr2->Draw(); zeroLine()->Draw();
    gRaw_osss_LR->Draw("P same");
    if (gSubA_osss_LR) gSubA_osss_LR->Draw("P same");
    if (gSubB_osss_LR) gSubB_osss_LR->Draw("P same");
    if (gSubC_osss_LR) gSubC_osss_LR->Draw("P same");
    auto* leg2 = new TLegend(0.240, 0.295, 0.981, 0.507, "", "brNDC");
    leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextFont(42); leg2->SetTextSize(0.036);
    leg2->AddEntry(gRaw_osss_LR,  "Raw c_{1}(OS#minusSS)", "p");
    leg2->AddEntry(gSubA_osss_LR, "Sub-A: c_{1}(all) scaler", "p");
    leg2->AddEntry(gSubB_osss_LR, "Sub-B: [c_{1}(OS)-c_{1}(SS)]_{sub}", "p");
    leg2->AddEntry(gSubC_osss_LR, "Sub-C: c_{1}(OS#minusSS) scaler [#equiv 0]", "p");
    leg2->Draw();
    lat.SetTextFont(42); lat.SetTextSize(0.036);
    lat.DrawLatex(0.22, 0.24, "Long-Range: |#Delta#eta| > 1.3");

    // ── Panel 3: c1 OS-SS SR — legend at BOTTOM ──
    c->cd(3); setPad(0.20, 0.03, 0.18, 0.08);
    auto* fr3 = makeFrame("frC3", "c_{1} (OS#minusSS, SR)");
    autoScale(fr3, {gRaw_osss_SR, gSubA_osss_SR, gSubB_osss_SR, gSubC_osss_SR});
    fr3->Draw(); zeroLine()->Draw();
    gRaw_osss_SR->Draw("P same");
    if (gSubA_osss_SR) gSubA_osss_SR->Draw("P same");
    if (gSubB_osss_SR) gSubB_osss_SR->Draw("P same");
    if (gSubC_osss_SR) gSubC_osss_SR->Draw("P same");
    auto* leg3 = new TLegend(0.240, 0.295, 0.981, 0.507, "", "brNDC");
    leg3->SetBorderSize(0); leg3->SetFillStyle(0); leg3->SetTextFont(42); leg3->SetTextSize(0.036);
    leg3->AddEntry(gRaw_osss_SR,  "Raw c_{1}(OS#minusSS)", "p");
    leg3->AddEntry(gSubA_osss_SR, "Sub-A: c_{1}(all) scaler", "p");
    leg3->AddEntry(gSubB_osss_SR, "Sub-B: [c_{1}(OS)-c_{1}(SS)]_{sub}", "p");
    leg3->AddEntry(gSubC_osss_SR, "Sub-C: c_{1}(OS#minusSS) scaler [#equiv 0]", "p");
    leg3->Draw();
    lat.DrawLatex(0.22, 0.24, "Short-Range: |#Delta#eta| < 1.3");

    c->Update();
    string base = Form("canvasC_c1_ref%d", refBin);
    c->SaveAs((base+".pdf").c_str()); c->SaveAs((base+".png").c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────────────────────────────────────
void plot_closure_LR_vs_SR_pp200(int n=2, int refBin=8) {
    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(3);

    drawCanvasA(n, refBin);         // panels 2&3 share global y-axis
    drawCanvasB(n, refBin, "LR");   // all 4 panels share global y-axis
    drawCanvasB(n, refBin, "SR");   // all 4 panels share global y-axis
    drawCanvasC(refBin);            // c1 sanity check, independent scales
}

