#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>
#include <map>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TColor.h"
#include "TString.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TGaxis.h"

using namespace std;

// --- Structures ---
struct SeriesRaw {
    string fname;
    string label;
    int color;
    int marker;
    double xshift;
    double scale;
};

struct SubTask {
    string vn_file;      // Signal file
    string v1_file;      // Scaler file (v1)
    string label;
    int color;
    int marker;
    double xshift;
};

// --- Helper: Get Reference Value (Bin 8 = 0-6 Mult) ---
pair<double, double> getRefValue(string fname) {
    ifstream in(fname.c_str());
    pair<double, double> result = {0.0, 0.0};

    if (!in.is_open()) return result; // Silent fail

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        double val;
        vector<double> row;
        while (ss >> val) row.push_back(val);

        if (row.size() < 7) continue;

        // PYTHIA: Bin 8 is the lowest multiplicity (0-6)
        if (std::abs(row[6] - 8.0) < 0.01) {
            result.first = row[4];  // y value
            result.second = row[5]; // y error
            break;
        }
    }
    in.close();
    return result;
}

// --- 1. Read RAW Graph (Unsubtracted) ---
TGraphErrors* readRawGraph(const SeriesRaw &s) {
    ifstream in(s.fname.c_str());
    if (!in.is_open()) {
        // Silent skip as requested
        return nullptr;
    }

    vector<double> vx, vy, vey;
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        double val;
        vector<double> row;
        while (ss >> val) row.push_back(val);

        if (row.size() < 7) continue;

        // Plot all bins 0-8 directly
        vx.push_back(row[6] + s.xshift); 
        vy.push_back(row[4] * s.scale);            
        vey.push_back(row[5] * s.scale);           
    }
    in.close();

    if (vx.empty()) return nullptr;

    TGraphErrors *g = new TGraphErrors(vx.size());
    for (size_t i = 0; i < vx.size(); i++) {
        g->SetPoint(i, vx[i], vy[i]);
        g->SetPointError(i, 0.0, vey[i]);
    }
    g->SetMarkerStyle(s.marker);
    g->SetMarkerSize(1.6);
    g->SetMarkerColor(s.color);
    g->SetLineColor(s.color);
    g->SetLineWidth(2);
    return g;
}

// --- 2. Create SUBTRACTED Graph ---
TGraphErrors* makeSubtractedGraph(const SubTask &task) {
    // 1. Get Reference Values
    pair<double, double> vn_ref = getRefValue(task.vn_file); 
    pair<double, double> v1_ref = getRefValue(task.v1_file); 

    double C = v1_ref.first;
    double sigma_C = v1_ref.second;
    double D = vn_ref.first;
    double sigma_D = vn_ref.second;

    // Safety check for div by zero
    if (abs(C) < 1e-9) return nullptr;

    // 2. Open Files
    ifstream in_vn(task.vn_file.c_str());
    ifstream in_v1(task.v1_file.c_str());
    if (!in_vn.is_open() || !in_v1.is_open()) return nullptr;

    // Map v1 for lookup
    std::map<int, pair<double, double>> v1_map;
    string line;
    while(getline(in_v1, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        vector<double> row; 
        while(ss >> row.emplace_back()); // simple read
        if(row.size() >= 7) v1_map[(int)row[6]] = {row[4], row[5]};
    }
    in_v1.close();

    vector<double> vx, vy, vey;
    while (getline(in_vn, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        double val;
        vector<double> row;
        while (ss >> val) row.push_back(val);

        if (row.size() < 7) continue;

        int binIdx = (int)row[6];
        double A = row[4];       // v_n raw
        double sigma_A = row[5]; // v_n raw error

        if (v1_map.find(binIdx) == v1_map.end()) continue;
        double B = v1_map[binIdx].first;        // v_1 raw
        double sigma_B = v1_map[binIdx].second; // v_1 raw error

        // Subtraction: Y = A - (B * D) / C
        double term = (B * D) / C;
        double Y = A - term;

        // Error Propagation
        double term_err_sq = pow((D/C)*sigma_B, 2) + 
                             pow((B/C)*sigma_D, 2) + 
                             pow((B*D)/(C*C)*sigma_C, 2);
        double sigma_Y = sqrt(pow(sigma_A, 2) + term_err_sq);

        vx.push_back((double)binIdx + task.xshift);
        vy.push_back(Y);
        vey.push_back(sigma_Y);
    }
    in_vn.close();

    if (vx.empty()) return nullptr;

    TGraphErrors *g = new TGraphErrors(vx.size());
    for (size_t i = 0; i < vx.size(); i++) {
        g->SetPoint(i, vx[i], vy[i]);
        g->SetPointError(i, 0.0, vey[i]);
    }
    g->SetMarkerStyle(task.marker);
    g->SetMarkerSize(1.6);
    g->SetMarkerColor(task.color);
    g->SetLineColor(task.color);
    g->SetLineWidth(2);
    return g;
}

// --- 3. Template Fit Reader ---
TGraphErrors* readTemplateGraph(string fname, int color, int marker, double xshift) {
    ifstream in(fname.c_str());
    if (!in.is_open()) return nullptr;

    vector<double> vx, vy, vey;
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        vector<double> row;
        double val;
        while (ss >> val) row.push_back(val);
        if (row.size() < 7) continue;

        vx.push_back(row[6] + xshift);
        vy.push_back(row[4]);
        vey.push_back(row[5]);
    }
    in.close();

    if (vx.empty()) return nullptr;
    TGraphErrors *g = new TGraphErrors(vx.size());
    for (size_t i = 0; i < vx.size(); i++) {
        g->SetPoint(i, vx[i], vy[i]);
        g->SetPointError(i, 0.0, vey[i]);
    }
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.6);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
    return g;
}

// --- MAIN ---
void plot_vn_both_sub_unsub_pp200(int n = 2, double scaleFactor = 0.5, bool twoColumnLegend = false) 
{
    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(6); // Force decimal

    TCanvas *c = new TCanvas("c_vn_pythia", "vn combined pythia", 29, 94, 820, 609);
    c->SetLeftMargin(0.16);
    c->SetRightMargin(0.041);
    c->SetTopMargin(0.031);
    c->SetBottomMargin(0.228); 
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);

    string ytitle = Form("c_{%d}=v_{%d}^{2}",n, n);
    string scaleStr = Form("(#times%g)", scaleFactor);

    // --- Define Series ---
    vector<SeriesRaw> rawTasks;
    vector<SubTask> subTasks;
    
    // We will build a unified list of graphs to draw and add to legend in specific order
    // Order requested:
    // 1. Raw Full-evnt.
    // 2. Raw Sub-Event
    // 3. Raw Fourier
    // 4. Sub Full-evnt.
    // 5. Sub Sub-Event
    // 6. Sub Fourier
    // 7. Template Fit

    if (n == 1) {
        // v1 special case (Raw v1 is plotted, Sub v1 is zero/check)
        rawTasks = {
            { "v12_fullEv_deta1p3_nofit", "Raw Full-evnt. v_{1}^{2}{2} #Delta#eta>1.3", kGray+2, 24, -0.25, 1.0 },
            { "v12_subEv_nofit",          "Raw Sub-evnt. V_{1,1} (#eta-gap=1.3)",        kCyan-2, 26, -0.20, 1.0 },
            { "v12pc_fourierfit",         "Raw Fourier fit to C(#Delta#phi)",           kBlue-7, 25, -0.15, 1.0 }
        };
        subTasks = {
            { "v12_fullEv_deta1p3_nofit", "v12_fullEv_deta1p3_nofit", "c1-sub Full-evnt. v_{1}^{2}{2} #Delta#eta>1.3",           kBlack,   20, -0.05 },
            { "v12_subEv_nofit",          "v12_subEv_nofit",          "c1-sub Sub-evnt. V_{1,1} (#eta-gap=1.3)",                  kBlue+1,  22, 0.00 },
            { "v12pc_fourierfit",         "v12pc_fourierfit",         "c1-sub Fourier fit to C(#Delta#phi)",                      kAzure-1, 33, 0.05 }
        };
        // Template usually doesn't exist for v1, but we keep the logic generic
    } 
    else { // n=2, 3
        string vn = Form("v%d", n); // v2, v3
        
        rawTasks = {
            { vn+"2_fullEv_deta1p3_nofit", "Raw Full-evnt. v_{"+to_string(n)+"}^{2}{2} #Delta#eta>1.3",                     kGray+2, 24, -0.25, 1.0 },
            { vn+"2_subEv_nofit",          "Raw Sub-evnt. V_{"+to_string(n)+","+to_string(n)+"} (#eta-gap=1.3)",            kCyan-2, 26, -0.20, 1.0 },
            { vn+"2pc_fourierfit",         "Raw Fourier fit to C(#Delta#phi), V_{"+to_string(n)+"#Delta}(2PC; #Delta#eta>1.3)", kBlue-7, 25, -0.15, 1.0 }
        };
        
        subTasks = {
            { vn+"2_fullEv_deta1p3_nofit", "v12_fullEv_deta1p3_nofit", "c1-sub Full-evnt. v_{"+to_string(n)+"}^{2}{2} #Delta#eta>1.3",           kBlack,   20, -0.05 },
            { vn+"2_subEv_nofit",          "v12_subEv_nofit",          "c1-sub Sub-evnt. V_{"+to_string(n)+","+to_string(n)+"} (#eta-gap=1.3)",  kGreen-2,  22, 0.00 },
            { vn+"2pc_fourierfit",         "v12pc_fourierfit",         "c1-sub Fourier fit to C(#Delta#phi), V_{"+to_string(n)+"#Delta}(2PC; #Delta#eta>1.3)", kAzure-1, 21, 0.05 }
        };
    }

    // --- Generate Graphs & Store in Order ---
    vector<TGraphErrors*> G_final;
    vector<string> labels_final;

    // 1. Raw Graphs
    for (auto &t : rawTasks) {
        auto g = readRawGraph(t);
        if (g) { G_final.push_back(g); labels_final.push_back(t.label); }
    }

    // 2. Subtracted Graphs
    for (auto &t : subTasks) {
        auto g = makeSubtractedGraph(t);
        if (g) { G_final.push_back(g); labels_final.push_back(t.label); }
    }

    // 3. Template Fit (Benchmark)
    if (n > 1) {
        string tempFile = Form("v%d2pc_templatefit", n);
        string tempLabel = Form("Template sub Fourier fit to C(#Delta#phi), V_{%d#Delta}(2PC; #Delta#eta>1.3)", n);
        auto g_temp = readTemplateGraph(tempFile, kOrange+7, 34, 0.15); // Cross marker
        if (g_temp) { G_final.push_back(g_temp); labels_final.push_back(tempLabel); }
    }

    if (G_final.empty()) {
        cerr << "No graphs produced. Check file names." << endl;
        return;
    }

    // --- Auto Range ---
    double ymin = 1e9, ymax = -1e9;
    for (auto g : G_final) {
        for (int i=0; i<g->GetN(); i++) {
            double x,y; g->GetPoint(i,x,y);
            if (y < -100) continue; // Skip outliers
            double ey = g->GetErrorY(i);
            if (y-ey < ymin) ymin = y-ey;
            if (y+ey > ymax) ymax = y+ey;
        }
    }
    if (ymin > ymax) { ymin = -0.001; ymax = 0.001; }
    double dy = ymax - ymin;
    if (dy == 0) dy = 1.0;
    ymin -= 0.1 * dy;
    ymax += 0.5 * dy;

    // --- Frame with Pythia Labels ---
    // Bins 0 to 8
    TH1F *frame = new TH1F("frame", "", 9, -0.5, 8.5);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->SetStats(0);
    frame->GetYaxis()->SetNoExponent(kTRUE);

    const char *labels[9] = {
        "49-60", // Bin 0
        "43-48", // Bin 1
        "37-42", // Bin 2
        "31-36", // Bin 3
        "25-30", // Bin 4
        "19-24", // Bin 5
        "13-18", // Bin 6
        "7-12",  // Bin 7
        "0-6"    // Bin 8 (Ref)
    };
    for (int i=1; i<=9; i++) frame->GetXaxis()->SetBinLabel(i, labels[i-1]);

    frame->GetXaxis()->SetTitle("Activity Class (Tracks |#eta|<1.5)");
    frame->GetXaxis()->SetBit(TAxis::kLabelsVert);
    frame->GetXaxis()->CenterTitle(true);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetTitleFont(42);

    frame->GetYaxis()->SetTitle(ytitle.c_str());
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleOffset(1.35);
    frame->GetYaxis()->SetTitleFont(42);

    frame->Draw();

    // --- Text Box ---
    TPaveText *pt = new TPaveText(0.19, 0.74, 0.55, 0.89, "brNDC");
    pt->SetBorderSize(0); pt->SetFillColor(0); pt->SetFillStyle(0); pt->SetTextAlign(12);
    TText *t1 = pt->AddText("PYTHIA 8.316, p+p 200 GeV"); t1->SetTextFont(62); t1->SetTextSize(0.045);
    TText *t2 = pt->AddText("Detroit Tune, Non-diffractive"); t2->SetTextFont(42); t2->SetTextSize(0.035);
    TText *t3 = pt->AddText("0.2 < p^{trig,asco}_{T} < 2.0 GeV/c, |#eta|<1.5"); t3->SetTextFont(42); t3->SetTextSize(0.035);
    pt->Draw();

    // --- Legend ---
    TLegend *leg;
    if (twoColumnLegend) {
        leg = new TLegend(0.20, 0.52, 0.75, 0.75, NULL, "brNDC");
        leg->SetNColumns(2);
        leg->SetTextSize(0.032);
        leg->SetMargin(0.15);
    } else {
        leg = new TLegend(0.5, 0.52, 0.70, 0.79, NULL, "brNDC");
        leg->SetNColumns(1);
        leg->SetTextSize(0.035);
    }
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

    for (size_t i = 0; i < G_final.size(); i++) {
        G_final[i]->Draw("P same");
        leg->AddEntry(G_final[i], labels_final[i].c_str(), "lep");
    }

    leg->Draw();
    c->Update();
}
