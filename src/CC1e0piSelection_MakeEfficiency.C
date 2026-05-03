#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Tree.h"

// helpers — include order matters!
#include "../include/Config.h"
#include "../include/CC1e0piSelection_Vars_Base.h"
#include "../include/CC1e0piSelection_Vars_Reco.h"
#include "../include/CC1e0piSelection_Vars_Derived.h"
#include "../include/CC1e0piSelection_Vars_NuGraph.h"
#include "../include/CC1e0piSelection_Cuts.h"
#include "../include/CC1e0piSelection_TruthCuts.h"
#include "../include/CC1e0piSelection_Efficiency.h"
#include "../include/CC1e0piSelection_Ntuple.h"

// ROOT
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TPad.h"
#include "TDirectory.h"

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

using namespace ana;

// -----------------------------------------------------------------------------
// Style
// -----------------------------------------------------------------------------
void SetPrettyStyle()
{
    gStyle->SetOptStat(0);
    const int font = 132;
    gStyle->SetTextFont(font);
    gStyle->SetLabelFont(font, "XYZ");
    gStyle->SetTitleFont(font, "XYZ");
    gStyle->SetLegendFont(font);
    gStyle->SetLabelSize(0.048, "XYZ");
    gStyle->SetTitleSize(0.062, "XYZ");
    gStyle->SetTitleOffset(0.95, "X");
    gStyle->SetTitleOffset(0.85, "Y");
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetTickLength(0.02, "X");
    gStyle->SetTickLength(0.02, "Y");
    gStyle->SetNdivisions(506, "XY");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetEndErrorSize(0);
    gStyle->SetLineScalePS(2.0);
}

// -----------------------------------------------------------------------------
// Color palette
// -----------------------------------------------------------------------------
std::vector<int> GetPrettyColors(bool triggerMode)
{
    std::vector<std::string> hex_trigger = {
        "#111111", "#8c1d18", "#f26b6c", "#f39c12",
        "#3f8f3f", "#6ddf5f", "#61d6d6", "#4257f5", "#c13fd4"
    };
    std::vector<std::string> hex_notrigger = {
        "#111111", "#8c1d18", "#f26b6c", "#f39c12",
        "#3f8f3f", "#61d6d6", "#4257f5", "#c13fd4"
    };
    const auto& hex = triggerMode ? hex_trigger : hex_notrigger;
    std::vector<int> cols;
    for (const auto& h : hex) cols.push_back(TColor::GetColor(h.c_str()));
    return cols;
}

// -----------------------------------------------------------------------------
// Efficiency summary plot — no shaded bands, with truth backdrop
// -----------------------------------------------------------------------------
void DrawEfficiencySummary(const std::string& canvasName,
                           const std::string& legendTitle,
                           TH1* hTrue,
                           const std::vector<Spectrum*>& spectra,
                           const std::vector<SelDef>& steps,
                           bool triggerMode,
                           TFile& fout)
{
    std::vector<int> cols = GetPrettyColors(triggerMode);

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 1350);
    c->cd();
    c->SetTicks(1, 1);

    // truth energy shape backdrop
    TH1* hTrueScaled = (TH1*)hTrue->Clone((canvasName + "_backdrop").c_str());
    hTrueScaled->SetDirectory(nullptr);
    hTrueScaled->Scale(1. / hTrueScaled->GetMaximum());
    hTrueScaled->SetFillColorAlpha(TColor::GetColor("#bdbdbd"), 0.28);
    hTrueScaled->SetLineColor(0);
    hTrueScaled->SetTitle("");
    hTrueScaled->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    hTrueScaled->GetYaxis()->SetTitle("Selection efficiency");
    hTrueScaled->GetYaxis()->SetRangeUser(0.0, 1.35);
    hTrueScaled->GetXaxis()->SetRangeUser(0.0, 4.0);
    hTrueScaled->Draw("HIST");

    TLine* unity = new TLine(0.0, 1.0, 4.0, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->SetLineColor(TColor::GetColor("#555555"));
    unity->Draw("SAME");

    TLegend* leg = new TLegend(0.15, 0.80, 0.85, 0.93);
    leg->SetNColumns(3);
    leg->SetTextAlign(22);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);

    TDirectory* subdir = fout.mkdir(("details_" + canvasName).c_str());

    for (unsigned int iSel = 0; iSel < steps.size(); ++iSel) {
        TH1* hSel = spectra[iSel]->ToTH1(spectra[iSel]->POT());
        hSel->SetDirectory(nullptr);

        const double nPass = hSel->GetEntries();
        const double nTot  = hTrue->GetEntries();
        const double frac  = (nTot > 0.0) ? nPass / nTot : 0.0;

        TEfficiency* eff = new TEfficiency(*hSel, *hTrue);
        TGraphAsymmErrors* g = eff->CreateGraph();
        g->SetName((canvasName + "_g_" + steps[iSel].suffix + "_" + std::to_string(iSel)).c_str());
        g->SetLineColor(cols[iSel]);
        g->SetMarkerColor(cols[iSel]);
        g->SetLineWidth(3);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.0);
        g->Draw("PZ SAME");

        subdir->cd();
        g->Write();
        fout.cd();

        TString legLabel = steps[iSel].label + Form(" (%.0f%%)", 100.0 * frac);
        leg->AddEntry(g, legLabel, "lf");

        std::cout << canvasName << " : "
                  << steps[iSel].label << "  "
                  << nPass << " / " << nTot
                  << "  ->  " << frac << std::endl;
    }

    leg->Draw();

    subdir->cd();
    hTrueScaled->Write();
    fout.cd();
    c->Write();

    delete c;
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
void CC1e0piSelection_MakeEfficiency()
{
    //const std::string listFile = "/exp/icarus/app/users/csaia/bnb_nominalflux_caf_small.list"; //BNB MC file list (500 files ~ 10k events) in v10_06_00_06p03 used for debug-only (too low stats for reliable results --> huge error bars)
    const std::string listFile = "/exp/icarus/app/users/csaia/bnb_nominalflux_caf.list"; //BNB MC file list (All files!) in v10_06_00_06p03 
    std::vector<std::string> fileList;
    std::ifstream fin(listFile);
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty() && line[0] != '#')
            fileList.push_back(line);
    }
    if (fileList.empty()) {
        std::cerr << "Error: no files loaded from " << listFile << std::endl;
        return;
    }

    SpectrumLoader NuLoader(fileList);

    // -------------------------------------------------------------------------
    // Efficiency spectra
    // -------------------------------------------------------------------------
    std::vector<double> TrueEnergyBinning = {
        0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 4.0
    };

    Spectrum* sTrueNeutrinoEnergy = new Spectrum(
        "E_{#nu} [GeV]", Binning::Custom(TrueEnergyBinning),
        NuLoader, kCC1e0p1Signal_TrueNeutrinoEnergy, kNoSpillCut);

    const unsigned int nTrig = SelectionSteps.size();
    std::vector<Spectrum*> sTrig(nTrig, nullptr);
    for (unsigned int iSel = 0; iSel < nTrig; ++iSel) {
        sTrig[iSel] = new Spectrum(
            "E_{#nu} [GeV]", Binning::Custom(TrueEnergyBinning),
            NuLoader,
            kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps[iSel].cut),
            kNoSpillCut);
    }

    const unsigned int nNoTrig = SelectionSteps_NoTrigger.size();
    std::vector<Spectrum*> sNoTrig(nNoTrig, nullptr);
    for (unsigned int iSel = 0; iSel < nNoTrig; ++iSel) {
        sNoTrig[iSel] = new Spectrum(
            "E_{#nu} [GeV]", Binning::Custom(TrueEnergyBinning),
            NuLoader,
            kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps_NoTrigger[iSel].cut),
            kNoSpillCut);
    }

    // -------------------------------------------------------------------------
    // Debug ntuple - loose pre-selection, all raw variables
    // Flat branches via CAFAna Tree, per-PFP vectors via separate Tree
    // -------------------------------------------------------------------------
    std::vector<std::string> branchNames;
    std::vector<Var> branchVars;
    for (const auto& b : NtupleBranches) {
        branchNames.push_back(b.name);
        branchVars.push_back(b.var);
    }

    std::vector<std::string> vecBranchNames;
    std::vector<MultiVar> vecBranchVars;
    for (const auto& b : NtupleVecBranches) {
        vecBranchNames.push_back(b.name);
        vecBranchVars.push_back(b.var);
    }

    // flat branches
    Tree debugTree("ntuple", branchNames, NuLoader, branchVars,
                   kNoSpillCut, kNtuplePreSelection, kNoShift, true, true);

    // per-PFP vector branches
    Tree debugTreeVec("ntuple_pfp", vecBranchNames, NuLoader, vecBranchVars,
                      kNoSpillCut, kNtuplePreSelection, kNoShift, true, true);

    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    NuLoader.Go();

    gROOT->SetBatch(kTRUE);
    SetPrettyStyle();

    TFile FOut("../output/root/CC1e0piSelection_EfficiencyTest.root", "RECREATE");

    // hTrueSignal
    const double TargetPOT = sTrueNeutrinoEnergy->POT();
    TH1* hTrue = sTrueNeutrinoEnergy->ToTH1(TargetPOT);
    hTrue->SetDirectory(nullptr);
    hTrue->SetName("hTrueSignal");
    hTrue->SetTitle(";E_{#nu} [GeV];Signal events");
    hTrue->SetMarkerStyle(20);
    hTrue->SetMarkerSize(1.0);
    hTrue->SetMarkerColor(TColor::GetColor("#111111"));
    hTrue->SetLineColor(TColor::GetColor("#888888"));
    hTrue->SetLineWidth(2);
    hTrue->GetXaxis()->SetTitleFont(132);
    hTrue->GetXaxis()->SetLabelFont(132);
    hTrue->GetYaxis()->SetTitleFont(132);
    hTrue->GetYaxis()->SetLabelFont(132);
    hTrue->GetXaxis()->SetTitleSize(0.062);
    hTrue->GetYaxis()->SetTitleSize(0.062);
    hTrue->GetXaxis()->SetLabelSize(0.048);
    hTrue->GetYaxis()->SetLabelSize(0.048);
    FOut.cd();
    hTrue->Write();

    // efficiency plots
    DrawEfficiencySummary("efficiencyCC1e0pi_trigger",
                          "BNB CV (trigger)",
                          hTrue, sTrig, SelectionSteps, true, FOut);

    DrawEfficiencySummary("efficiencyCC1e0pi_notrigger",
                          "BNB CV (no trigger)",
                          hTrue, sNoTrig, SelectionSteps_NoTrigger, false, FOut);

    // debug ntuple
    TDirectory* ntupleDir = FOut.mkdir("ntuple");
    debugTree.SaveTo(ntupleDir);

    TDirectory* ntuplePFPDir = FOut.mkdir("ntuple_pfp");
    debugTreeVec.SaveTo(ntuplePFPDir);

    FOut.Close();

    return;
}
