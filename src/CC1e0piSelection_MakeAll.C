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
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TBox.h"
#include "TDirectory.h"

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

using namespace ana;

// =============================================================================
// Style
// =============================================================================
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

// =============================================================================
// Color palettes
// =============================================================================
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

std::vector<int> GetViridisColors(int n)
{
    gStyle->SetPalette(kViridis);
    const int nPalette = TColor::GetNumberOfColors();
    std::vector<int> cols;
    for (int i = 0; i < n; ++i) {
        double frac = (n > 1) ? (double)i / (n - 1) : 0.5;
        cols.push_back(TColor::GetColorPalette((int)(frac * (nPalette - 1))));
    }
    return cols;
}

// =============================================================================
// Efficiency summary
// =============================================================================
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

// =============================================================================
// Purity summary
// =============================================================================
void DrawPuritySummary(const std::string& canvasName,
                       const PlotDef& plotDef,
                       TH1* hSignal,
                       const std::vector<TH1*>& hAllVec,
                       const std::vector<TH1*>& hSelectedVec,
                       const std::vector<SelDef>& steps,
                       bool triggerMode,
                       TDirectory* parentDir)
{
    std::vector<int> cols = GetPrettyColors(triggerMode);

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 1350);
    c->cd();
    c->SetTicks(1, 1);

    TH1* hSigScaled = (TH1*)hSignal->Clone((canvasName + "_backdrop").c_str());
    hSigScaled->SetDirectory(nullptr);
    hSigScaled->Scale(1. / hSigScaled->GetMaximum());
    hSigScaled->SetFillColorAlpha(TColor::GetColor("#bdbdbd"), 0.28);
    hSigScaled->SetLineColor(0);
    hSigScaled->SetTitle("");
    hSigScaled->GetXaxis()->SetTitle(plotDef.label.c_str());
    hSigScaled->GetXaxis()->SetTitleFont(132);
    hSigScaled->GetXaxis()->SetLabelFont(132);
    hSigScaled->GetYaxis()->SetTitle("Selection purity");
    hSigScaled->GetYaxis()->SetTitleFont(132);
    hSigScaled->GetYaxis()->SetLabelFont(132);
    hSigScaled->GetYaxis()->SetRangeUser(0.0, 1.35);
    hSigScaled->Draw("HIST");

    TLine* unity = new TLine(hSigScaled->GetXaxis()->GetXmin(), 1.0,
                             hSigScaled->GetXaxis()->GetXmax(), 1.0);
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

    TDirectory* subdir = parentDir->mkdir(("details_" + canvasName).c_str());

    for (unsigned int iSel = 0; iSel < steps.size(); ++iSel) {
        TH1* hAll      = hAllVec[iSel];
        TH1* hSelected = hSelectedVec[iSel];

        const double nPass = hSelected->GetEntries();
        const double nTot  = hAll->GetEntries();
        const double frac  = (nTot > 0.0) ? nPass / nTot : 0.0;

        std::string effName = canvasName + "_eff_" + steps[iSel].suffix + "_" + std::to_string(iSel);
        TEfficiency* eff = new TEfficiency(effName.c_str(), "",
                                           plotDef.bins.NBins(), &plotDef.bins.Edges()[0]);
        for (int i = 1; i <= plotDef.bins.NBins(); ++i) {
            if (hAll->GetBinContent(i) == 0) continue;
            eff->SetTotalEvents(i, hAll->GetBinContent(i));
            eff->SetPassedEvents(i, hSelected->GetBinContent(i));
        }

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
        parentDir->cd();

        TString legLabel = steps[iSel].label + Form(" (%.0f%%)", 100.0 * frac);
        leg->AddEntry(g, legLabel, "lf");

        std::cout << canvasName << " : "
                  << steps[iSel].label << "  "
                  << nPass << " / " << nTot
                  << "  ->  " << frac << std::endl;
    }

    leg->Draw();
    subdir->cd();
    hSigScaled->Write();
    parentDir->cd();
    c->Write();
    delete c;
}

// =============================================================================
// Selection stacked plot
// =============================================================================
void DrawSelectionPlot(const PlotDef& plotDef,
                       const std::vector<Spectrum*>& spectra,
                       const std::vector<int>& cols,
                       TDirectory* parentDir)
{
    const std::string& suffix = plotDef.suffix;

    TCanvas* c = new TCanvas(suffix.c_str(), suffix.c_str(), 1800, 1350);
    THStack* hs = new THStack(suffix.c_str(), plotDef.label.c_str());
    TLegend* leg = new TLegend(0.65, 0.5, 0.85, 0.85, "BNB CV");

    const double TargetPOT = spectra[0]->POT();
    TH1* hAll = spectra[0]->ToTH1(TargetPOT);

    double yMax = 0.;
    for (int i = 1; i <= hAll->GetNbinsX(); ++i) {
        double y = hAll->GetBinContent(i) + hAll->GetBinError(i);
        if (y > yMax) yMax = y;
    }

    for (unsigned int jSel = 1; jSel < spectra.size(); ++jSel) {
        TH1* h = spectra[jSel]->ToTH1(spectra[jSel]->POT());
        h->SetFillColor(cols[jSel - 1]);
        h->SetLineWidth(0);
        h->SetLineColor(kBlack);

        std::string hName = suffix + "_" + InteractionTypes[jSel].suffix;
        h->SetName(hName.c_str());
        parentDir->cd();
        h->Write();

        leg->AddEntry(h, InteractionTypes[jSel].label.c_str(), "f");
        hs->Add(h);
    }

    hs->SetMaximum(yMax + 0.1 * yMax);
    hs->Draw("HIST");

    std::string title = std::string(";") + plotDef.label + std::string(";") +
                        Form("Slices / %.1e POT", TargetPOT);
    hs->SetTitle(title.c_str());
    gPad->Modified();
    gPad->Update();

    // error band on total
    gStyle->SetHatchesLineWidth(3);
    gStyle->SetHatchesSpacing(1.5);
    for (int i = 1; i <= hAll->GetNbinsX(); ++i) {
        double xlow = hAll->GetBinLowEdge(i);
        double xup  = xlow + hAll->GetBinWidth(i);
        double y    = hAll->GetBinContent(i);
        double err  = hAll->GetBinError(i);
        TBox* box = new TBox(xlow, y - err, xup, y + err);
        box->SetFillStyle(3003);
        box->SetFillColor(InteractionTypes[0].color);
        box->SetLineColor(InteractionTypes[0].color);
        box->Draw("SAME");
    }
    hAll->SetLineColor(InteractionTypes[0].color);
    hAll->SetLineWidth(2);
    hAll->Draw("HIST SAME");

    std::string hAllName = suffix + "_total";
    hAll->SetName(hAllName.c_str());
    parentDir->cd();
    hAll->Write();

    leg->SetTextSize(0.04);
    leg->Draw();
    parentDir->cd();
    c->Write();
    delete c;
}

// =============================================================================
// Main
// =============================================================================
void CC1e0piSelection_MakeAll()
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
    for (unsigned int iSel = 0; iSel < nTrig; ++iSel)
        sTrig[iSel] = new Spectrum("E_{#nu} [GeV]", Binning::Custom(TrueEnergyBinning),
            NuLoader,
            kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps[iSel].cut),
            kNoSpillCut);

    const unsigned int nNoTrig = SelectionSteps_NoTrigger.size();
    std::vector<Spectrum*> sNoTrig(nNoTrig, nullptr);
    for (unsigned int iSel = 0; iSel < nNoTrig; ++iSel)
        sNoTrig[iSel] = new Spectrum("E_{#nu} [GeV]", Binning::Custom(TrueEnergyBinning),
            NuLoader,
            kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(SelectionSteps_NoTrigger[iSel].cut),
            kNoSpillCut);

    // -------------------------------------------------------------------------
    // Purity spectra
    // -------------------------------------------------------------------------
    const unsigned int kNVar     = SelectionPlots.size();
    const unsigned int kNSelTrig = SelectionSteps.size();
    const unsigned int kNSelNoTrig = SelectionSteps_NoTrigger.size();

    std::vector<Spectrum*> spectra_Signal(kNVar, nullptr);
    std::vector<std::vector<Spectrum*>> spectra_AllTrig    (kNVar, std::vector<Spectrum*>(kNSelTrig,   nullptr));
    std::vector<std::vector<Spectrum*>> spectra_SignalTrig (kNVar, std::vector<Spectrum*>(kNSelTrig,   nullptr));
    std::vector<std::vector<Spectrum*>> spectra_AllNoTrig  (kNVar, std::vector<Spectrum*>(kNSelNoTrig, nullptr));
    std::vector<std::vector<Spectrum*>> spectra_SignalNoTrig(kNVar, std::vector<Spectrum*>(kNSelNoTrig, nullptr));

    for (unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        spectra_Signal[iVar] = new Spectrum(SelectionPlots[iVar].label,
            SelectionPlots[iVar].bins, NuLoader, SelectionPlots[iVar].var,
            kNoSpillCut, kTrueCC1e0pi);

        for (unsigned int jSel = 0; jSel < kNSelTrig; ++jSel) {
            spectra_AllTrig[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label,
                SelectionPlots[iVar].bins, NuLoader, SelectionPlots[iVar].var,
                kNoSpillCut, SelectionSteps[jSel].cut);
            spectra_SignalTrig[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label,
                SelectionPlots[iVar].bins, NuLoader, SelectionPlots[iVar].var,
                kNoSpillCut, SelectionSteps[jSel].cut && kTrueCC1e0pi);
        }

        for (unsigned int jSel = 0; jSel < kNSelNoTrig; ++jSel) {
            spectra_AllNoTrig[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label,
                SelectionPlots[iVar].bins, NuLoader, SelectionPlots[iVar].var,
                kNoSpillCut, SelectionSteps_NoTrigger[jSel].cut);
            spectra_SignalNoTrig[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label,
                SelectionPlots[iVar].bins, NuLoader, SelectionPlots[iVar].var,
                kNoSpillCut, SelectionSteps_NoTrigger[jSel].cut && kTrueCC1e0pi);
        }
    }

    // -------------------------------------------------------------------------
    // Selection spectra
    // -------------------------------------------------------------------------
    const unsigned int kNSel = InteractionTypes.size();
    std::vector<std::vector<Spectrum*>> spectra_Selection(kNVar, std::vector<Spectrum*>(kNSel, nullptr));

    for (unsigned int iVar = 0; iVar < kNVar; ++iVar)
        for (unsigned int jSel = 0; jSel < kNSel; ++jSel)
            spectra_Selection[iVar][jSel] = new Spectrum(SelectionPlots[iVar].label,
                SelectionPlots[iVar].bins, NuLoader, SelectionPlots[iVar].var,
                kNoSpillCut, kAutomaticSelection_NoTrigger && InteractionTypes[jSel].cut);

    // -------------------------------------------------------------------------
    // Ntuple - flat + per-PFP vector branches
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
    Tree ntuple("ntuple", branchNames, NuLoader, branchVars,
                kNoSpillCut, kNtuplePreSelection, kNoShift, true, true);

    // per-PFP vector branches
    Tree ntuplePFP("ntuple_pfp", vecBranchNames, NuLoader, vecBranchVars,
                   kNoSpillCut, kNtuplePreSelection, kNoShift, true, true);

    // -------------------------------------------------------------------------
    // Single Go() call — reads files once for everything
    // -------------------------------------------------------------------------
    NuLoader.Go();

    gROOT->SetBatch(kTRUE);
    SetPrettyStyle();

    TFile FOut("../output/root/CC1e0piSelection_All.root", "RECREATE");

    // =========================================================================
    // Efficiency outputs
    // =========================================================================
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

    DrawEfficiencySummary("efficiencyCC1e0pi_trigger",   "BNB CV (trigger)",
                          hTrue, sTrig,   SelectionSteps,          true,  FOut);
    DrawEfficiencySummary("efficiencyCC1e0pi_notrigger", "BNB CV (no trigger)",
                          hTrue, sNoTrig, SelectionSteps_NoTrigger, false, FOut);

    // =========================================================================
    // Purity outputs
    // =========================================================================
    TDirectory* purityDir = FOut.mkdir("purity");

    for (unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        const double pot = spectra_Signal[iVar]->POT();
        TH1* hSignal = spectra_Signal[iVar]->ToTH1(pot);
        hSignal->SetDirectory(nullptr);
        hSignal->SetName(("hSignal_" + SelectionPlots[iVar].suffix).c_str());
        hSignal->SetTitle((";" + SelectionPlots[iVar].label + ";Signal events").c_str());
        hSignal->SetMarkerStyle(20);
        hSignal->SetMarkerSize(1.0);
        hSignal->SetMarkerColor(TColor::GetColor("#111111"));
        hSignal->SetLineColor(TColor::GetColor("#888888"));
        hSignal->SetLineWidth(2);
        hSignal->GetXaxis()->SetTitleFont(132);
        hSignal->GetXaxis()->SetLabelFont(132);
        hSignal->GetYaxis()->SetTitleFont(132);
        hSignal->GetYaxis()->SetLabelFont(132);
        purityDir->cd();
        hSignal->Write();

        std::vector<TH1*> hAllTrig(kNSelTrig), hSelTrig(kNSelTrig);
        for (unsigned int jSel = 0; jSel < kNSelTrig; ++jSel) {
            hAllTrig[jSel] = spectra_AllTrig[iVar][jSel]->ToTH1(spectra_AllTrig[iVar][jSel]->POT());
            hAllTrig[jSel]->SetDirectory(nullptr);
            hSelTrig[jSel] = spectra_SignalTrig[iVar][jSel]->ToTH1(spectra_SignalTrig[iVar][jSel]->POT());
            hSelTrig[jSel]->SetDirectory(nullptr);
        }

        std::vector<TH1*> hAllNoTrig(kNSelNoTrig), hSelNoTrig(kNSelNoTrig);
        for (unsigned int jSel = 0; jSel < kNSelNoTrig; ++jSel) {
            hAllNoTrig[jSel] = spectra_AllNoTrig[iVar][jSel]->ToTH1(spectra_AllNoTrig[iVar][jSel]->POT());
            hAllNoTrig[jSel]->SetDirectory(nullptr);
            hSelNoTrig[jSel] = spectra_SignalNoTrig[iVar][jSel]->ToTH1(spectra_SignalNoTrig[iVar][jSel]->POT());
            hSelNoTrig[jSel]->SetDirectory(nullptr);
        }

        purityDir->cd();
        DrawPuritySummary("purityCC1e0pi_trigger_"   + SelectionPlots[iVar].suffix,
                          SelectionPlots[iVar], hSignal,
                          hAllTrig,   hSelTrig,   SelectionSteps,          true,  purityDir);
        DrawPuritySummary("purityCC1e0pi_notrigger_" + SelectionPlots[iVar].suffix,
                          SelectionPlots[iVar], hSignal,
                          hAllNoTrig, hSelNoTrig, SelectionSteps_NoTrigger, false, purityDir);
    }

    // =========================================================================
    // Selection outputs
    // =========================================================================
    TDirectory* selectionDir = FOut.mkdir("selection");
    const int nStackTypes = (int)kNSel - 1;
    std::vector<int> viridisCols = GetViridisColors(nStackTypes);

    for (unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        selectionDir->cd();
        DrawSelectionPlot(SelectionPlots[iVar], spectra_Selection[iVar], viridisCols, selectionDir);
    }

    // =========================================================================
    // Ntuple output
    // =========================================================================
    TDirectory* ntupleDir = FOut.mkdir("ntuple");
    ntuple.SaveTo(ntupleDir);

    TDirectory* ntuplePFPDir = FOut.mkdir("ntuple_pfp");
    ntuplePFP.SaveTo(ntuplePFPDir);

    FOut.Close();

    return;
}
