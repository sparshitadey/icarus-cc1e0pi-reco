#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <vector>
#include <math.h>

// local — must be included after all Vars headers
#include "Config.h"

namespace ana {

    // =========================================================================
    // SelDef struct — used by SelectionSteps and InteractionTypes
    // =========================================================================
    struct SelDef {
        std::string suffix = "";
        std::string label  = "";
        Cut cut   = kNoCut;
        int color = kBlack;
    };

    // =========================================================================
    // Spill-level cuts
    // =========================================================================
    const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* sr) {
        double minTime = 0., maxTime = 0.;
        for (const auto& match : sr->crtpmt_matches) {
            if ( sr->hdr.ismc) { minTime =  0.0; maxTime =  9.8; }
            if (!sr->hdr.ismc) { minTime = -0.4; maxTime = 10.5; }
            if (match.flashGateTime > minTime &&
                match.flashGateTime < maxTime &&
                match.flashClassification == 0)
                return true;
        }
        return false;
    });

    // =========================================================================
    // Slice-level pre-selection cuts
    // =========================================================================
    const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
        return !(slc->is_clear_cosmic);
    });

    const Cut kVertexInFV([](const caf::SRSliceProxy* slc) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        return kIsInFV(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    });

    // =========================================================================
    // Flash match cuts — use Config.h constants
    // =========================================================================
    const Cut kTrigFlashMatch([](const caf::SRSliceProxy* slc) {
        return (slc->barycenterFM.deltaZ_Trigger >= FM_DELTA_Z_MIN &&
                slc->barycenterFM.deltaZ_Trigger <= FM_DELTA_Z_MAX);
    });

    const Cut kFlashMatch([](const caf::SRSliceProxy* slc) {
        return (slc->barycenterFM.deltaZ    >= FM_DELTA_Z_MIN &&
                slc->barycenterFM.deltaZ    <= FM_DELTA_Z_MAX &&
                slc->barycenterFM.flashTime >  FM_FLASH_T_MIN &&
                slc->barycenterFM.flashTime <  FM_FLASH_T_MAX);
    });

    // =========================================================================
    // Muon veto
    // =========================================================================
    const Cut kMuonVeto([](const caf::SRSliceProxy* slc) {
        return !kHaveMuonCandidate(slc);
    });

    // =========================================================================
    // Electron identification cuts — use Config.h constants
    // =========================================================================
    const Cut kLargestRecoShower_EnergyCut([](const caf::SRSliceProxy* slc) {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return false;
        if (std::isnan(slc->reco.pfp[idx].shw.plane[2].energy)) return false;
        return slc->reco.pfp[idx].shw.plane[2].energy > SHOWER_ENERGY_MIN;
    });

    const Cut kLargestRecoShower_dEdxCut([](const caf::SRSliceProxy* slc) {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return false;
        if (std::isnan(slc->reco.pfp[idx].shw.plane[2].dEdx)) return false;
        return slc->reco.pfp[idx].shw.plane[2].dEdx < SHOWER_DEDX_MAX;
    });

    const Cut kLargestRecoShower_AvailabledEdxCut([](const caf::SRSliceProxy* slc) {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return false;
        if (std::isnan(kLargestRecoShower_AvailabledEdx(slc))) return false;
        return kLargestRecoShower_AvailabledEdx(slc) < SHOWER_DEDX_MAX;
    });

    const Cut kLargestRecoShower_OpenAngleCut([](const caf::SRSliceProxy* slc) {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return false;
        if (std::isnan(slc->reco.pfp[idx].shw.open_angle)) return false;
        double openAngle = 180. * slc->reco.pfp[idx].shw.open_angle / M_PI;
        return openAngle < SHOWER_ANGLE_MAX;
    });

    const Cut kLargestRecoShower_ConvGapCut([](const caf::SRSliceProxy* slc) {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return false;
        if (std::isnan(slc->reco.pfp[idx].shw.conversion_gap)) return false;
        return slc->reco.pfp[idx].shw.conversion_gap < SHOWER_CONVGAP_MAX;
    });

    // =========================================================================
    // Np0pi topology cuts
    // =========================================================================
    const Cut kNSelectedProtons([](const caf::SRSliceProxy* slc) {
        if (kLargestRecoShowerIdx(slc) == -1) return false;
        return kNSelectedProtonsIdx(slc).size() > 0;
    });

    const Cut kNoOtherParticle([](const caf::SRSliceProxy* slc) {
        if (kLargestRecoShowerIdx(slc) == -1) return false;
        int NOther = kNSelectedProtonsIdx_NOtherParticles(slc);
        if (NOther == -1) return false;
        return NOther == 0;
    });

    // =========================================================================
    // Combined selections
    // =========================================================================
    const Cut kAutomaticSelection =
        kNotClearCosmic && kVertexInFV && kTrigFlashMatch &&
        kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut &&
        kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut &&
        kNSelectedProtons && kNoOtherParticle && kMuonVeto;

    const Cut kAutomaticSelection_NoTrigger =
        kNotClearCosmic && kVertexInFV && kFlashMatch &&
        kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut &&
        kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut &&
        kNSelectedProtons && kNoOtherParticle && kMuonVeto;

    // loose pre-selection for ntuple output
    const Cut kNtuplePreSelection =
        kNotClearCosmic && kVertexInFV && kFlashMatch;

    const Cut kPreSelection =
        kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut;

    const Cut kPreSelection_NoTrigger =
        kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut;

    const Cut kCleaning =
        kNotClearCosmic && kVertexInFV && kTrigFlashMatch;

    // =========================================================================
    // Selection cutflow steps
    // =========================================================================
    std::vector<SelDef> SelectionSteps = {
        {"presel",      "Presel.",      kNotClearCosmic && kVertexInFV,                                                                                                                                                                                                         kBlack},
        {"flash",       "FM",           kNotClearCosmic && kVertexInFV && kTrigFlashMatch,                                                                                                                                                                                      kRed+2},
        {"shower",      "Electron ID",  kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut,                                                                                                                                                      kRed-7},
        {"showercuts1", "dE/dx",        kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut,                                                                                                                        kOrange-3},
        {"showercuts2", "Angle",        kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut,                                                                                     kGreen-2},
        {"showercuts3", "Gap",          kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut,                                                     kGreen-7},
        {"proton",      "Proton ID",    kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons,                                kCyan-7},
        {"nothingelse", "0#pi",         kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons && kNoOtherParticle,            kBlue-4},
        {"muonveto",    "LE-#mu veto",  kNotClearCosmic && kVertexInFV && kTrigFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons && kNoOtherParticle && kMuonVeto, kMagenta-3},
    };

    std::vector<SelDef> SelectionSteps_NoTrigger = {
        {"presel",      "Presel.",      kNotClearCosmic && kVertexInFV,                                                                                                                                                                                                         kBlack},
        {"flash",       "FM",           kNotClearCosmic && kVertexInFV && kFlashMatch,                                                                                                                                                                                          kRed+2},
        {"shower",      "Electron ID",  kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut,                                                                                                                                                          kRed-7},
        {"showercuts1", "dE/dx",        kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut,                                                                                                                            kOrange-3},
        {"showercuts2", "Angle, gap",   kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut,                                                        kGreen-2},
        {"proton",      "Proton ID",    kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons,                                   kCyan-7},
        {"nothingelse", "0#pi",         kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons && kNoOtherParticle,               kBlue-4},
        {"muonveto",    "LE-#mu veto",  kNotClearCosmic && kVertexInFV && kFlashMatch && kLargestRecoShower_EnergyCut && kLargestRecoShower_dEdxCut && kLargestRecoShower_OpenAngleCut && kLargestRecoShower_ConvGapCut && kNSelectedProtons && kNoOtherParticle && kMuonVeto,  kMagenta-3},
    };

} // namespace ana
