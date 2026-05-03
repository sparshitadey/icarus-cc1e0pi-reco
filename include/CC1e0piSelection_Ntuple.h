#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <vector>
#include <math.h>
#include <numeric>
#include <algorithm>

// root stuff
#include "TVector3.h"

// local — must be included after all Vars, Cuts, TruthCuts headers
#include "Config.h"

namespace ana {

    // =========================================================================
    // Helper: safe NaN check returning sentinel
    // =========================================================================
    inline double safe(double v, double sentinel = -5.) {
        return std::isnan(v) ? sentinel : v;
    }

    // =========================================================================
    // Truth variables
    // =========================================================================
    const Var kNtuple_IsSignal([](const caf::SRSliceProxy* slc) -> double {
        if (std::isnan(slc->truth.position.x) ||
            std::isnan(slc->truth.position.y) ||
            std::isnan(slc->truth.position.z)) return 0.;
        return kTrueCC1e0pi(slc) ? 1. : 0.;
    });

    const Var kNtuple_TrueNuEnergy([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0 || std::isnan(slc->truth.E)) return -5.;
        return slc->truth.E;
    });

    const Var kNtuple_TrueNuPDG([](const caf::SRSliceProxy* slc) -> double {
        return slc->truth.index >= 0 ? (double)slc->truth.pdg : -9999.;
    });

    const Var kNtuple_TrueIsCC([](const caf::SRSliceProxy* slc) -> double {
        return slc->truth.index >= 0 ? (slc->truth.iscc ? 1. : 0.) : -1.;
    });

    const Var kNtuple_TrueNuMode([](const caf::SRSliceProxy* slc) -> double {
        return slc->truth.index >= 0 ? (double)slc->truth.genie_mode : -9999.;
    });

    const Var kNtuple_TrueNuIndex([](const caf::SRSliceProxy* slc) -> double {
        return slc->truth.index;
    });

    const Var kNtuple_TrueVtxX([](const caf::SRSliceProxy* slc) -> double { return safe(slc->truth.position.x, -9999.); });
    const Var kNtuple_TrueVtxY([](const caf::SRSliceProxy* slc) -> double { return safe(slc->truth.position.y, -9999.); });
    const Var kNtuple_TrueVtxZ([](const caf::SRSliceProxy* slc) -> double { return safe(slc->truth.position.z, -9999.); });

    const Var kNtuple_TrueQ2([](const caf::SRSliceProxy* slc) -> double {
        return slc->truth.index >= 0 ? safe(slc->truth.Q2) : -5.;
    });
    const Var kNtuple_TrueW([](const caf::SRSliceProxy* slc) -> double {
        return -5.; // slc->truth.W not available in this SRProxy version
    });
    const Var kNtuple_TrueX([](const caf::SRSliceProxy* slc) -> double {
        return -5.; // slc->truth.bjorkX not available in this SRProxy version
    });
    const Var kNtuple_TrueY([](const caf::SRSliceProxy* slc) -> double {
        return -5.; // slc->truth.bjorkY not available in this SRProxy version
    });
    const Var kNtuple_TruePt([](const caf::SRSliceProxy* slc) -> double {
        return -5.; // slc->truth.pt not available in this SRProxy version
    });
    const Var kNtuple_TruePtX([](const caf::SRSliceProxy* slc) -> double {
        return -9999.; // slc->truth.momentum not available in this SRProxy version
    });
    const Var kNtuple_TruePtY([](const caf::SRSliceProxy* slc) -> double {
        return -9999.; // slc->truth.momentum not available in this SRProxy version
    });

    const Var kNtuple_TrueElectronKE([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -5.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 11)
                return slc->truth.prim[ip].startE - slc->truth.prim[ip].endE;
        return -5.;
    });

    const Var kNtuple_TrueElectronDirX([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -9999.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 11) return safe(slc->truth.prim[ip].start.x, -9999.);
        return -9999.;
    });
    const Var kNtuple_TrueElectronDirY([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -9999.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 11) return safe(slc->truth.prim[ip].start.y, -9999.);
        return -9999.;
    });
    const Var kNtuple_TrueElectronDirZ([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -9999.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 11) return safe(slc->truth.prim[ip].start.z, -9999.);
        return -9999.;
    });

    const Var kNtuple_TrueNVisProtons([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 2212 &&
                (slc->truth.prim[ip].startE - slc->truth.prim[ip].endE) > VISIBILTY_THRESHOLD_P &&
                kIsInContained(slc->truth.prim[ip].end.x, slc->truth.prim[ip].end.y, slc->truth.prim[ip].end.z))
                ++n;
        return n;
    });

    const Var kNtuple_TrueLeadingProtonKE([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -5.;
        double maxKE = -5.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 2212) {
                double ke = slc->truth.prim[ip].startE - slc->truth.prim[ip].endE;
                if (ke > maxKE) maxKE = ke;
            }
        return maxKE;
    });

    const Var kNtuple_TrueSubLeadingProtonKE([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -5.;
        std::vector<double> kes;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 2212)
                kes.push_back(slc->truth.prim[ip].startE - slc->truth.prim[ip].endE);
        if (kes.size() < 2) return -5.;
        std::sort(kes.begin(), kes.end(), std::greater<double>());
        return kes[1];
    });

    const Var kNtuple_TrueNSubThresholdProtons([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 2212) {
                double ke = slc->truth.prim[ip].startE - slc->truth.prim[ip].endE;
                if (ke > 0 && ke <= VISIBILTY_THRESHOLD_P) ++n;
            }
        return n;
    });

    const Var kNtuple_TrueNSubThresholdPions([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 211) {
                double ke = slc->truth.prim[ip].startE - slc->truth.prim[ip].endE;
                if (ke > 0 && ke <= VISIBILTY_THRESHOLD_PI) ++n;
            }
        return n;
    });

    const Var kNtuple_TrueNNeutrons([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 2112) ++n;
        return n;
    });

    const Var kNtuple_TrueNPi0([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 111) ++n;
        return n;
    });

    const Var kNtuple_TrueNPhotons([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 22) ++n;
        return n;
    });

    const Var kNtuple_TrueNVisPions([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -1.;
        int n = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 211 &&
                (slc->truth.prim[ip].startE - slc->truth.prim[ip].endE) >= VISIBILTY_THRESHOLD_PI)
                ++n;
        return n;
    });

    const Var kNtuple_TrueLeadingPionKE([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -5.;
        double maxKE = -5.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (abs(slc->truth.prim[ip].pdg) == 211) {
                double ke = slc->truth.prim[ip].startE - slc->truth.prim[ip].endE;
                if (ke > maxKE) maxKE = ke;
            }
        return maxKE;
    });

    const Var kNtuple_TrueLeadingPhotonE([](const caf::SRSliceProxy* slc) -> double {
        if (slc->truth.index < 0) return -5.;
        double maxE = -5.;
        for (int ip(0); ip < slc->truth.nprim; ++ip)
            if (slc->truth.prim[ip].pdg == 22 && slc->truth.prim[ip].startE > maxE)
                maxE = slc->truth.prim[ip].startE;
        return maxE;
    });

    // =========================================================================
    // Reconstructed vertex
    // =========================================================================
    const Var kNtuple_RecoVtxX([](const caf::SRSliceProxy* slc) -> double { return safe(slc->vertex.x, -9999.); });
    const Var kNtuple_RecoVtxY([](const caf::SRSliceProxy* slc) -> double { return safe(slc->vertex.y, -9999.); });
    const Var kNtuple_RecoVtxZ([](const caf::SRSliceProxy* slc) -> double { return safe(slc->vertex.z, -9999.); });
    const Var kNtuple_RecoVtxVsTruth([](const caf::SRSliceProxy* slc) -> double { return kVertex_vsTruth(slc); });
    const Var kNtuple_RecoVtxInFV([](const caf::SRSliceProxy* slc) -> double {
        if (std::isnan(slc->vertex.x)) return -1.;
        return kIsInFV(slc->vertex.x, slc->vertex.y, slc->vertex.z) ? 1. : 0.;
    });

    // =========================================================================
    // Shower variables — all planes, all raw values
    // =========================================================================
    const Var kNtuple_ShowerEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.plane[2].energy);
    });
    const Var kNtuple_ShowerEnergyPlane0([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.plane[0].energy);
    });
    const Var kNtuple_ShowerEnergyPlane1([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.plane[1].energy);
    });
    const Var kNtuple_ShowerdEdx([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.plane[2].dEdx);
    });
    const Var kNtuple_ShowerdEdxPlane0([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.plane[0].dEdx);
    });
    const Var kNtuple_ShowerdEdxPlane1([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.plane[1].dEdx);
    });
    const Var kNtuple_ShowerAvailabledEdx([](const caf::SRSliceProxy* slc) -> double {
        return kLargestRecoShower_AvailabledEdx(slc);
    });
    const Var kNtuple_ShowerNHitsPlane0([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return slc->reco.pfp[idx].shw.plane[0].nHits;
    });
    const Var kNtuple_ShowerNHitsPlane1([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return slc->reco.pfp[idx].shw.plane[1].nHits;
    });
    const Var kNtuple_ShowerNHitsPlane2([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return slc->reco.pfp[idx].shw.plane[2].nHits;
    });
    const Var kNtuple_ShowerOpenAngle([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].shw.open_angle)) return -5.;
        return 180. * slc->reco.pfp[idx].shw.open_angle / M_PI;
    });
    const Var kNtuple_ShowerConvGap([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.conversion_gap);
    });
    const Var kNtuple_ShowerLength([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.len);
    });
    const Var kNtuple_ShowerTrackScore([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].trackScore);
    });
    const Var kNtuple_ShowerStartX([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -9999.;
        return safe(slc->reco.pfp[idx].shw.start.x, -9999.);
    });
    const Var kNtuple_ShowerStartY([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -9999.;
        return safe(slc->reco.pfp[idx].shw.start.y, -9999.);
    });
    const Var kNtuple_ShowerStartZ([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -9999.;
        return safe(slc->reco.pfp[idx].shw.start.z, -9999.);
    });
    const Var kNtuple_ShowerDirX([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -9999.;
        return safe(slc->reco.pfp[idx].shw.dir.x, -9999.);
    });
    const Var kNtuple_ShowerDirY([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -9999.;
        return safe(slc->reco.pfp[idx].shw.dir.y, -9999.);
    });
    const Var kNtuple_ShowerDirZ([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -9999.;
        return safe(slc->reco.pfp[idx].shw.dir.z, -9999.);
    });
    const Var kNtuple_ShowerStartVsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].shw.start.x) || std::isnan(slc->truth.position.x)) return -5.;
        TVector3 shwStart(slc->reco.pfp[idx].shw.start.x, slc->reco.pfp[idx].shw.start.y, slc->reco.pfp[idx].shw.start.z);
        TVector3 trueVtx(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
        return (shwStart - trueVtx).Mag();
    });
    const Var kNtuple_ShowerTruePDG([](const caf::SRSliceProxy* slc) -> double { return kLargestRecoShower_TruePdg(slc); });
    const Var kNtuple_ShowerHitPurity([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.truth.pur);
    });
    const Var kNtuple_ShowerHitEfficiency([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.truth.eff);
    });
    const Var kNtuple_ShowerTrueStartEndEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].shw.truth.p.startE - slc->reco.pfp[idx].shw.truth.p.endE);
    });


#ifdef USE_NUGRAPH
    // NuGraph scores on leading shower
    const Var kNtuple_ShowerNGShrFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.shr_frac);
    });
    const Var kNtuple_ShowerNGHipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.hip_frac);
    });
    const Var kNtuple_ShowerNGMipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.mip_frac);
    });
    const Var kNtuple_ShowerNGMhlFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.mhl_frac);
    });
    const Var kNtuple_ShowerNGDifFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.dif_frac);
    });
    const Var kNtuple_ShowerNGBkgFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.bkg_frac);
    });
    const Var kNtuple_ShowerNGSemCat([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc); if (idx == -1) return -5.;
        return safe(slc->reco.pfp[idx].ngscore.sem_cat);
    });

#endif // USE_NUGRAPH
    // =========================================================================
    // Flash match
    // =========================================================================
    const Var kNtuple_FM_DeltaZ([](const caf::SRSliceProxy* slc) -> double { return slc->barycenterFM.deltaZ; });
    const Var kNtuple_FM_DeltaZ_Trigger([](const caf::SRSliceProxy* slc) -> double { return slc->barycenterFM.deltaZ_Trigger; });
    const Var kNtuple_FM_FlashTime([](const caf::SRSliceProxy* slc) -> double { return slc->barycenterFM.flashTime; });

    // =========================================================================
    // Proton variables
    // =========================================================================
    const Var kNtuple_NProtons([](const caf::SRSliceProxy* slc) -> double {
        return (double)kNSelectedProtonsIdx(slc).size();
    });
    const Var kNtuple_NOtherParticles([](const caf::SRSliceProxy* slc) -> double {
        int n = kNSelectedProtonsIdx_NOtherParticles(slc);
        return n == -1 ? -5. : (double)n;
    });
    const Var kNtuple_LeadingProtonMomentum([](const caf::SRSliceProxy* slc) -> double { return kLeadingProtonMomentum(slc); });
    const Var kNtuple_SubLeadingProtonMomentum([](const caf::SRSliceProxy* slc) -> double { return kSubLeadingProtonMomentum(slc); });

    // helper lambda to get leading proton index by momentum
    auto getLeadProtonIdx = [](const caf::SRSliceProxy* slc) -> int {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -1;
        double maxMom = 0.; int leadIdx = -1;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); leadIdx = i; }
        }
        return leadIdx;
    };

    const Var kNtuple_LeadingProtonChi2Proton([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = slc->reco.pfp[i].trk.chi2pid[2].chi2_proton; }
        }
        return val;
    });

    const Var kNtuple_LeadingProtonChi2Muon([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = slc->reco.pfp[i].trk.chi2pid[2].chi2_muon; }
        }
        return val;
    });

    const Var kNtuple_LeadingProtonLength([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = safe(slc->reco.pfp[i].trk.len); }
        }
        return val;
    });


#ifdef USE_NUGRAPH
    const Var kNtuple_LeadingProtonNGHipFrac([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = safe(slc->reco.pfp[i].ngscore.hip_frac); }
        }
        return val;
    });

    const Var kNtuple_LeadingProtonNGMipFrac([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = safe(slc->reco.pfp[i].ngscore.mip_frac); }
        }
        return val;
    });

    const Var kNtuple_LeadingProtonNGSemCat([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = safe(slc->reco.pfp[i].ngscore.sem_cat); }
        }
        return val;
    });

#endif // USE_NUGRAPH

    const Var kNtuple_LeadingProtonTruePDG([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -9999.;
        double maxMom = 0.; double val = -9999.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = (double)slc->reco.pfp[i].trk.truth.p.pdg; }
        }
        return val;
    });

    const Var kNtuple_LeadingProtonHitPurity([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = safe(slc->reco.pfp[i].trk.truth.pur); }
        }
        return val;
    });

    const Var kNtuple_LeadingProtonHitEfficiency([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double maxMom = 0.; double val = -5.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); val = safe(slc->reco.pfp[i].trk.truth.eff); }
        }
        return val;
    });

    // =========================================================================
    // Muon veto raw variables
    // =========================================================================
    const Var kNtuple_LongestTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
        double highestLength = 0.; double val = -5.;
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.len)) continue;
            if (slc->reco.pfp[i].trackScore < 0.5) continue;
            if (slc->reco.pfp[i].trk.len > highestLength && slc->reco.pfp[i].parent_is_primary)
                { highestLength = slc->reco.pfp[i].trk.len; val = slc->reco.pfp[i].trk.chi2pid[2].chi2_muon; }
        }
        return val;
    });
    const Var kNtuple_LongestTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
        double highestLength = 0.; double val = -5.;
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.len)) continue;
            if (slc->reco.pfp[i].trackScore < 0.5) continue;
            if (slc->reco.pfp[i].trk.len > highestLength && slc->reco.pfp[i].parent_is_primary)
                { highestLength = slc->reco.pfp[i].trk.len; val = slc->reco.pfp[i].trk.chi2pid[2].chi2_proton; }
        }
        return val;
    });
    const Var kNtuple_LongestTrackLength([](const caf::SRSliceProxy* slc) -> double {
        double val = -5.;
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.len)) continue;
            if (slc->reco.pfp[i].trackScore < 0.5) continue;
            if (slc->reco.pfp[i].trk.len > val && slc->reco.pfp[i].parent_is_primary)
                val = slc->reco.pfp[i].trk.len;
        }
        return val;
    });
    const Var kNtuple_LongestTrackScore([](const caf::SRSliceProxy* slc) -> double {
        double highestLength = 0.; double val = -5.;
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.len)) continue;
            if (slc->reco.pfp[i].trk.len > highestLength && slc->reco.pfp[i].parent_is_primary)
                { highestLength = slc->reco.pfp[i].trk.len; val = slc->reco.pfp[i].trackScore; }
        }
        return val;
    });

#ifdef USE_NUGRAPH
#ifdef USE_NUGRAPH
    const Var kNtuple_LongestTrackNGSemCat([](const caf::SRSliceProxy* slc) -> double {
        double highestLength = 0.; double val = -5.;
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.len)) continue;
            if (slc->reco.pfp[i].trk.len > highestLength && slc->reco.pfp[i].parent_is_primary)
                { highestLength = slc->reco.pfp[i].trk.len; val = safe(slc->reco.pfp[i].ngscore.sem_cat); }
        }
        return val;
    });
#endif // USE_NUGRAPH
#endif // USE_NUGRAPH

    // =========================================================================
    // Reconstructed neutrino kinematics
    // =========================================================================
    const Var kNtuple_RecoNuEnergy([](const caf::SRSliceProxy* slc) -> double { return kRecoNeutrino_CC0piEnergy(slc); });
    const Var kNtuple_RecoNuEnergyVsTruth([](const caf::SRSliceProxy* slc) -> double { return kRecoNeutrino_CC0piEnergy_VsTruth(slc); });
    const Var kNtuple_RecoHadronEnergy([](const caf::SRSliceProxy* slc) -> double { return kRecoNeutrino_CC0piHadronEnergy(slc); });
    const Var kNtuple_RecoInelasticity([](const caf::SRSliceProxy* slc) -> double { return kRecoNeutrino_CC0piInelasticity(slc); });
    const Var kNtuple_RecoTransverseMomentum([](const caf::SRSliceProxy* slc) -> double { return kRecoNeutrino_CC0piTransverseMomentum(slc); });

    const Var kNtuple_RecoOpeningAngleEP([](const caf::SRSliceProxy* slc) -> double {
        const int shwIdx = kLargestRecoShowerIdx(slc);
        if (shwIdx == -1) return -5.;
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        if (std::isnan(slc->reco.pfp[shwIdx].shw.dir.x)) return -5.;
        TVector3 eDir(slc->reco.pfp[shwIdx].shw.dir.x, slc->reco.pfp[shwIdx].shw.dir.y, slc->reco.pfp[shwIdx].shw.dir.z);
        double maxMom = 0.; int leadP = -1;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) continue;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            if (p.Mag() > maxMom) { maxMom = p.Mag(); leadP = i; }
        }
        if (leadP == -1 || std::isnan(slc->reco.pfp[leadP].trk.dir.x)) return -5.;
        TVector3 pDir(slc->reco.pfp[leadP].trk.dir.x, slc->reco.pfp[leadP].trk.dir.y, slc->reco.pfp[leadP].trk.dir.z);
        return 180. * eDir.Angle(pDir) / M_PI;
    });

    // =========================================================================
    // Slice-level variables
    // =========================================================================
    const Var kNtuple_SliceNPFPs([](const caf::SRSliceProxy* slc) -> double { return slc->reco.npfp; });
    const Var kNtuple_SliceIsClearCosmic([](const caf::SRSliceProxy* slc) -> double { return slc->is_clear_cosmic ? 1. : 0.; });
    const Var kNtuple_SliceNuScore([](const caf::SRSliceProxy* slc) -> double { return safe(slc->nu_score); });
    const Var kNtuple_SliceCrumbsScore([](const caf::SRSliceProxy* slc) -> double { return safe(slc->crumbs_result.score); });


#ifdef USE_NUGRAPH
    // NuGraph event-level
    const Var kNtuple_NGFilterFraction([](const caf::SRSliceProxy* slc) -> double { return safe(slc->ng_filt_pass_frac); });
    const Var kNtuple_NGVtxHipHitsInd1([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[0].ng_vtx_hip_hits; });
    const Var kNtuple_NGVtxHipHitsInd2([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[1].ng_vtx_hip_hits; });
    const Var kNtuple_NGVtxHipHitsColl([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[2].ng_vtx_hip_hits; });
    const Var kNtuple_NGShrHitsInd1([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[0].shr_hits; });
    const Var kNtuple_NGShrHitsInd2([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[1].shr_hits; });
    const Var kNtuple_NGShrHitsColl([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[2].shr_hits; });
    const Var kNtuple_NGUnclShrHitsInd1([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[0].unclustered_shr_hits; });
    const Var kNtuple_NGUnclShrHitsInd2([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[1].unclustered_shr_hits; });
    const Var kNtuple_NGUnclShrHitsColl([](const caf::SRSliceProxy* slc) -> double { return slc->ng_plane[2].unclustered_shr_hits; });

#endif // USE_NUGRAPH
    // =========================================================================
    // Fixed cut pass/fail flags
    // =========================================================================
    const Var kNtuple_PassNotClearCosmic([](const caf::SRSliceProxy* slc) -> double { return kNotClearCosmic(slc) ? 1. : 0.; });
    const Var kNtuple_PassVertexInFV([](const caf::SRSliceProxy* slc) -> double { return kVertexInFV(slc) ? 1. : 0.; });
    const Var kNtuple_PassFlashMatch([](const caf::SRSliceProxy* slc) -> double { return kFlashMatch(slc) ? 1. : 0.; });
    const Var kNtuple_PassTrigFlashMatch([](const caf::SRSliceProxy* slc) -> double { return kTrigFlashMatch(slc) ? 1. : 0.; });

    // =========================================================================
    // Per-PFP vector branches — full NuGraph vs Pandora comparison
    // =========================================================================
    const MultiVar kNtuple_PFP_TrackScore([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].trackScore)); return v;
    });
    const MultiVar kNtuple_PFP_IsTrack([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(slc->reco.pfp[i].trackScore >= 0.5 ? 1. : 0.); return v;
    });
    const MultiVar kNtuple_PFP_ParentIsPrimary([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(slc->reco.pfp[i].parent_is_primary ? 1. : 0.); return v;
    });

#ifdef USE_NUGRAPH
    const MultiVar kNtuple_PFP_NGSemCat([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.sem_cat)); return v;
    });
    const MultiVar kNtuple_PFP_NGShrFrac([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.shr_frac)); return v;
    });
    const MultiVar kNtuple_PFP_NGHipFrac([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.hip_frac)); return v;
    });
    const MultiVar kNtuple_PFP_NGMipFrac([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.mip_frac)); return v;
    });
    const MultiVar kNtuple_PFP_NGMhlFrac([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.mhl_frac)); return v;
    });
    const MultiVar kNtuple_PFP_NGDifFrac([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.dif_frac)); return v;
    });
    const MultiVar kNtuple_PFP_NGBkgFrac([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].ngscore.bkg_frac)); return v;
    });
#endif // USE_NUGRAPH

        const MultiVar kNtuple_PFP_NHitsPlane0([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(slc->reco.pfp[i].shw.plane[0].nHits); return v;
    });
    const MultiVar kNtuple_PFP_NHitsPlane1([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(slc->reco.pfp[i].shw.plane[1].nHits); return v;
    });
    const MultiVar kNtuple_PFP_NHitsPlane2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(slc->reco.pfp[i].shw.plane[2].nHits); return v;
    });
    const MultiVar kNtuple_PFP_ShwEnergyPlane2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.plane[2].energy)); return v;
    });
    const MultiVar kNtuple_PFP_ShwDedxPlane2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.plane[2].dEdx)); return v;
    });
    const MultiVar kNtuple_PFP_ShwOpenAngle([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v;
        for (unsigned int i = 0; i < slc->reco.npfp; i++)
            v.push_back(std::isnan(slc->reco.pfp[i].shw.open_angle) ? -5. : 180. * slc->reco.pfp[i].shw.open_angle / M_PI);
        return v;
    });
    const MultiVar kNtuple_PFP_ShwConvGap([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.conversion_gap)); return v;
    });
    const MultiVar kNtuple_PFP_ShwLength([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.len)); return v;
    });
    const MultiVar kNtuple_PFP_TrkChi2Proton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].trk.chi2pid[2].chi2_proton)); return v;
    });
    const MultiVar kNtuple_PFP_TrkChi2Muon([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].trk.chi2pid[2].chi2_muon)); return v;
    });
    const MultiVar kNtuple_PFP_TrkLength([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].trk.len)); return v;
    });
    const MultiVar kNtuple_PFP_TrkRangeMomProton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].trk.rangeP.p_proton)); return v;
    });
    const MultiVar kNtuple_PFP_TrkRangeMomPion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].trk.rangeP.p_pion)); return v;
    });
    const MultiVar kNtuple_PFP_TruePDG([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back((double)slc->reco.pfp[i].shw.truth.p.pdg); return v;
    });
    const MultiVar kNtuple_PFP_HitPurity([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.truth.pur)); return v;
    });
    const MultiVar kNtuple_PFP_HitEfficiency([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.truth.eff)); return v;
    });
    const MultiVar kNtuple_PFP_TrueStartEndEnergy([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v;
        for (unsigned int i = 0; i < slc->reco.npfp; i++)
            v.push_back(safe(slc->reco.pfp[i].shw.truth.p.startE - slc->reco.pfp[i].shw.truth.p.endE));
        return v;
    });
    const MultiVar kNtuple_PFP_TrueDirX([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.truth.p.genp.x, -9999.)); return v;
    });
    const MultiVar kNtuple_PFP_TrueDirY([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.truth.p.genp.y, -9999.)); return v;
    });
    const MultiVar kNtuple_PFP_TrueDirZ([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> v; for (unsigned int i = 0; i < slc->reco.npfp; i++) v.push_back(safe(slc->reco.pfp[i].shw.truth.p.genp.z, -9999.)); return v;
    });

    // =========================================================================
    // Branch structs
    // =========================================================================
    struct NtupleBranch {
        std::string name;
        Var         var;
    };

    struct NtupleVecBranch {
        std::string name;
        MultiVar    var;
    };

    // flat branches
    std::vector<NtupleBranch> NtupleBranches = {
        // truth
        {"is_signal",                    kNtuple_IsSignal},
        {"true_nu_energy",               kNtuple_TrueNuEnergy},
        {"true_nu_pdg",                  kNtuple_TrueNuPDG},
        {"true_iscc",                    kNtuple_TrueIsCC},
        {"true_nu_mode",                 kNtuple_TrueNuMode},
        {"true_nu_index",                kNtuple_TrueNuIndex},
        {"true_vtx_x",                   kNtuple_TrueVtxX},
        {"true_vtx_y",                   kNtuple_TrueVtxY},
        {"true_vtx_z",                   kNtuple_TrueVtxZ},
        {"true_q2",                      kNtuple_TrueQ2},
        {"true_w",                       kNtuple_TrueW},
        {"true_bjorken_x",               kNtuple_TrueX},
        {"true_bjorken_y",               kNtuple_TrueY},
        {"true_pt",                      kNtuple_TruePt},
        {"true_pt_x",                    kNtuple_TruePtX},
        {"true_pt_y",                    kNtuple_TruePtY},
        {"true_electron_ke",             kNtuple_TrueElectronKE},
        {"true_electron_dir_x",          kNtuple_TrueElectronDirX},
        {"true_electron_dir_y",          kNtuple_TrueElectronDirY},
        {"true_electron_dir_z",          kNtuple_TrueElectronDirZ},
        {"true_n_vis_protons",           kNtuple_TrueNVisProtons},
        {"true_leading_proton_ke",       kNtuple_TrueLeadingProtonKE},
        {"true_subleading_proton_ke",    kNtuple_TrueSubLeadingProtonKE},
        {"true_n_subthreshold_protons",  kNtuple_TrueNSubThresholdProtons},
        {"true_n_subthreshold_pions",    kNtuple_TrueNSubThresholdPions},
        {"true_n_neutrons",              kNtuple_TrueNNeutrons},
        {"true_n_pi0",                   kNtuple_TrueNPi0},
        {"true_n_photons",               kNtuple_TrueNPhotons},
        {"true_n_vis_pions",             kNtuple_TrueNVisPions},
        {"true_leading_pion_ke",         kNtuple_TrueLeadingPionKE},
        {"true_leading_photon_e",        kNtuple_TrueLeadingPhotonE},
        // reco vertex
        {"reco_vtx_x",                   kNtuple_RecoVtxX},
        {"reco_vtx_y",                   kNtuple_RecoVtxY},
        {"reco_vtx_z",                   kNtuple_RecoVtxZ},
        {"reco_vtx_vs_truth",            kNtuple_RecoVtxVsTruth},
        {"reco_vtx_in_fv",               kNtuple_RecoVtxInFV},
        // shower
        {"shower_energy",                kNtuple_ShowerEnergy},
        {"shower_energy_plane0",         kNtuple_ShowerEnergyPlane0},
        {"shower_energy_plane1",         kNtuple_ShowerEnergyPlane1},
        {"shower_dedx",                  kNtuple_ShowerdEdx},
        {"shower_dedx_plane0",           kNtuple_ShowerdEdxPlane0},
        {"shower_dedx_plane1",           kNtuple_ShowerdEdxPlane1},
        {"shower_avail_dedx",            kNtuple_ShowerAvailabledEdx},
        {"shower_nhits_plane0",          kNtuple_ShowerNHitsPlane0},
        {"shower_nhits_plane1",          kNtuple_ShowerNHitsPlane1},
        {"shower_nhits_plane2",          kNtuple_ShowerNHitsPlane2},
        {"shower_open_angle",            kNtuple_ShowerOpenAngle},
        {"shower_conv_gap",              kNtuple_ShowerConvGap},
        {"shower_length",                kNtuple_ShowerLength},
        {"shower_track_score",           kNtuple_ShowerTrackScore},
        {"shower_start_x",               kNtuple_ShowerStartX},
        {"shower_start_y",               kNtuple_ShowerStartY},
        {"shower_start_z",               kNtuple_ShowerStartZ},
        {"shower_dir_x",                 kNtuple_ShowerDirX},
        {"shower_dir_y",                 kNtuple_ShowerDirY},
        {"shower_dir_z",                 kNtuple_ShowerDirZ},
        {"shower_start_vs_truth",        kNtuple_ShowerStartVsTruth},
        {"shower_true_pdg",              kNtuple_ShowerTruePDG},
        {"shower_hit_purity",            kNtuple_ShowerHitPurity},
        {"shower_hit_efficiency",      kNtuple_ShowerHitEfficiency},
        {"shower_true_start_end_energy", kNtuple_ShowerTrueStartEndEnergy},
#ifdef USE_NUGRAPH
        {"shower_ng_shr_frac",           kNtuple_ShowerNGShrFrac},
        {"shower_ng_hip_frac",           kNtuple_ShowerNGHipFrac},
        {"shower_ng_mip_frac",           kNtuple_ShowerNGMipFrac},
        {"shower_ng_mhl_frac",           kNtuple_ShowerNGMhlFrac},
        {"shower_ng_dif_frac",           kNtuple_ShowerNGDifFrac},
        {"shower_ng_bkg_frac",           kNtuple_ShowerNGBkgFrac},
        {"shower_ng_sem_cat",            kNtuple_ShowerNGSemCat},
#endif // USE_NUGRAPH
                // flash match
        {"fm_delta_z",                   kNtuple_FM_DeltaZ},
        {"fm_delta_z_trigger",           kNtuple_FM_DeltaZ_Trigger},
        {"fm_flash_time",                kNtuple_FM_FlashTime},
        // protons
        {"n_protons",                    kNtuple_NProtons},
        {"n_other_particles",            kNtuple_NOtherParticles},
        {"lead_proton_momentum",         kNtuple_LeadingProtonMomentum},
        {"sublead_proton_momentum",      kNtuple_SubLeadingProtonMomentum},
        {"lead_proton_chi2_proton",      kNtuple_LeadingProtonChi2Proton},
        {"lead_proton_chi2_muon",        kNtuple_LeadingProtonChi2Muon},
        {"lead_proton_length",           kNtuple_LeadingProtonLength},
#ifdef USE_NUGRAPH
        {"lead_proton_ng_hip_frac",      kNtuple_LeadingProtonNGHipFrac},
        {"lead_proton_ng_mip_frac",      kNtuple_LeadingProtonNGMipFrac},
        {"lead_proton_ng_sem_cat",       kNtuple_LeadingProtonNGSemCat},
#endif // USE_NUGRAPH
                {"lead_proton_true_pdg",         kNtuple_LeadingProtonTruePDG},
        {"lead_proton_hit_purity",       kNtuple_LeadingProtonHitPurity},
        {"lead_proton_hit_efficiency", kNtuple_LeadingProtonHitEfficiency},
        // muon veto
        {"longest_trk_chi2_muon",        kNtuple_LongestTrackChi2Muon},
        {"longest_trk_chi2_proton",      kNtuple_LongestTrackChi2Proton},
        {"longest_trk_length",           kNtuple_LongestTrackLength},
        {"longest_trk_track_score",      kNtuple_LongestTrackScore},
#ifdef USE_NUGRAPH
        {"longest_trk_ng_sem_cat",       kNtuple_LongestTrackNGSemCat},
#endif // USE_NUGRAPH
        // reco kinematics
        {"reco_nu_energy",               kNtuple_RecoNuEnergy},
        {"reco_nu_energy_vs_truth",      kNtuple_RecoNuEnergyVsTruth},
        {"reco_hadron_energy",           kNtuple_RecoHadronEnergy},
        {"reco_inelasticity",            kNtuple_RecoInelasticity},
        {"reco_transverse_momentum",     kNtuple_RecoTransverseMomentum},
        {"reco_opening_angle_ep",        kNtuple_RecoOpeningAngleEP},
        // slice
        {"slice_n_pfps",                 kNtuple_SliceNPFPs},
        {"slice_is_clear_cosmic",        kNtuple_SliceIsClearCosmic},
        {"slice_nu_score",               kNtuple_SliceNuScore},
        {"slice_crumbs_score",           kNtuple_SliceCrumbsScore},
#ifdef USE_NUGRAPH
        // NuGraph event-level
        {"ng_filter_fraction",           kNtuple_NGFilterFraction},
        {"ng_vtx_hip_hits_ind1",         kNtuple_NGVtxHipHitsInd1},
        {"ng_vtx_hip_hits_ind2",         kNtuple_NGVtxHipHitsInd2},
        {"ng_vtx_hip_hits_coll",         kNtuple_NGVtxHipHitsColl},
        {"ng_shr_hits_ind1",             kNtuple_NGShrHitsInd1},
        {"ng_shr_hits_ind2",             kNtuple_NGShrHitsInd2},
        {"ng_shr_hits_coll",             kNtuple_NGShrHitsColl},
        {"ng_uncl_shr_hits_ind1",        kNtuple_NGUnclShrHitsInd1},
        {"ng_uncl_shr_hits_ind2",        kNtuple_NGUnclShrHitsInd2},
        {"ng_uncl_shr_hits_coll",        kNtuple_NGUnclShrHitsColl},
#endif // USE_NUGRAPH
                // fixed cut flags
        {"pass_not_clear_cosmic",        kNtuple_PassNotClearCosmic},
        {"pass_vertex_fv",               kNtuple_PassVertexInFV},
        {"pass_flash_match",             kNtuple_PassFlashMatch},
        {"pass_trig_flash_match",        kNtuple_PassTrigFlashMatch},
    };

    // vector (per-PFP) branches
    std::vector<NtupleVecBranch> NtupleVecBranches = {
        {"pfp_track_score",           kNtuple_PFP_TrackScore},
        {"pfp_is_track",              kNtuple_PFP_IsTrack},
        {"pfp_parent_is_primary",     kNtuple_PFP_ParentIsPrimary},
#ifdef USE_NUGRAPH
        {"pfp_ng_sem_cat",            kNtuple_PFP_NGSemCat},
        {"pfp_ng_shr_frac",           kNtuple_PFP_NGShrFrac},
        {"pfp_ng_hip_frac",           kNtuple_PFP_NGHipFrac},
        {"pfp_ng_mip_frac",           kNtuple_PFP_NGMipFrac},
        {"pfp_ng_mhl_frac",           kNtuple_PFP_NGMhlFrac},
        {"pfp_ng_dif_frac",           kNtuple_PFP_NGDifFrac},
        {"pfp_ng_bkg_frac",           kNtuple_PFP_NGBkgFrac},
#endif // USE_NUGRAPH
                {"pfp_nhits_plane0",          kNtuple_PFP_NHitsPlane0},
        {"pfp_nhits_plane1",          kNtuple_PFP_NHitsPlane1},
        {"pfp_nhits_plane2",          kNtuple_PFP_NHitsPlane2},
        {"pfp_shw_energy_plane2",     kNtuple_PFP_ShwEnergyPlane2},
        {"pfp_shw_dedx_plane2",       kNtuple_PFP_ShwDedxPlane2},
        {"pfp_shw_open_angle",        kNtuple_PFP_ShwOpenAngle},
        {"pfp_shw_conv_gap",          kNtuple_PFP_ShwConvGap},
        {"pfp_shw_length",            kNtuple_PFP_ShwLength},
        {"pfp_trk_chi2_proton",       kNtuple_PFP_TrkChi2Proton},
        {"pfp_trk_chi2_muon",         kNtuple_PFP_TrkChi2Muon},
        {"pfp_trk_length",            kNtuple_PFP_TrkLength},
        {"pfp_trk_range_mom_proton",  kNtuple_PFP_TrkRangeMomProton},
        {"pfp_trk_range_mom_pion",    kNtuple_PFP_TrkRangeMomPion},
        {"pfp_true_pdg",              kNtuple_PFP_TruePDG},
        {"pfp_hit_purity",            kNtuple_PFP_HitPurity},
        {"pfp_hit_efficiency",      kNtuple_PFP_HitEfficiency},
        {"pfp_true_start_end_energy", kNtuple_PFP_TrueStartEndEnergy},
        {"pfp_true_dir_x",            kNtuple_PFP_TrueDirX},
        {"pfp_true_dir_y",            kNtuple_PFP_TrueDirY},
        {"pfp_true_dir_z",            kNtuple_PFP_TrueDirZ},
    };

} // namespace ana
#ifdef USE_NUGRAPH
#endif // USE_NUGRAPH


