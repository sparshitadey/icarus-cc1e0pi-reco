#pragma once

// Only compile this file when USE_NUGRAPH is defined
#ifdef USE_NUGRAPH

// sbnana stuff
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// c++ stuff
#include <vector>
#include <math.h>
#include <algorithm>

// local — must be included after Vars_Base.h, Vars_Reco.h, Vars_Derived.h
#include "Config.h"
#include "CC1e0piSelection_Vars_Derived.h"

namespace ana {

    // =========================================================================
    // Event-level NuGraph variables
    // =========================================================================
    const Var kNuGraph_FilterFraction([](const caf::SRSliceProxy *slc) -> double {
        if (std::isnan(slc->ng_filt_pass_frac)) return -5.;
        return slc->ng_filt_pass_frac;
    });

    // PFP counts by semantic category
    const Var kNuGraph_NShowerPFPs([](const caf::SRSliceProxy *slc) -> double {
        int kNPFPs(0);
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (!std::isnan(slc->reco.pfp[i].ngscore.sem_cat) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 2)) {
                int bestPlaneHits = slc->reco.pfp[i].shw.plane[2].nHits > slc->reco.pfp[i].shw.plane[1].nHits
                                  ? slc->reco.pfp[i].shw.plane[2].nHits
                                  : slc->reco.pfp[i].shw.plane[1].nHits;
                if (bestPlaneHits > 3) kNPFPs += 1;
            }
        }
        return kNPFPs;
    });

    const Var kNuGraph_NHIPPFPs([](const caf::SRSliceProxy *slc) -> double {
        int kNPFPs(0);
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (!std::isnan(slc->reco.pfp[i].ngscore.sem_cat) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 1)) {
                int bestPlaneHits = slc->reco.pfp[i].shw.plane[2].nHits > slc->reco.pfp[i].shw.plane[1].nHits
                                  ? slc->reco.pfp[i].shw.plane[2].nHits
                                  : slc->reco.pfp[i].shw.plane[1].nHits;
                if (bestPlaneHits > 3) kNPFPs += 1;
            }
        }
        return kNPFPs;
    });

    const Var kNuGraph_NMIPPFPs([](const caf::SRSliceProxy *slc) -> double {
        int kNPFPs(0);
        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (!std::isnan(slc->reco.pfp[i].ngscore.sem_cat) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 0)) {
                int bestPlaneHits = slc->reco.pfp[i].shw.plane[2].nHits > slc->reco.pfp[i].shw.plane[1].nHits
                                  ? slc->reco.pfp[i].shw.plane[2].nHits
                                  : slc->reco.pfp[i].shw.plane[1].nHits;
                if (bestPlaneHits > 3) kNPFPs += 1;
            }
        }
        return kNPFPs;
    });

    // =========================================================================
    // HIP tagging at vertex
    // =========================================================================
    const Var kProton_NuGraph_Ind1HIPTag([](const caf::SRSliceProxy* slc) -> int {
        return slc->ng_plane[0].ng_vtx_hip_hits;
    });

    const Var kProton_NuGraph_Ind2HIPTag([](const caf::SRSliceProxy* slc) -> int {
        return slc->ng_plane[1].ng_vtx_hip_hits;
    });

    const Var kProton_NuGraph_CollHIPTag([](const caf::SRSliceProxy* slc) -> int {
        return slc->ng_plane[2].ng_vtx_hip_hits;
    });

    const Var kProton_NuGraph_HIPTag([](const caf::SRSliceProxy* slc) -> int {
        int Ind1 = kProton_NuGraph_Ind1HIPTag(slc) >= 0 ? kProton_NuGraph_Ind1HIPTag(slc) : 0;
        int Ind2 = kProton_NuGraph_Ind2HIPTag(slc) >= 0 ? kProton_NuGraph_Ind2HIPTag(slc) : 0;
        int Coll = kProton_NuGraph_CollHIPTag(slc) >= 0 ? kProton_NuGraph_CollHIPTag(slc) : 0;
        return Ind1 + Ind2 + Coll;
    });

    const Var kProton_NuGraph_MaxHIPTag([](const caf::SRSliceProxy* slc) -> int {
        return std::max({kProton_NuGraph_Ind1HIPTag(slc),
                         kProton_NuGraph_Ind2HIPTag(slc),
                         kProton_NuGraph_CollHIPTag(slc)});
    });

    // =========================================================================
    // Shower hit counts per plane
    // =========================================================================
    const Var kNuGraph_Ind1ShowerHits([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[0].shr_hits) || slc->ng_plane[0].shr_hits < 0) return -5;
        return slc->ng_plane[0].shr_hits;
    });

    const Var kNuGraph_Ind2ShowerHits([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[1].shr_hits) || slc->ng_plane[1].shr_hits < 0) return -5;
        return slc->ng_plane[1].shr_hits;
    });

    const Var kNuGraph_CollShowerHits([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[2].shr_hits) || slc->ng_plane[2].shr_hits < 0) return -5;
        return slc->ng_plane[2].shr_hits;
    });

    const Var kNuGraph_Ind1ShowerHits_Unclustered([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[0].unclustered_shr_hits) || slc->ng_plane[0].unclustered_shr_hits < 0) return -5;
        return slc->ng_plane[0].unclustered_shr_hits;
    });

    const Var kNuGraph_Ind2ShowerHits_Unclustered([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[1].unclustered_shr_hits) || slc->ng_plane[1].unclustered_shr_hits < 0) return -5;
        return slc->ng_plane[1].unclustered_shr_hits;
    });

    const Var kNuGraph_CollShowerHits_Unclustered([](const caf::SRSliceProxy *slc) -> int {
        if (std::isnan(slc->ng_plane[2].unclustered_shr_hits) || slc->ng_plane[2].unclustered_shr_hits < 0) return -5;
        return slc->ng_plane[2].unclustered_shr_hits;
    });

    // =========================================================================
    // Leading shower NuGraph score fractions
    // All depend on kLargestRecoShowerIdx (Vars_Reco.h)
    // =========================================================================
    const Var kLargestRecoShower_NuGraph_ShowerFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].ngscore.shr_frac)) return -5.;
        return slc->reco.pfp[idx].ngscore.shr_frac;
    });

    const Var kLargestRecoShower_NuGraph_HipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].ngscore.hip_frac)) return -5.;
        return slc->reco.pfp[idx].ngscore.hip_frac;
    });

    const Var kLargestRecoShower_NuGraph_MipFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].ngscore.mip_frac)) return -5.;
        return slc->reco.pfp[idx].ngscore.mip_frac;
    });

    const Var kLargestRecoShower_NuGraph_MhlFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].ngscore.mhl_frac)) return -5.;
        return slc->reco.pfp[idx].ngscore.mhl_frac;
    });

    const Var kLargestRecoShower_NuGraph_DifFrac([](const caf::SRSliceProxy* slc) -> double {
        const int idx = kLargestRecoShowerIdx(slc);
        if (idx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[idx].ngscore.dif_frac)) return -5.;
        return slc->reco.pfp[idx].ngscore.dif_frac;
    });

    // =========================================================================
    // Leading proton NuGraph score fractions
    // =========================================================================
    const Var kLeadingProton_NuGraph_HipFrac([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        std::vector<double> protonMomenta, hipFracs;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            protonMomenta.push_back(p.Mag());
            hipFracs.push_back(slc->reco.pfp[i].ngscore.hip_frac);
        }
        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](double i1, double i2){ return protonMomenta[i1] > protonMomenta[i2]; });
        return hipFracs[idx[0]];
    });

    // =========================================================================
    // Append NuGraph plots to SelectionPlots and SelectionPlots_LowStat
    // This runs at include time — the vectors are already defined in Vars_Derived.h
    // =========================================================================
    struct NuGraphPlotsAppender {
        NuGraphPlotsAppender() {
            // append to SelectionPlots
            std::vector<PlotDef> extra = {
                {"ngfiltfrac",       "NuGraph2 S/N",              Binning::Simple(40, 0,    1),  kNuGraph_FilterFraction},
                {"ngtaggedshws",     "NG2-tagged showers [#]",    Binning::Simple(8,  0,    8),  kNuGraph_NShowerPFPs},
                {"ngshwfrac",        "NuGraph2 shw_frac",         Binning::Simple(40, 0,    1),  kLargestRecoShower_NuGraph_ShowerFrac},
                {"nghipfrac",        "NuGraph2 hip_frac",         Binning::Simple(40, 0,    1),  kLargestRecoShower_NuGraph_HipFrac},
                {"ngmipfrac",        "NuGraph2 mip_frac",         Binning::Simple(40, 0,    1),  kLargestRecoShower_NuGraph_MipFrac},
                {"ngmhlfrac",        "NuGraph2 mhl_frac",         Binning::Simple(40, 0,    1),  kLargestRecoShower_NuGraph_MhlFrac},
                {"ngdiffrac",        "NuGraph2 dif_frac",         Binning::Simple(40, 0,    1),  kLargestRecoShower_NuGraph_DifFrac},
            };
            for (auto& p : extra) SelectionPlots.push_back(p);

            // append to SelectionPlots_LowStat
            std::vector<PlotDef> extraLS = {
                {"ngshwfrac",        "NG2 e^{#pm} shw_frac",      Binning::Simple(20, 0,    1),  kLargestRecoShower_NuGraph_ShowerFrac},
                {"nghipfrac",        "NG2 e^{#pm} hip_frac",      Binning::Simple(20, 0,    1),  kLargestRecoShower_NuGraph_HipFrac},
                {"ngmipfrac",        "NG2 e^{#pm} mip_frac",      Binning::Simple(20, 0,    1),  kLargestRecoShower_NuGraph_MipFrac},
                {"ngmhlfrac",        "NG2 e^{#pm} mhl_frac",      Binning::Simple(20, 0,    1),  kLargestRecoShower_NuGraph_MhlFrac},
                {"ngdiffrac",        "NG2 e^{#pm} dif_frac",      Binning::Simple(20, 0,    1),  kLargestRecoShower_NuGraph_DifFrac},
                {"nghipfracleadp",   "NG2 p_{1} hip_frac",        Binning::Simple(20, 0,    1),  kLeadingProton_NuGraph_HipFrac},
                {"nghiptagind1",     "Ind-1 ng_vtx_hip_hits [#]", Binning::Simple(20, 0,   40),  kProton_NuGraph_Ind1HIPTag},
                {"nghiptagind2",     "Ind-2 ng_vtx_hip_hits [#]", Binning::Simple(20, 0,   40),  kProton_NuGraph_Ind2HIPTag},
                {"nghiptagcoll",     "Coll ng_vtx_hip_hits [#]",  Binning::Simple(20, 0,   40),  kProton_NuGraph_CollHIPTag},
                {"nghiptag",         "All ng_vtx_hip_hits [#]",   Binning::Simple(20, 0,   40),  kProton_NuGraph_HIPTag},
                {"nghiptagmax",      "Max ng_vtx_hip_hits [#]",   Binning::Simple(20, 0,   40),  kProton_NuGraph_MaxHIPTag},
                {"shrhitsind1",      "Ind1 shr_hits [#]",         Binning::Simple(20, 0, 1200),  kNuGraph_Ind1ShowerHits},
                {"shrhitsind2",      "Ind2 shr_hits [#]",         Binning::Simple(20, 0, 1200),  kNuGraph_Ind2ShowerHits},
                {"shrhitscoll",      "Coll shr_hits [#]",         Binning::Simple(20, 0, 1200),  kNuGraph_CollShowerHits},
                {"unclshrhitsind1",  "Ind1 uncl. shr_hits [#]",  Binning::Simple(15, 0,  300),  kNuGraph_Ind1ShowerHits_Unclustered},
                {"unclshrhitsind2",  "Ind2 uncl. shr_hits [#]",  Binning::Simple(15, 0,  300),  kNuGraph_Ind2ShowerHits_Unclustered},
                {"unclshrhitscoll",  "Coll uncl. shr_hits [#]",  Binning::Simple(15, 0,  300),  kNuGraph_CollShowerHits_Unclustered},
            };
            for (auto& p : extraLS) SelectionPlots_LowStat.push_back(p);
        }
    };

    // instantiating this appends the NuGraph plots at include time
    static NuGraphPlotsAppender sNuGraphPlotsAppender;

} // namespace ana

#endif // USE_NUGRAPH
