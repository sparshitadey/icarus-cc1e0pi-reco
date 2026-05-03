#pragma once

// sbnana stuff
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// c++ stuff
#include <vector>
#include <math.h>

// root stuff
#include "TVector3.h"

// local — Vars_Base.h must come first (defines kIsPFPPionLike etc.)
#include "Config.h"
#include "CC1e0piSelection_Vars_Base.h"

namespace ana {

    // =========================================================================
    // Muon veto
    // =========================================================================
#ifdef USE_NUGRAPH
    // NuGraph: identify muon candidates via semantic category (MIP = 0)
    const Var kHaveMuonCandidate([](const caf::SRSliceProxy *slc) -> bool {
        bool haveMuonCandidate = false;
        double highestLength = 0.;

        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.start.x) || std::isnan(slc->reco.pfp[i].trk.len)) continue;

            if ((slc->reco.pfp[i].trk.len > highestLength) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 0) &&
                (slc->reco.pfp[i].trk.end.x * slc->vertex.x > 0.) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                highestLength = slc->reco.pfp[i].trk.len;
                haveMuonCandidate = true;
            }
        }
        return haveMuonCandidate;
    });

#else
    // Legacy: track score + chi2 PID
    const Var kHaveMuonCandidate([](const caf::SRSliceProxy *slc) -> bool {
        bool haveMuonCandidate = false;
        double highestLength = 0.;

        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (std::isnan(slc->reco.pfp[i].trk.start.x) || std::isnan(slc->reco.pfp[i].trk.len)) continue;

            if (slc->reco.pfp[i].trackScore < 0.5) continue;

            if ((slc->reco.pfp[i].trk.len > highestLength) &&
                (slc->reco.pfp[i].trk.chi2pid[2].chi2_muon   < MUON_CHI2_MUON_MAX) &&
                (slc->reco.pfp[i].trk.chi2pid[2].chi2_proton > MUON_CHI2_PROTON_MIN) &&
                (slc->reco.pfp[i].trk.end.x * slc->vertex.x > 0.) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                highestLength = slc->reco.pfp[i].trk.len;
                haveMuonCandidate = true;
            }
        }
        return haveMuonCandidate;
    });
#endif

    // =========================================================================
    // Leading shower (electron candidate) index
    // =========================================================================
#ifdef USE_NUGRAPH
    // NuGraph: shower-like PFPs have sem_cat == 2
    const Var kLargestRecoShowerIdx([](const caf::SRSliceProxy* slc) -> int {
        int electronIdx(-1);
        double maxNHits(-1);

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            auto const& shw = slc->reco.pfp[i].shw;
            int bestPlaneIdx = shw.plane[1].nHits > shw.plane[2].nHits ? 1 : 2;

            if ((shw.plane[bestPlaneIdx].nHits > maxNHits) &&
                (slc->reco.pfp[i].ngscore.sem_cat == 2) &&
                (slc->reco.pfp[i].parent_is_primary)) {
                electronIdx = i;
                maxNHits = shw.plane[bestPlaneIdx].nHits;
            }
        }
        return electronIdx;
    });

#else
    // Legacy: shower-like = trackScore < 0.5
    const Var kLargestRecoShowerIdx([](const caf::SRSliceProxy* slc) -> int {
        int electronIdx(-1);
        double maxNHits(-1);

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            auto const& pfp = slc->reco.pfp[i];
            auto const& shw = pfp.shw;
            int bestPlaneIdx = shw.plane[1].nHits > shw.plane[2].nHits ? 1 : 2;

            if (pfp.trackScore >= 0.5) continue; // track-like
            if (!pfp.parent_is_primary) continue;

            if (shw.plane[bestPlaneIdx].nHits > maxNHits) {
                maxNHits = shw.plane[bestPlaneIdx].nHits;
                electronIdx = i;
            }
        }
        return electronIdx;
    });
#endif

    // =========================================================================
    // Hit share — depends on reconstruction mode
    // =========================================================================
#ifdef USE_NUGRAPH
    // NuGraph: count shower hits only from shower-like PFPs (sem_cat 2/3/4)
    const Var kLargestRecoShower_BestPlaneShowerHitShare([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        int NShowerHits = 0;
        int bestPlaneIdx = 2;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if ((slc->reco.pfp[i].ngscore.sem_cat == 2) ||
                (slc->reco.pfp[i].ngscore.sem_cat == 3) ||
                (slc->reco.pfp[i].ngscore.sem_cat == 4)) {
                bestPlaneIdx = slc->reco.pfp[i].shw.plane[2].nHits > slc->reco.pfp[i].shw.plane[1].nHits ? 2 : 1;
                NShowerHits += slc->reco.pfp[i].shw.plane[bestPlaneIdx].nHits;
            }
        }
        if (NShowerHits > 0)
            return (float)slc->reco.pfp[largestShwIdx].shw.plane[bestPlaneIdx].nHits / NShowerHits;
        return -5.;
    });

#else
    // Legacy: count hits from all PFPs
    const Var kLargestRecoShower_BestPlaneShowerHitShare([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        int NShowerHits = 0;
        int bestPlaneIdx = 2;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            bestPlaneIdx = slc->reco.pfp[i].shw.plane[2].nHits > slc->reco.pfp[i].shw.plane[1].nHits ? 2 : 1;
            NShowerHits += slc->reco.pfp[i].shw.plane[bestPlaneIdx].nHits;
        }
        if (NShowerHits > 0)
            return (float)slc->reco.pfp[largestShwIdx].shw.plane[bestPlaneIdx].nHits / NShowerHits;
        return -5.;
    });
#endif

    // =========================================================================
    // Proton + other particle selection
    // =========================================================================
#ifdef USE_NUGRAPH
    // NuGraph: classify by sem_cat (0=MIP, 1=HIP, 2=shower)
    const MultiVar kNSelectedProtonsIdx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return selectedProtonIdx;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int)largestShwIdx) continue;

            if (slc->reco.pfp[i].ngscore.sem_cat == 0) {         // MIP → pion-like
                if (kIsPFPPionLike(slc, i)) NOtherParticles += 1;
            }
            else if (slc->reco.pfp[i].ngscore.sem_cat == 1) {    // HIP → proton-like
                if (kIsPFPProtonLike(slc, i)) selectedProtonIdx.push_back(i);
            }
            else if (slc->reco.pfp[i].ngscore.sem_cat == 2) {    // shower
                if (kIsPFPShowerLike(slc, i)) NOtherParticles += 1;
            }
        }
        return selectedProtonIdx;
    });

    const Var kNSelectedProtonsIdx_NOtherParticles([](const caf::SRSliceProxy* slc) -> int {
        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -1;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int)largestShwIdx) continue;

            if (slc->reco.pfp[i].ngscore.sem_cat == 0) {
                if (kIsPFPPionLike(slc, i)) NOtherParticles += 1;
            }
            else if (slc->reco.pfp[i].ngscore.sem_cat == 1) {
                if (kIsPFPProtonLike(slc, i)) selectedProtonIdx.push_back(i);
            }
            else if (slc->reco.pfp[i].ngscore.sem_cat == 2) {
                if (kIsPFPShowerLike(slc, i)) NOtherParticles += 1;
            }
        }
        return NOtherParticles;
    });

#else
    // Legacy: classify by trackScore thresholds
    const MultiVar kNSelectedProtonsIdx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
        std::vector<double> selectedProtonIdx;
        int NOtherParticles(0);
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return selectedProtonIdx;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int)largestShwIdx) continue;

            if (slc->reco.pfp[i].trackScore > 0.5) {
                if (kIsPFPPionLike(slc, i))        NOtherParticles += 1;
                else if (kIsPFPProtonLike(slc, i)) selectedProtonIdx.push_back(i);
            }
            else if (slc->reco.pfp[i].trackScore > 0.45) {
                if (kIsPFPProtonLike(slc, i))       selectedProtonIdx.push_back(i);
                else if (kIsPFPShowerLike(slc, i))  NOtherParticles += 1;
            }
            else {
                if (kIsPFPShowerLike(slc, i)) NOtherParticles += 1;
            }
        }
        return selectedProtonIdx;
    });

    const Var kNSelectedProtonsIdx_NOtherParticles([](const caf::SRSliceProxy* slc) -> int {
        std::vector<int> selectedProtonIdx;
        int NOtherParticles(0);
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -1;

        for (unsigned int i = 0; i < slc->reco.npfp; i++) {
            if (i == (unsigned int)largestShwIdx) continue;

            if (slc->reco.pfp[i].trackScore > 0.5) {
                if (kIsPFPPionLike(slc, i))        NOtherParticles += 1;
                else if (kIsPFPProtonLike(slc, i)) selectedProtonIdx.push_back(i);
            }
            else if (slc->reco.pfp[i].trackScore > 0.45) {
                if (kIsPFPProtonLike(slc, i))       selectedProtonIdx.push_back(i);
                else if (kIsPFPShowerLike(slc, i))  NOtherParticles += 1;
            }
            else {
                if (kIsPFPShowerLike(slc, i)) NOtherParticles += 1;
            }
        }
        return NOtherParticles;
    });
#endif

} // namespace ana
