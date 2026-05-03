#pragma once

// sbnana stuff
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// c++ stuff
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>

// root stuff
#include "TVector3.h"

// local
#include "Config.h"

namespace ana {

    // =========================================================================
    // Geometry helpers — use margin constants from Config.h
    // =========================================================================
    bool kIsInFV(double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z)) return false;
        return (( ( x < -61.94 - FV_MARGIN_X && x > -358.49 + FV_MARGIN_X ) ||
                  ( x >  61.94 + FV_MARGIN_X && x <  358.49 - FV_MARGIN_X )) &&
                  ( y > -181.86 + FV_MARGIN_Y && y < 134.96  - FV_MARGIN_Y ) &&
                  ( z > -894.95 + FV_MARGIN_Z_UP && z < 894.95 - FV_MARGIN_Z_DN ));
    }

    bool kIsInContained(double ex, double ey, double ez) {
        if (std::isnan(ex) || std::isnan(ey) || std::isnan(ez)) return false;
        return (( ( ex < -61.94 - CONT_MARGIN && ex > -358.49 + CONT_MARGIN ) ||
                  ( ex >  61.94 + CONT_MARGIN && ex <  358.49 - CONT_MARGIN )) &&
                  ( ey > -181.86 + CONT_MARGIN && ey < 134.96  - CONT_MARGIN ) &&
                  ( ez > -894.95 + CONT_MARGIN && ez < 894.95  - CONT_MARGIN ));
    }

    // =========================================================================
    // General event variables
    // =========================================================================
    const Var kCounting([](const caf::SRSliceProxy *slc) -> int {
        return 1;
    });

    const Var kVertex_vsTruth([](const caf::SRSliceProxy *slc) -> double {
        TVector3 VertexReco(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 VertexTrue(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
        return (VertexReco - VertexTrue).Mag();
    });

    const Var kTrue_NVisProtons([](const caf::SRSliceProxy* slc) -> int {
        int nVisProtons = 0;
        for (int ip(0); ip < slc->truth.nprim; ++ip) {
            if ((slc->truth.prim[ip].pdg == 2212) &&
                ((slc->truth.prim[ip].startE - slc->truth.prim[ip].endE) > VISIBILTY_THRESHOLD_P))
                nVisProtons += 1;
        }
        return nVisProtons;
    });

    // =========================================================================
    // Flash match variables
    // =========================================================================
    const Var kBarycenterFM_DeltaZ([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ;
    });

    const Var kBarycenterFM_DeltaZ_Trigger([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.deltaZ_Trigger;
    });

    const Var kBarycenterFM_FlashTime([](const caf::SRSliceProxy *slc) -> double {
        return slc->barycenterFM.flashTime;
    });

    // =========================================================================
    // Particle ID helpers
    // Mode-independent — Vars_Reco.h decides which PFPs to pass to these
    // =========================================================================
    bool kIsPFPPionLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.end.x)   || std::isnan(slc->reco.pfp[iPFP].trk.end.y)   || std::isnan(slc->reco.pfp[iPFP].trk.end.z))   return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_pion,
                               slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_pion,
                               slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_pion);
        double K = sqrt(pow(0.139570, 2) + pow(startMomentum.Mag(), 2));

        return (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton > PION_CHI2_MIN) &&
               kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
               ((recoStart - recoVertex).Mag() < 10) &&
               (K >= VISIBILTY_THRESHOLD_PI);
    }

    bool kIsPFPShowerLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].shw.start.x) || std::isnan(slc->reco.pfp[iPFP].shw.start.y) || std::isnan(slc->reco.pfp[iPFP].shw.start.z)) return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 recoStart(slc->reco.pfp[iPFP].shw.start.x, slc->reco.pfp[iPFP].shw.start.y, slc->reco.pfp[iPFP].shw.start.z);

        return ((recoStart - recoVertex).Mag() < 50) &&
               (slc->reco.pfp[iPFP].shw.plane[2].energy >= VISIBILTY_THRESHOLD_PI);
    }

    bool kIsPFPProtonLike(const caf::SRSliceProxy* slc, unsigned int iPFP) {
        if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.start.x) || std::isnan(slc->reco.pfp[iPFP].trk.start.y) || std::isnan(slc->reco.pfp[iPFP].trk.start.z)) return false;
        if (std::isnan(slc->reco.pfp[iPFP].trk.end.x)   || std::isnan(slc->reco.pfp[iPFP].trk.end.y)   || std::isnan(slc->reco.pfp[iPFP].trk.end.z))   return false;

        TVector3 recoVertex(slc->vertex.x, slc->vertex.y, slc->vertex.z);
        TVector3 recoStart(slc->reco.pfp[iPFP].trk.start.x, slc->reco.pfp[iPFP].trk.start.y, slc->reco.pfp[iPFP].trk.start.z);
        TVector3 startMomentum(slc->reco.pfp[iPFP].trk.dir.x * slc->reco.pfp[iPFP].trk.rangeP.p_proton,
                               slc->reco.pfp[iPFP].trk.dir.y * slc->reco.pfp[iPFP].trk.rangeP.p_proton,
                               slc->reco.pfp[iPFP].trk.dir.z * slc->reco.pfp[iPFP].trk.rangeP.p_proton);
        double K = sqrt(pow(0.9383, 2) + pow(startMomentum.Mag(), 2));

        return (slc->reco.pfp[iPFP].trk.chi2pid[2].chi2_proton <= PROTON_CHI2_MAX) &&
               kIsInContained(slc->reco.pfp[iPFP].trk.end.x, slc->reco.pfp[iPFP].trk.end.y, slc->reco.pfp[iPFP].trk.end.z) &&
               ((recoStart - recoVertex).Mag() < 10) &&
               (K >= VISIBILTY_THRESHOLD_P);
    }

    // =========================================================================
    // PlotDef struct — defined here so Vars_Derived.h can use it
    // =========================================================================
    struct PlotDef {
        std::string suffix = "";
        std::string label  = "";
        Binning bins = Binning::Simple(3, 0, 3);
        Var var = kCounting;
    };

} // namespace ana
