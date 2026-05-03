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

// local — must be included after Vars_Base.h and Vars_Reco.h
#include "Config.h"

namespace ana {

    // =========================================================================
    // Shower property vars
    // All depend on kLargestRecoShowerIdx (defined in Vars_Reco.h)
    // =========================================================================
    const Var kLargestRecoShower_CollEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5.;
        return slc->reco.pfp[largestShwIdx].shw.plane[2].energy;
    });

    const Var kLargestRecoShower_TruePdg([](const caf::SRSliceProxy* slc) -> int {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.pdg)) return -5;
        return slc->reco.pfp[largestShwIdx].shw.truth.p.pdg;
    });

    const Var kLargestRecoShower_ColldEdx([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx)) return -5.;
        return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx;
    });

    const Var kLargestRecoShower_AvailabledEdx([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        if (!std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx) &&
             slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx > 0)
            return slc->reco.pfp[largestShwIdx].shw.plane[2].dEdx;
        if (!std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[1].dEdx) &&
             slc->reco.pfp[largestShwIdx].shw.plane[1].dEdx > 0)
            return slc->reco.pfp[largestShwIdx].shw.plane[1].dEdx;
        if (!std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[0].dEdx) &&
             slc->reco.pfp[largestShwIdx].shw.plane[0].dEdx > 0)
            return slc->reco.pfp[largestShwIdx].shw.plane[0].dEdx;
        return -5.;
    });

    const Var kLargestRecoShower_TrackScore([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[largestShwIdx].trackScore)) return -5.;
        return slc->reco.pfp[largestShwIdx].trackScore;
    });

    const Var kLargestRecoShower_OpenAngle([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        return 180. * slc->reco.pfp[largestShwIdx].shw.open_angle / M_PI;
    });

    const Var kLargestRecoShower_ConvGap([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        return slc->reco.pfp[largestShwIdx].shw.conversion_gap;
    });

    const Var kLargestRecoShower_TrueEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.startE) ||
            std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.endE)) return -5.;
        return slc->reco.pfp[largestShwIdx].shw.truth.p.startE -
               slc->reco.pfp[largestShwIdx].shw.truth.p.endE;
    });

    const Var kLargestRecoShower_CollEnergy_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        double recoEnergy = kLargestRecoShower_CollEnergy(slc);
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.startE) ||
            std::isnan(slc->reco.pfp[largestShwIdx].shw.truth.p.endE)) return -5.;
        double trueEnergy = slc->reco.pfp[largestShwIdx].shw.truth.p.startE -
                            slc->reco.pfp[largestShwIdx].shw.truth.p.endE;
        return (recoEnergy - trueEnergy) / trueEnergy;
    });

    // =========================================================================
    // Proton property vars
    // Depend on kNSelectedProtonsIdx (defined in Vars_Reco.h)
    // =========================================================================
    const Var kNSelectedProtons_N([](const caf::SRSliceProxy* slc) -> double {
        return kNSelectedProtonsIdx(slc).size();
    });

    const Var kLeadingProtonMomentum([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        std::vector<double> protonMomenta;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            protonMomenta.push_back(p.Mag());
        }
        std::sort(protonMomenta.begin(), protonMomenta.end(), std::greater<>());
        return protonMomenta[0];
    });

    const Var kSubLeadingProtonMomentum([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty() || selectedProtonIdx.size() < 2) return -5.;
        std::vector<double> protonMomenta;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            protonMomenta.push_back(p.Mag());
        }
        std::sort(protonMomenta.begin(), protonMomenta.end(), std::greater<>());
        return protonMomenta[1];
    });

    const Var kLeadingProton_TrackScore([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        std::vector<double> protonMomenta, trackScores;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            protonMomenta.push_back(p.Mag());
            trackScores.push_back(slc->reco.pfp[i].trackScore);
        }
        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](double i1, double i2){ return protonMomenta[i1] > protonMomenta[i2]; });
        return trackScores[idx[0]];
    });

    const Var kSubLeadingProton_TrackScore([](const caf::SRSliceProxy* slc) -> double {
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty() || selectedProtonIdx.size() < 2) return -5.;
        std::vector<double> protonMomenta, trackScores;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            protonMomenta.push_back(p.Mag());
            trackScores.push_back(slc->reco.pfp[i].trackScore);
        }
        std::vector<unsigned int> idx(protonMomenta.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](double i1, double i2){ return protonMomenta[i1] > protonMomenta[i2]; });
        return trackScores[idx[1]];
    });

    // =========================================================================
    // Neutrino reconstructed energy variables
    // =========================================================================
    const Var kRecoNeutrino_CC0piEnergy([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5.;
        double E_e = slc->reco.pfp[largestShwIdx].shw.plane[2].energy;
        double E_p = 0.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            E_p += sqrt(pow(0.9383, 2) + pow(p.Mag(), 2)) - 0.9383;
        }
        return E_e + E_p;
    });

    const Var kRecoNeutrino_CC0piEnergy_VsTruth([](const caf::SRSliceProxy* slc) -> double {
        if (kLargestRecoShowerIdx(slc) == -1) return -5.;
        if (kNSelectedProtonsIdx(slc).empty()) return -5.;
        if (std::isnan(slc->truth.E)) return -5.;
        return (kRecoNeutrino_CC0piEnergy(slc) - slc->truth.E) / slc->truth.E;
    });

    const Var kRecoNeutrino_CC0piHadronEnergy([](const caf::SRSliceProxy* slc) -> double {
        if (kLargestRecoShowerIdx(slc) == -1) return -5.;
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        double E_p = 0.;
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            TVector3 p(slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                       slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
            E_p += sqrt(pow(0.9383, 2) + pow(p.Mag(), 2)) - 0.9383;
        }
        return E_p;
    });

    const Var kRecoNeutrino_CC0piInelasticity([](const caf::SRSliceProxy* slc) -> double {
        if (kLargestRecoShowerIdx(slc) == -1) return -5.;
        if (kNSelectedProtonsIdx(slc).empty()) return -5.;
        return kRecoNeutrino_CC0piHadronEnergy(slc) / kRecoNeutrino_CC0piEnergy(slc);
    });

    const Var kRecoNeutrino_CC0piTransverseMomentum([](const caf::SRSliceProxy* slc) -> double {
        const int largestShwIdx = kLargestRecoShowerIdx(slc);
        if (largestShwIdx == -1) return -5.;
        std::vector<double> selectedProtonIdx = kNSelectedProtonsIdx(slc);
        if (selectedProtonIdx.empty()) return -5.;
        if (std::isnan(slc->reco.pfp[largestShwIdx].shw.plane[2].energy)) return -5.;
        double P_e = sqrt(pow(slc->reco.pfp[largestShwIdx].shw.plane[2].energy, 2) - pow(0.510998e-3, 2));
        TVector3 pE(P_e * slc->reco.pfp[largestShwIdx].shw.dir.x,
                    P_e * slc->reco.pfp[largestShwIdx].shw.dir.y,
                    P_e * slc->reco.pfp[largestShwIdx].shw.dir.z);
        TVector3 pP(0., 0., 0.);
        for (auto i : selectedProtonIdx) {
            if (std::isnan(slc->reco.pfp[i].trk.rangeP.p_proton)) return -5.;
            pP.SetXYZ(pP.X() + slc->reco.pfp[i].trk.dir.x * slc->reco.pfp[i].trk.rangeP.p_proton,
                      pP.Y() + slc->reco.pfp[i].trk.dir.y * slc->reco.pfp[i].trk.rangeP.p_proton,
                      pP.Z() + slc->reco.pfp[i].trk.dir.z * slc->reco.pfp[i].trk.rangeP.p_proton);
        }
        TVector3 pTot = pE + pP;
        return sqrt(pow(pTot.X(), 2) + pow(pTot.Y(), 2));
    });

    // =========================================================================
    // SelectionPlots and SelectionPlots_LowStat
    // These are the mode-independent core plots — NuGraph-specific plots
    // are appended by Vars_NuGraph.h when USE_NUGRAPH is defined
    // =========================================================================
    std::vector<PlotDef> SelectionPlots = {
        {"count",                "Counts [#]",                                         Binning::Simple(3,  0,    3),  kCounting},

        // electron variables
        {"collenergy",           "E_{Coll} [GeV]",                                     Binning::Simple(30, 0,    3),  kLargestRecoShower_CollEnergy},
        {"colldedx",             "dE/dx_{Coll} [MeV/cm]",                              Binning::Simple(40, 0,    9),  kLargestRecoShower_ColldEdx},
        {"availdedx",            "dE/dx_{Coll, Ind} [MeV/cm]",                         Binning::Simple(40, 0,    9),  kLargestRecoShower_AvailabledEdx},
        {"trackscore",           "Track score",                                        Binning::Simple(50, 0,    1),  kLargestRecoShower_TrackScore},
        {"openangle",            "Opening angle [deg.]",                               Binning::Simple(40, 0,   30),  kLargestRecoShower_OpenAngle},
        {"convgap",              "Conversion gap [cm]",                                Binning::Simple(40, 0,   10),  kLargestRecoShower_ConvGap},
        {"hitshare",             "Hit share",                                          Binning::Simple(40, 0,    1),  kLargestRecoShower_BestPlaneShowerHitShare},

        // slice-level variables
        {"muonrej",              "#mu veto",                                           Binning::Simple(2,  0,    2),  kHaveMuonCandidate},

        // proton variables
        {"leadproton",           "P_{p_{1}} [GeV/c]",                                  Binning::Simple(30, 0,    2),  kLeadingProtonMomentum},
        {"subleadproton",        "P_{p_{2}} [GeV/c]",                                  Binning::Simple(30, 0,    2),  kSubLeadingProtonMomentum},
        {"leadpts",              "Track score (p_{1})",                                Binning::Simple(50, 0,    1),  kLeadingProton_TrackScore},
        {"subleadpts",           "Track score (p_{2})",                                Binning::Simple(50, 0,    1),  kSubLeadingProton_TrackScore},

        // neutrino variables
        {"reconuenergy",         "E^{reco}_{#nu} [GeV]",                               Binning::Simple(30, 0,    3),  kRecoNeutrino_CC0piEnergy},
        {"nuenergyres",          "(E^{reco}_{#nu}-E^{true}_{#nu})/E^{true}_{#nu}",     Binning::Simple(40, -1, 0.5),  kRecoNeutrino_CC0piEnergy_VsTruth},
        {"collenergyres",        "(E^{reco}_{e}-E^{true}_{e})/E^{true}_{e}",           Binning::Simple(40, -1, 0.5),  kLargestRecoShower_CollEnergy_VsTruth},
        {"inelasticity",         "y = E^{reco}_{had.}/E^{reco}_{#nu}",                 Binning::Simple(40, 0,    1),  kRecoNeutrino_CC0piInelasticity},
        {"tranvmomentum",        "P_{T} [GeV/c]",                                      Binning::Simple(30, 0,    3),  kRecoNeutrino_CC0piTransverseMomentum},

        // light information
        {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",               Binning::Simple(40, 0,  150),  kBarycenterFM_DeltaZ_Trigger},
        {"barycenterfmdeltaz",   "Barycenter-FM #DeltaZ [cm]",                         Binning::Simple(15, 0,  150),  kBarycenterFM_DeltaZ},
        {"barycenterfmtime",     "Barycenter-FM time [#mus]",                          Binning::Simple(40, -1,  14),  kBarycenterFM_FlashTime},
    };

    std::vector<PlotDef> SelectionPlots_LowStat = {
        {"count",                "Counts [#]",                                         Binning::Simple(3,  0,    3),  kCounting},

        // electron variables
        {"collenergy",           "E_{Coll} [GeV]",                                     Binning::Simple(25, 0,  2.5),  kLargestRecoShower_CollEnergy},
        {"colldedx",             "dE/dx_{Coll} [MeV/cm]",                              Binning::Simple(40, 0,   11),  kLargestRecoShower_ColldEdx},
        {"availdedx",            "dE/dx_{Coll, Ind} [MeV/cm]",                         Binning::Simple(40, 0,   11),  kLargestRecoShower_AvailabledEdx},
        {"trackscore",           "Track score",                                        Binning::Simple(25, 0.2, 0.8), kLargestRecoShower_TrackScore},
        {"openangle",            "Opening angle [deg.]",                               Binning::Simple(25, 0,   20),  kLargestRecoShower_OpenAngle},
        {"convgap",              "Conversion gap [cm]",                                Binning::Simple(25, 0,    8),  kLargestRecoShower_ConvGap},
        {"hitshare",             "Hit share",                                          Binning::Simple(20, 0,    1),  kLargestRecoShower_BestPlaneShowerHitShare},

        // proton variables
        {"leadproton",           "P_{p_{1}} [GeV/c]",                                  Binning::Simple(15, 0,  1.5),  kLeadingProtonMomentum},
        {"subleadproton",        "P_{p_{2}} [GeV/c]",                                  Binning::Simple(15, 0,  1.5),  kSubLeadingProtonMomentum},

        // neutrino variables
        {"reconuenergy",         "E^{reco}_{#nu} [GeV]",                               Binning::Simple(20, 0,    3),  kRecoNeutrino_CC0piEnergy},
        {"collenergyres",        "(E^{reco}_{e}-E^{true}_{e})/E^{true}_{e}",           Binning::Simple(40, -1, 0.5),  kLargestRecoShower_CollEnergy_VsTruth},
        {"inelasticity",         "y = E^{reco}_{had.}/E^{reco}_{#nu}",                 Binning::Simple(20, 0,    1),  kRecoNeutrino_CC0piInelasticity},
        {"tranvmomentum",        "P_{T} [GeV/c]",                                      Binning::Simple(20, 0,    2),  kRecoNeutrino_CC0piTransverseMomentum},

        // light information
        {"barycenterfmdeltaztr", "Barycenter-FM #DeltaZ (trigger) [cm]",               Binning::Simple(15, 0,  150),  kBarycenterFM_DeltaZ_Trigger},
        {"barycenterfmdeltaz",   "Barycenter-FM #DeltaZ [cm]",                         Binning::Simple(15, 0,  150),  kBarycenterFM_DeltaZ},
        {"barycenterfmtime",     "Barycenter-FM time [#mus]",                          Binning::Simple(40, -4,  14),  kBarycenterFM_FlashTime},
    };

} // namespace ana
