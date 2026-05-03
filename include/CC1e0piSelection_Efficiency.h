#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <vector>
#include <math.h>

// local — must be included after all Vars, Cuts, TruthCuts headers
#include "Config.h"

namespace ana {

    // =========================================================================
    // True neutrino energy of all signal neutrinos in the spill
    // Denominator for efficiency calculations
    // =========================================================================
    const SpillMultiVar kCC1e0p1Signal_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> trueEnergies;
        for (auto const& nu : sr->mc.nu) {
            if (kIsTrueCC1e0pi(nu))
                trueEnergies.push_back(nu.E);
        }
        return trueEnergies;
    });

    // =========================================================================
    // True electron KE of all signal neutrinos in the spill
    // =========================================================================
    const SpillMultiVar kCC1e0p1Signal_TrueElectronEnergy([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> trueEnergies;
        for (auto const& nu : sr->mc.nu) {
            if (kIsTrueCC1e0pi(nu)) {
                for (int ip(0); ip < nu.nprim; ++ip) {
                    if (abs(nu.prim[ip].pdg) == 11) {
                        trueEnergies.push_back(nu.prim[ip].startE - nu.prim[ip].endE);
                        break;
                    }
                }
            }
        }
        return trueEnergies;
    });

    // =========================================================================
    // No pile-up variants
    // =========================================================================
    const SpillMultiVar kCC1e0p1Signal_NoPileup_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> trueEnergies;
        for (auto const& nu : sr->mc.nu) {
            if (kIsTrueCC1e0pi(nu) && (sr->mc.nnu == 1))
                trueEnergies.push_back(nu.E);
        }
        return trueEnergies;
    });

    // =========================================================================
    // Generic CC nu in FV (for cross-checks)
    // =========================================================================
    const SpillMultiVar kNuCC_TrueNeutrinoEnergy([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> trueEnergies;
        for (auto const& nu : sr->mc.nu) {
            if (nu.iscc &&
                kIsInFV(nu.position.x, nu.position.y, nu.position.z) &&
                sr->mc.nnu > 0)
                trueEnergies.push_back(nu.E);
        }
        return trueEnergies;
    });

    // =========================================================================
    // Factory: true neutrino energy for signal neutrinos that were
    // reconstructed and pass a given selection cut step
    //
    // Logic:
    //   - loop over true signal neutrinos in the spill
    //   - for each, look for a reconstructed slice that:
    //       * matches the true neutrino index (pile-up safe)
    //       * is truth-tagged as signal (kTrueCC1e0pi)
    //       * passes the reco cut step
    //   - count each true neutrino at most once (first matching slice wins)
    // =========================================================================
    SpillMultiVar kCC1e0p1Signal_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (!kIsTrueCC1e0pi(nu)) continue;

                bool counted = false;
                for (auto const& islc : sr->slc) {
                    if (counted) break;
                    if (islc.truth.index == nu.index &&
                        kTrueCC1e0pi(&islc) &&
                        cut(&islc)) {
                        selectedEnergies.push_back(islc.truth.E);
                        counted = true;
                    }
                }
            }

            return selectedEnergies;
        });
    }

    // =========================================================================
    // Factory: true electron KE for selected signal neutrinos
    // =========================================================================
    SpillMultiVar kCC1e0p1Signal_TrueElectronEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (!kIsTrueCC1e0pi(nu)) continue;

                bool counted = false;
                for (auto const& islc : sr->slc) {
                    if (counted) break;
                    if (islc.truth.index == nu.index &&
                        kTrueCC1e0pi(&islc) &&
                        cut(&islc)) {
                        for (int ip(0); ip < nu.nprim; ++ip) {
                            if (abs(nu.prim[ip].pdg) == 11) {
                                selectedEnergies.push_back(nu.prim[ip].startE - nu.prim[ip].endE);
                                counted = true;
                                break;
                            }
                        }
                    }
                }
            }

            return selectedEnergies;
        });
    }

    // =========================================================================
    // Factory: true neutrino energy for generic CC nu in FV
    // =========================================================================
    SpillMultiVar kNuCC_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (!nu.iscc) continue;
                if (!kIsInFV(nu.position.x, nu.position.y, nu.position.z)) continue;
                if (sr->mc.nnu == 0) continue;

                bool counted = false;
                for (auto const& islc : sr->slc) {
                    if (counted) break;
                    if (islc.truth.index == nu.index &&
                        kIsCC(&islc) &&
                        kIsNuinFV(&islc) &&
                        cut(&islc)) {
                        selectedEnergies.push_back(islc.truth.E);
                        counted = true;
                    }
                }
            }

            return selectedEnergies;
        });
    }

    // =========================================================================
    // Factory: no pile-up variant
    // =========================================================================
    SpillMultiVar kCC1e0p1Signal_NoPileup_TrueNeutrinoEnergy_MakeSelectionStep(const Cut& cut)
    {
        return SpillMultiVar([cut](const caf::SRSpillProxy* sr) -> std::vector<double> {
            std::vector<double> selectedEnergies;

            for (auto const& nu : sr->mc.nu) {
                if (!kIsTrueCC1e0pi(nu)) continue;
                if (sr->mc.nnu != 1) continue;

                bool counted = false;
                for (auto const& islc : sr->slc) {
                    if (counted) break;
                    if (islc.truth.index == nu.index &&
                        kTrueCC1e0pi(&islc) &&
                        cut(&islc)) {
                        selectedEnergies.push_back(islc.truth.E);
                        counted = true;
                    }
                }
            }

            return selectedEnergies;
        });
    }

} // namespace ana
