// Reference copy: builds only inside icaruscode.
// Source: github.com/sparshitadey/icaruscode, branch feature/sdey_NuGraph2_TruthAnalyzer_v04p04, commit 617e12ef

/**
 * @file    icaruscode/TPC/NuGraph/ICARUSNuGraphAnalyzer_module.cc
 * @brief   Analyzer to check on TPC hits and their relationship with NuGraph2 predictions.
 * @author  Riccardo Triozzi ( triozzi@pd.infn.it ), Giuseppe Cerati ( cerati@fnal.gov )
 * @date    September 9, 2025
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TTree.h"
#include "art_root_io/TFileService.h"

#include <limits>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <memory>

#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

// --- truth-matching additions (ported from ICARUSHDF5Maker_module.cc) ---
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "larcore/CoreUtils/ServiceUtil.h" // for lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// --- overlay truth-matching additions ---
// Live BackTrackerService::ChannelToTrackIDEs returned zero matches on a
// BNB+overlay sample despite SimChannel/MCParticle labels all confirmed
// correct against the file (see session notes) -- overlay's timing/channel
// mixing appears to break live backtracking's assumptions. This production
// already ships a precomputed Assns<Hit, MCParticle, BackTrackerHitMatchingData>
// as a workaround for exactly this; HitTruthAssnsLabel lets the analyzer use
// it instead of live backtracking when set.
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

class ICARUSNuGraphAnalyzer;

using std::vector;

// --- true-semantic-label additions ---
// Port of pynuml's SimpleLabels/StandardLabels (github.com/nugraph/nugraph,
// pynuml/pynuml/labels/{standard,simple}.py) -- the scheme NuGraph's own
// training truth is generated with. SimpleLabels collapses StandardLabels'
// {pion,muon}->MIP and {kaon,hadron}->HIP, keeping shower/michel/diffuse
// and an "invisible" 6th class. Per pynuml's process/hitgraph.py, a hit's
// truth_semantic is -1 ONLY when its dominant g4_id matched no particle at
// all (fillna(-1) after a left-merge on g4_id -- a true orphan/noise hit);
// a hit matched to a particle labeled "invisible" keeps that label's own
// index (5) rather than being folded into -1. hitgraph.py's metadata then
// drops "invisible" from the 5 trained semantic_classes, i.e. NuGraph is
// never trained to predict it -- but the truth value itself stays 5, not -1.
//
// KNOWN DEVIATIONS from the Python source, worth validating before trusting:
//  - Python uses the `particle` package for pdgid.charge()/is_baryon()/
//    is_nucleus(). Here, "neutral" is a hardcoded list of the PDG codes that
//    actually appear in G4/LArTPC output (gamma, neutral mesons/hyperons,
//    neutrons, neutrinos) rather than a general charge-database lookup --
//    ROOT's TDatabasePDG would do this more generally, but pulling in
//    ROOT::Physics for one boolean check wasn't worth the link-dependency
//    (it also turned out not to resolve cleanly for this MODULE-type
//    plugin on this build, for reasons not worth chasing further). Baryon/
//    nucleus detection uses the standard PDG Monte Carlo numbering scheme
//    directly, same as before, unaffected by this change.
//  - gamma_threshold (0.02 GeV) and hadron_threshold (0.2 GeV) are
//    StandardLabels' constructor defaults. If ICARUS training overrode
//    these, update kGammaThreshold/kHadronThreshold below to match.
//  - Any PDG/process combination the Python raises an Exception on is
//    logged via mf::LogWarning and bucketed as "invisible" here instead of
//    aborting the job.
namespace {

enum SemanticLabel { kMIP = 0, kHIP = 1, kShower = 2, kMichel = 3, kDiffuse = 4, kInvisible = 5 };

constexpr double kGammaThreshold  = 0.02; // GeV/c, momentum cut for gamma/conversion electrons
constexpr double kHadronThreshold = 0.2;  // GeV/c, proton HIP threshold
constexpr double kProtonDiffuseCut = 0.0015; // GeV/c, hIoni knock-on proton -> diffuse

// Neutral PDG codes actually reachable from G4 output in a LAr TPC. Not a
// general charge lookup -- see deviation note above. Only matters for the
// CoupledTransportation check just below; any baryon/nucleus that shows up
// here (neutron, Lambda, Xi0, ...) gets re-classified by PdgIsBaryonOrNucleus
// regardless of what this function returns, so it's not load-bearing there.
bool PdgIsNeutral(int pdg)
{
  int a = std::abs(pdg);
  switch (a) {
    case 22:   // gamma
    case 111:  // pi0
    case 113:  // rho0
    case 130:  // K0L
    case 310:  // K0S
    case 311:  // K0
    case 2112: // neutron
    case 3122: // Lambda
    case 3212: // Sigma0
    case 3322: // Xi0
    case 12: case 14: case 16: // neutrinos (rarely tracked as G4 particles, but harmless to include)
      return true;
    default:
      return false;
  }
}

// Standard PDG Monte Carlo numbering: ground-state baryons are 4-digit
// codes q1 q2 q3 (2J+1) with 1<=qi<=6; ions/nuclei are 10-digit codes.
bool PdgIsBaryonOrNucleus(int pdg)
{
  int a = std::abs(pdg);
  if (a >= 1000000000) return true; // nucleus
  if (a < 1000 || a > 9999) return false;
  int d4 = (a / 1000) % 10, d3 = (a / 100) % 10, d2 = (a / 10) % 10;
  return (d4 >= 1 && d4 <= 6) && (d3 >= 1 && d3 <= 6) && (d2 >= 1 && d2 <= 6);
}

// Mirrors pynuml's s(part, particles) -> (sl, slc). slc < 0 means "no
// forced propagation" (Python's None): children compute their own label.
// slc >= 0 means children inherit that label directly (sticky subtree).
struct LabelResult { int sl; int slc; };

LabelResult ComputeParticleLabel(int pdg, std::string const& sp, std::string const& ep,
                                  double mom, int parentPdg)
{
  int apdg = std::abs(pdg);
  int sl = -1, slc = -1;

  // Neutral particle that exited the volume without interacting. NOTE: in
  // the Python source this does NOT return early -- it falls through to the
  // baryon/nucleus/tau overrides below, same as the particle_processor
  // dispatch does. A neutral baryon (e.g. a neutron) hitting this branch
  // gets re-classified to diffuse by that override, same as an interacting
  // one would. An earlier version of this function returned immediately
  // here, which silently misclassified exiting neutrons/neutral hyperons as
  // "invisible" instead of "diffuse" -- fixed.
  if (PdgIsNeutral(pdg) && ep == "CoupledTransportation") {
    sl = kInvisible; // slc stays unset/-1 (Python's None) in this branch too
  }
  else {
    switch (apdg) {
      case 211: case 221: case 331: case 223: // pion-like -> MIP
        sl = kMIP; break;
      case 13: // muon -> MIP
        sl = kMIP; break;
      case 321: case 411: // kaon, D meson -> HIP
        sl = kHIP; break;
      case 111: case 311: case 310: case 130: case 113: // neutral mesons
        sl = kInvisible; break;
      case 11: // electron/positron
        if (sp == "primary") { sl = kShower; slc = kShower; }
        else if (std::abs(parentPdg) == 13 &&
                 (sp == "muMinusCaptureAtRest" || sp == "muPlusCaptureAtRest" || sp == "Decay")) {
          sl = kMichel; slc = kMichel;
        }
        else if (sp == "conv" || ep == "conv" || sp == "compt" || ep == "compt") {
          if (mom >= kGammaThreshold) { sl = kShower; slc = kShower; }
          else { sl = kDiffuse; slc = kDiffuse; }
        }
        else if (sp == "muIoni" || sp == "hIoni" || sp == "eIoni") {
          if (sp == "muIoni") sl = kMIP;
          else if (sp == "hIoni") {
            if (std::abs(parentPdg) == 2212) {
              sl = kHIP;
              if (mom <= kProtonDiffuseCut) sl = kDiffuse;
            }
            else sl = kMIP;
          }
          else sl = kDiffuse; // eIoni
        }
        else if (sp == "eBrem" || ep == "phot" || ep == "photonNuclear" || ep == "eIoni") {
          sl = kDiffuse;
        }
        else if (ep == "StepLimiter" || ep == "annihil" || ep == "eBrem" ||
                 sp == "hBertiniCaptureAtRest" || ep == "FastScintillation" ||
                 sp == "muPairProd" || sp == "phot") {
          sl = kDiffuse; slc = kDiffuse;
        }
        else {
          mf::LogWarning("ICARUSNuGraphAnalyzer")
            << "unhandled electron process combo: start=" << sp << " end=" << ep
            << " -- bucketing as diffuse";
          sl = kDiffuse;
        }
        break;
      case 22: // gamma
        if (sp == "conv" || ep == "conv" || sp == "compt" || ep == "compt") {
          if (mom >= kGammaThreshold) { sl = kShower; slc = kShower; }
          else { sl = kDiffuse; slc = kDiffuse; }
        }
        else if (sp == "eBrem" || ep == "phot" || ep == "photonNuclear") {
          sl = kDiffuse;
        }
        else {
          mf::LogWarning("ICARUSNuGraphAnalyzer")
            << "unhandled gamma process combo: start=" << sp << " end=" << ep
            << " -- bucketing as diffuse";
          sl = kDiffuse;
        }
        break;
      default:
        sl = -1; // resolved below by baryon/nucleus/tau checks, or falls through to "unrecognized"
    }
  }

  if (PdgIsBaryonOrNucleus(pdg)) {
    sl = (apdg == 2212 && mom >= kHadronThreshold) ? kHIP : kDiffuse;
    slc = -1;
  }
  if (apdg == 15) { sl = kHIP; slc = -1; } // charged tau, treated as highly-ionising

  if (sl < 0) {
    mf::LogWarning("ICARUSNuGraphAnalyzer")
      << "unrecognized particle for semantic labeling: pdg=" << pdg
      << " start=" << sp << " end=" << ep << " -- bucketing as invisible";
    sl = kInvisible;
  }
  return {sl, slc};
}

} // namespace

// --- overlay truth-matching additions ---
// Common representation for "a particle contributing to this hit's charge,
// with what fraction" -- fed by EITHER live ChannelToTrackIDEs (sim::TrackIDE
// has the same two fields) OR the precomputed BackTrackerHitMatchingData
// Assns, so the dominant-contributor/edep-tree logic below doesn't need to
// know or care which source produced it.
struct Contributor { int trackID; float energyFrac; };

class ICARUSNuGraphAnalyzer : public art::EDAnalyzer {
public:
  explicit ICARUSNuGraphAnalyzer(fhicl::ParameterSet const& p);

  ICARUSNuGraphAnalyzer(ICARUSNuGraphAnalyzer const&) = delete;
  ICARUSNuGraphAnalyzer(ICARUSNuGraphAnalyzer&&) = delete;
  ICARUSNuGraphAnalyzer& operator=(ICARUSNuGraphAnalyzer const&) = delete;
  ICARUSNuGraphAnalyzer& operator=(ICARUSNuGraphAnalyzer&&) = delete;

  void analyze(art::Event const& e) override;

private:
  TTree *_treeHit, *_treeEvt, *_treeEdep;
  int _run, _subrun, _event, _id, _wire, _plane, _tpc, _cryo;
  float _x_filter, _MIP, _HIP, _shower, _michel, _diffuse, _time;
  int _islc, _icluster, _ipfp, _has_ng_output;
  float _vtx_x, _vtx_y, _vtx_z;
  std::string fNGLabel, fSliceLabel, fPandoraLabel;

  // --- truth-matching additions ---
  art::InputTag fTruthLabel;
  int _truth_g4id, _truth_pdg, _truth_isnu;
  float _truth_energyFrac;
  int _evt_is_cc, _evt_nu_pdg;
  float _evt_nu_energy, _evt_lep_energy;
  int _edep_hitid, _edep_g4id, _edep_pdg, _edep_isnu;
  float _edep_energyFrac;
  int _truth_semantic;
  std::unordered_map<int, int> fSemanticByG4Id; // filled once per event by BuildSemanticLabels()

  void BuildSemanticLabels(art::Event const& e);

  // --- overlay truth-matching additions ---
  // Empty (default): use live BackTrackerService::ChannelToTrackIDEs, same
  // as always. Non-empty: use this Assns<Hit,MCParticle,BackTrackerHitMatchingData>
  // product instead (e.g. "mcassociationsGausCryoE") -- needed for samples
  // where live backtracking doesn't work, e.g. BNB+overlay.
  art::InputTag fHitTruthAssnsLabel;
  std::vector<Contributor> GetContributors(
      art::Ptr<recob::Hit> const& hit, size_t hitIdx,
      cheat::BackTrackerService const* bt, detinfo::DetectorClocksData const& clockData,
      art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>* hitTruthAssns) const;

  // --- slice-level truth additions ---
  TTree* _treeSlc;
  int _slc_islc, _slc_nhits, _slc_nnuhits, _slc_istruenuslice;
  float _slc_purity, _slc_completeness;

  // --- PFParticle-level additions ---
  TTree* _treePfp;
  int _pfp_islc, _pfp_id, _pfp_pdg, _pfp_isprimary, _pfp_nhits, _pfp_nhits_scored;
  float _pfp_mean_MIP, _pfp_mean_HIP, _pfp_mean_shower, _pfp_mean_michel, _pfp_mean_diffuse;
  int _pfp_truth_pdg, _pfp_truth_isnu_nhits, _pfp_nhits_truthmatched;
  float _pfp_truth_electron_frac, _pfp_purity;
};

ICARUSNuGraphAnalyzer::ICARUSNuGraphAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fNGLabel{p.get<std::string>("NuGraphLabel", "NGMultiSliceCryoE")},
  fSliceLabel{p.get<std::string>("SliceLabel", "NCCSlicesCryoE")},
  fPandoraLabel{p.get<std::string>("PandoraLabel", "pandoraGausCryoE")},
  fTruthLabel{p.get<art::InputTag>("TruthLabel", "generator")},
  fHitTruthAssnsLabel{p.get<art::InputTag>("HitTruthAssnsLabel", art::InputTag())}
{
  art::ServiceHandle<art::TFileService> tfs;
  _treeHit = tfs->make<TTree>("NuGraphHitOutput", "NuGraphHitOutput");
  _treeHit->Branch("run", &_run, "run/I");
  _treeHit->Branch("subrun", &_subrun, "subrun/I");
  _treeHit->Branch("event", &_event, "event/I");
  _treeHit->Branch("id", &_id, "id/I");
  _treeHit->Branch("wire", &_wire, "wire/I");
  _treeHit->Branch("plane", &_plane, "plane/I");
  _treeHit->Branch("tpc", &_tpc, "tpc/I");
  _treeHit->Branch("cryo", &_cryo, "cryo/I");
  _treeHit->Branch("x_filter", &_x_filter, "x_filter/F");
  _treeHit->Branch("MIP", &_MIP, "MIP/F");
  _treeHit->Branch("HIP", &_HIP, "HIP/F");
  _treeHit->Branch("shower", &_shower, "shower/F");
  _treeHit->Branch("michel", &_michel, "michel/F");
  _treeHit->Branch("diffuse", &_diffuse, "diffuse/F");
  _treeHit->Branch("time", &_time, "time/F");
  _treeHit->Branch("islc", &_islc, "islc/I");
  // Cluster's key (position) in the fPandoraLabel Cluster collection -- like
  // islc/ipfp, only meaningful within a given (run,subrun,event), not a
  // globally stable ID. -1 if the hit has no associated cluster.
  _treeHit->Branch("icluster", &_icluster, "icluster/I");
  _treeHit->Branch("ipfp", &_ipfp, "ipfp/I");
  // 1 if this hit had valid NuGraph filter+semantic output, 0 otherwise.
  // Filter on this before trusting x_filter/MIP/HIP/shower/michel/diffuse --
  // those five are set to -1 (not left stale) when this is 0.
  _treeHit->Branch("has_ng_output", &_has_ng_output, "has_ng_output/I");
  // --- truth-matching additions ---
  _treeHit->Branch("truth_g4id", &_truth_g4id, "truth_g4id/I");
  _treeHit->Branch("truth_pdg", &_truth_pdg, "truth_pdg/I");
  _treeHit->Branch("truth_energyFrac", &_truth_energyFrac, "truth_energyFrac/F");
  _treeHit->Branch("truth_isnu", &_truth_isnu, "truth_isnu/I");
  // -1=no g4_id match (orphan/noise hit), 0=MIP, 1=HIP, 2=shower, 3=Michel,
  // 4=diffuse, 5=invisible (matched a particle, but one SimpleLabels itself
  // excludes from the 5 trained classes -- distinct from -1, don't conflate
  // them), -2=g4_id matched but the particle-hierarchy walk never reached it
  // (see comment in analyze() -- watch for this one, it signals a broken
  // parent chain rather than ordinary noise).
  _treeHit->Branch("truth_semantic", &_truth_semantic, "truth_semantic/I");

  _treeEvt = tfs->make<TTree>("NuGraphEventOutput", "NuGraphEventOutput");
  _treeEvt->Branch("run", &_run, "run/I");
  _treeEvt->Branch("subrun", &_subrun, "subrun/I");
  _treeEvt->Branch("event", &_event, "event/I");
  _treeEvt->Branch("vtx_x", &_vtx_x, "vtx_x/F");
  _treeEvt->Branch("vtx_y", &_vtx_y, "vtx_y/F");
  _treeEvt->Branch("vtx_z", &_vtx_z, "vtx_z/F");
  // --- truth-matching additions ---
  _treeEvt->Branch("is_cc", &_evt_is_cc, "is_cc/I");
  _treeEvt->Branch("nu_pdg", &_evt_nu_pdg, "nu_pdg/I");
  _treeEvt->Branch("nu_energy", &_evt_nu_energy, "nu_energy/F");
  _treeEvt->Branch("lep_energy", &_evt_lep_energy, "lep_energy/F");

  // --- truth-matching additions ---
  // One row per (hit, contributing g4 particle) -- full energy-fraction
  // breakdown, mirrors ICARUSHDF5Maker's energyDepNtuple. Use this if you need
  // more than the single dominant contributor stored in NuGraphHitOutput.
  _treeEdep = tfs->make<TTree>("NuGraphEdepOutput", "NuGraphEdepOutput");
  _treeEdep->Branch("run", &_run, "run/I");
  _treeEdep->Branch("subrun", &_subrun, "subrun/I");
  _treeEdep->Branch("event", &_event, "event/I");
  _treeEdep->Branch("hit_id", &_edep_hitid, "hit_id/I");
  _treeEdep->Branch("g4id", &_edep_g4id, "g4id/I");
  _treeEdep->Branch("pdg", &_edep_pdg, "pdg/I");
  _treeEdep->Branch("energyFrac", &_edep_energyFrac, "energyFrac/F");
  _treeEdep->Branch("isnu", &_edep_isnu, "isnu/I");

  // --- slice-level truth additions ---
  // One row per (event, slice). "True nu slice" flag and purity/completeness
  // all use the ANY-contributor-from-nu criterion (hitIsNuAny in analyze()),
  // matching ICARUSTrueNuSliceHitsProducer_module.cc's own slice-selection
  // logic exactly (ChannelToTrackIDEs + NeutrinoSet(), break on first nu
  // contributor found). This is deliberately NOT the same criterion as
  // NuGraphHitOutput's truth_isnu, which reflects only the dominant
  // contributor -- the two can disagree on individual hits in busy/cosmic
  // events, so don't expect n_nu_hits here to equal a count of truth_isnu==1
  // over the same hits.
  _treeSlc = tfs->make<TTree>("NuGraphSliceOutput", "NuGraphSliceOutput");
  _treeSlc->Branch("run", &_run, "run/I");
  _treeSlc->Branch("subrun", &_subrun, "subrun/I");
  _treeSlc->Branch("event", &_event, "event/I");
  _treeSlc->Branch("islc", &_slc_islc, "islc/I");
  _treeSlc->Branch("n_hits", &_slc_nhits, "n_hits/I");
  _treeSlc->Branch("n_nu_hits", &_slc_nnuhits, "n_nu_hits/I");
  _treeSlc->Branch("purity", &_slc_purity, "purity/F");
  _treeSlc->Branch("completeness", &_slc_completeness, "completeness/F");
  _treeSlc->Branch("is_true_nu_slice", &_slc_istruenuslice, "is_true_nu_slice/I");

  // --- PFParticle-level additions ---
  // One row per (event, PFParticle). Hits are associated via the standard
  // Pandora PFParticle->Cluster->Hit chain (see the PFP association block in
  // analyze()). NG score columns are means over this PFP's hits that had
  // valid NuGraph output (n_hits_scored may be < n_hits if some hits failed
  // the same NG filter/semantic-validity check used for NuGraphHitOutput).
  // truth_pdg is a hit-count majority vote over each hit's dominant-
  // contributor truth_pdg (NuGraphHitOutput convention, not the any-
  // contributor one used for slice selection).
  _treePfp = tfs->make<TTree>("NuGraphPFPOutput", "NuGraphPFPOutput");
  _treePfp->Branch("run", &_run, "run/I");
  _treePfp->Branch("subrun", &_subrun, "subrun/I");
  _treePfp->Branch("event", &_event, "event/I");
  _treePfp->Branch("islc", &_pfp_islc, "islc/I");
  _treePfp->Branch("ipfp", &_pfp_id, "ipfp/I");
  _treePfp->Branch("pandora_pdg", &_pfp_pdg, "pandora_pdg/I");
  _treePfp->Branch("is_primary", &_pfp_isprimary, "is_primary/I");
  _treePfp->Branch("n_hits", &_pfp_nhits, "n_hits/I");
  _treePfp->Branch("n_hits_scored", &_pfp_nhits_scored, "n_hits_scored/I");
  _treePfp->Branch("mean_MIP", &_pfp_mean_MIP, "mean_MIP/F");
  _treePfp->Branch("mean_HIP", &_pfp_mean_HIP, "mean_HIP/F");
  _treePfp->Branch("mean_shower", &_pfp_mean_shower, "mean_shower/F");
  _treePfp->Branch("mean_michel", &_pfp_mean_michel, "mean_michel/F");
  _treePfp->Branch("mean_diffuse", &_pfp_mean_diffuse, "mean_diffuse/F");
  _treePfp->Branch("truth_pdg", &_pfp_truth_pdg, "truth_pdg/I");
  _treePfp->Branch("truth_isnu_nhits", &_pfp_truth_isnu_nhits, "truth_isnu_nhits/I");
  // Denominator for `purity` below -- hits with no truth match at all are
  // excluded from the majority vote, so this is <= n_hits.
  _treePfp->Branch("n_hits_truthmatched", &_pfp_nhits_truthmatched, "n_hits_truthmatched/I");
  _treePfp->Branch("truth_electron_frac", &_pfp_truth_electron_frac, "truth_electron_frac/F");
  _treePfp->Branch("purity", &_pfp_purity, "purity/F");
}

void ICARUSNuGraphAnalyzer::analyze(art::Event const& e)
{
  // --- truth-matching additions ---
  auto const* bt = art::ServiceHandle<cheat::BackTrackerService>().get();
  auto const* pi = lar::providerFrom<cheat::ParticleInventoryService>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  // --- true-semantic-label additions ---
  BuildSemanticLabels(e);

  // Event-level neutrino truth (mirrors ICARUSHDF5Maker's eventNtupleNu block).
  // Skipped (not filled) if there isn't exactly one MCTruth, e.g. pile-up samples.
  auto truthHandle = e.getHandle<std::vector<simb::MCTruth>>(fTruthLabel);
  if (truthHandle.isValid() && truthHandle->size() == 1) {
    simb::MCNeutrino const& nutruth = truthHandle->at(0).GetNeutrino();
    _evt_is_cc     = (nutruth.CCNC() == simb::kCC);
    _evt_nu_pdg    = nutruth.Nu().PdgCode();
    _evt_nu_energy = nutruth.Nu().E();
    _evt_lep_energy = nutruth.Lepton().E();
    _vtx_x = nutruth.Nu().Vx();
    _vtx_y = nutruth.Nu().Vy();
    _vtx_z = nutruth.Nu().Vz();
  }
  else {
    _evt_is_cc = _evt_nu_pdg = -1;
    _evt_nu_energy = _evt_lep_energy = -1.;
    _vtx_x = _vtx_y = _vtx_z = std::numeric_limits<float>::quiet_NaN();
  }
  _event  = e.event();
  _subrun = e.subRun();
  _run    = e.run();
  _treeEvt->Fill();

  // Get slices.
  const std::vector<art::Ptr<recob::Slice>> slices = e.getProduct<vector<art::Ptr<recob::Slice>>>(fSliceLabel);
  art::FindManyP<recob::Hit> findMHitsFromSlice(slices, e, fPandoraLabel);

  // Map hits to slices.
  std::vector<art::Ptr<recob::Hit>> allHits;
  std::map<unsigned int, int> hitToSliceID;
  for (size_t islc = 0; islc < slices.size(); ++islc) {
    art::Ptr<recob::Slice> slice = slices[islc];
    const std::vector<art::Ptr<recob::Hit>> hitsInSlice = findMHitsFromSlice.at(islc);
    for (auto const& h : hitsInSlice) {
      hitToSliceID[h.key()] = slice.key();
      allHits.emplace_back(h);
    }
  }

  art::InputTag fNGFilterLabel{fNGLabel, "filter"};
  art::InputTag fNGSemanticLabel{fNGLabel, "semantic"};

  art::FindOneP<anab::FeatureVector<1>> find1FilterFromHit(allHits, e, fNGFilterLabel);
  art::FindOneP<anab::FeatureVector<5>> find1SemanticFromHit(allHits, e, fNGSemanticLabel);

  // --- overlay truth-matching additions ---
  // Built once per event, against allHits, same as the NG FindOneP's above.
  // Null (not built) when fHitTruthAssnsLabel is unset -> GetContributors()
  // falls back to live BackTrackerService, i.e. no change from before.
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hitTruthAssns;
  if (!fHitTruthAssnsLabel.label().empty()) {
    art::Handle<art::Assns<recob::Hit, simb::MCParticle, anab::BackTrackerHitMatchingData>> assnsHandle;
    if (e.getByLabel(fHitTruthAssnsLabel, assnsHandle) && assnsHandle.isValid()) {
      hitTruthAssns = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
          allHits, e, fHitTruthAssnsLabel);
    }
    else {
      mf::LogWarning("ICARUSNuGraphAnalyzer")
        << "HitTruthAssnsLabel '" << fHitTruthAssnsLabel.encode() << "' configured but not found "
        << "in this event -- falling back to live BackTrackerService for this event only";
    }
  }

  // --- PFParticle-level additions ---
  // Standard Pandora PFParticle->Cluster->Hit chain (no direct PFP->Hit
  // association exists in the data model; Cluster is the hop in between).
  // hitKeyToPfpIdx maps a hit's art::Ptr key to its PFParticle's position in
  // pfpVec (not PFParticle::Self() -- that's used as the tree's "ipfp" value
  // once we're filling NuGraphPFPOutput below). hitKeyToClusterKey is the
  // same walk's free byproduct -- fills the previously-dead `icluster`
  // branch on NuGraphHitOutput (cluster.key(), i.e. its position in the
  // fPandoraLabel Cluster collection -- not a globally stable ID, same
  // caveat as islc/ipfp: only meaningful within a given event).
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  std::unordered_map<unsigned int, int> hitKeyToPfpIdx;
  std::unordered_map<unsigned int, int> hitKeyToClusterKey;
  {
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    if (e.getByLabel(fPandoraLabel, pfpHandle) && pfpHandle.isValid()) {
      art::fill_ptr_vector(pfpVec, pfpHandle);
      art::FindManyP<recob::Cluster> pfpClusterAssns(pfpHandle, e, fPandoraLabel);
      auto clusterHandle = e.getHandle<std::vector<recob::Cluster>>(fPandoraLabel);
      if (clusterHandle.isValid()) {
        art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, e, fPandoraLabel);
        for (size_t ip = 0; ip < pfpVec.size(); ++ip) {
          for (auto const& cl : pfpClusterAssns.at(ip)) {
            for (auto const& h : clusterHitAssns.at(cl.key())) {
              hitKeyToPfpIdx[h.key()] = ip;
              hitKeyToClusterKey[h.key()] = int(cl.key());
            }
          }
        }
      }
      else {
        mf::LogWarning("ICARUSNuGraphAnalyzer")
          << "no Cluster collection at " << fPandoraLabel
          << " -- PFP hit association (ipfp, icluster, NuGraphPFPOutput) will be empty this event";
      }
    }
  }

  struct PfpAccum {
    int nhits = 0, nhitsScored = 0, nhitsIsnu = 0, nhitsElectron = 0;
    double sumMIP = 0, sumHIP = 0, sumShower = 0, sumMichel = 0, sumDiffuse = 0;
    std::unordered_map<int, int> truthPdgVotes;
    int islc = -1;
  };
  std::unordered_map<int, PfpAccum> pfpAccum; // keyed by index into pfpVec

  // --- slice-level truth additions ---
  std::unordered_map<int, int> sliceNHits, sliceNNuHits;

  for (size_t hitIdx = 0; hitIdx < allHits.size(); hitIdx++) {
    art::Ptr<recob::Hit> hit = allHits[hitIdx];

    // Fill event information.
    _event  = e.event();
    _subrun = e.subRun();
    _run    = e.run();
    _id     = hit.key();

    // Fill slice ID.
    auto itSlice = hitToSliceID.find(hit.key());
    _islc = (itSlice != hitToSliceID.end()) ? itSlice->second : -1;

    // --- truth-matching additions ---
    // Find the dominant true particle contributing to this hit (highest
    // energy fraction). Contributors come from GetContributors() -- either
    // live ChannelToTrackIDEs (default) or the precomputed Assns path when
    // fHitTruthAssnsLabel is set (see overlay truth-matching additions).
    _truth_g4id = -1;
    _truth_pdg = 0;
    _truth_energyFrac = -1.;
    _truth_isnu = -1;
    bool hitIsNuAny = false; // matches ICARUSTrueNuSliceHitsProducer: ANY contributor from nu, not just dominant
    std::vector<Contributor> const contributors =
      GetContributors(hit, hitIdx, bt, clockData, hitTruthAssns.get());
    for (auto const& c : contributors) {
      int tid = c.trackID;
      if (pi->TrackIdToParticle_P(abs(tid))) tid = abs(tid);

      // Full breakdown -> edep tree (one row per contributor).
      _edep_hitid = hit.key();
      _edep_g4id = tid;
      _edep_energyFrac = c.energyFrac;
      _edep_pdg = 0;
      _edep_isnu = -1;
      if (auto const* p = pi->TrackIdToParticle_P(abs(tid))) {
        _edep_pdg = p->PdgCode();
        auto mct = pi->TrackIdToMCTruth_P(abs(tid));
        _edep_isnu = (mct.isAvailable() ? mct->NeutrinoSet() : 0);
      }
      _treeEdep->Fill();
      if (_edep_isnu == 1) hitIsNuAny = true;

      // Track dominant contributor -> summary fields on the hit tree.
      if (c.energyFrac > _truth_energyFrac) {
        _truth_g4id = tid;
        _truth_energyFrac = c.energyFrac;
        _truth_pdg = _edep_pdg;
        _truth_isnu = _edep_isnu;
      }
    }

    // --- true-semantic-label additions ---
    // Matches pynuml's hitgraph.py convention: -1 means the hit's g4_id
    // never matched any particle row at all (true orphan/noise hit --
    // `hits['semantic_label'].fillna(-1)` after a left-merge). A hit whose
    // dominant particle *did* match but was labeled "invisible" keeps that
    // label's literal index (5), same as training truth -- NOT folded into
    // -1. A THIRD case -- g4_id present but BuildSemanticLabels' walk never
    // reached it (broken parent chain, e.g. an ancestor missing from
    // ParticleList) -- is what hitgraph.py calls an "orphaned hit" and
    // discards the WHOLE EVENT for; we can't discard the event from inside
    // an analyzer, so it's flagged as -2 here instead of silently folding
    // into -1. If -2 shows up at any real rate, that's a sign of a genuine
    // particle-hierarchy gap worth investigating, not ordinary noise.
    _truth_semantic = -1;
    if (_truth_g4id >= 0) {
      auto slit = fSemanticByG4Id.find(_truth_g4id);
      _truth_semantic = (slit != fSemanticByG4Id.end()) ? slit->second : -2;
    }

    // --- PFParticle-level additions ---
    auto pfpIt = hitKeyToPfpIdx.find(hit.key());
    int ipfpIdx = (pfpIt != hitKeyToPfpIdx.end()) ? pfpIt->second : -1;
    _ipfp = (ipfpIdx >= 0) ? int(pfpVec[ipfpIdx]->Self()) : -1;
    auto clIt = hitKeyToClusterKey.find(hit.key());
    _icluster = (clIt != hitKeyToClusterKey.end()) ? clIt->second : -1;
    if (ipfpIdx >= 0) {
      PfpAccum& acc = pfpAccum[ipfpIdx];
      acc.nhits++;
      if (acc.islc < 0) acc.islc = _islc;
      // Only vote when the hit actually has a truth match -- otherwise
      // untruth-matched hits (_truth_pdg defaults to 0) would count as a
      // "vote" for PDG 0 as if that were a real particle type, inflating
      // purity for PFPs with lots of unmatched hits.
      if (_truth_g4id >= 0) acc.truthPdgVotes[_truth_pdg]++;
      if (_truth_isnu == 1) acc.nhitsIsnu++;
      if (_truth_pdg == 11) acc.nhitsElectron++;
    }

    // --- slice-level truth additions ---
    // Uses hitIsNuAny (any contributor from nu), matching
    // ICARUSTrueNuSliceHitsProducer's slice-selection criterion exactly --
    // NOT _truth_isnu (which reflects only the dominant contributor and is
    // what NuGraphHitOutput reports per hit).
    sliceNHits[_islc]++;
    if (hitIsNuAny) sliceNNuHits[_islc]++;

    const art::Ptr<anab::FeatureVector<1>> filter = find1FilterFromHit.at(hitIdx);
    const art::Ptr<anab::FeatureVector<5>> semantic = find1SemanticFromHit.at(hitIdx);

    if (filter.isNull() || semantic.isNull()) {
      std::cout << "Hit #" << hitIdx << " with key " << hit.key() << " had no valid NuGraphOutput. Skipping." << "\n";

      // Explicitly reset -- these are class members, so without this a row
      // for a hit with no NG output would silently carry the PREVIOUS hit's
      // real scores. Filter on has_ng_output==1 downstream, don't trust
      // these five columns being nonzero as a proxy for validity.
      _has_ng_output = 0;
      _x_filter = _MIP = _HIP = _shower = _michel = _diffuse = -1.;

      // Fill hit description.
      _wire  = hit->WireID().Wire;
      _plane = hit->WireID().Plane;
      _tpc   = hit->WireID().TPC;
      _cryo  = hit->WireID().Cryostat;
      _time  = hit->PeakTime();

      _treeHit->Fill();
      continue;
    }
    else {
      _has_ng_output = 1;
      // Fill NuGraph predictions.
      _x_filter = filter->at(0);
      _MIP      = semantic->at(0);
      _HIP      = semantic->at(1);
      _shower   = semantic->at(2);
      _michel   = semantic->at(3);
      _diffuse  = semantic->at(4);

      // --- PFParticle-level additions ---
      if (ipfpIdx >= 0) {
        PfpAccum& acc = pfpAccum[ipfpIdx];
        acc.nhitsScored++;
        acc.sumMIP += _MIP; acc.sumHIP += _HIP; acc.sumShower += _shower;
        acc.sumMichel += _michel; acc.sumDiffuse += _diffuse;
      }

      // Fill hit description.
      _wire  = hit->WireID().Wire;
      _plane = hit->WireID().Plane;
      _tpc   = hit->WireID().TPC;
      _cryo  = hit->WireID().Cryostat;
      _time  = hit->PeakTime();

      _treeHit->Fill();
    }
  }

  // --- slice-level truth additions ---
  int trueNuSlice = -1, maxNuHits = 0;
  for (auto const& kv : sliceNNuHits) {
    if (kv.second > maxNuHits) { maxNuHits = kv.second; trueNuSlice = kv.first; }
  }
  int totalNuHits = 0;
  for (auto const& kv : sliceNNuHits) totalNuHits += kv.second;

  for (auto const& kv : sliceNHits) {
    _slc_islc = kv.first;
    _slc_nhits = kv.second;
    _slc_nnuhits = sliceNNuHits.count(kv.first) ? sliceNNuHits.at(kv.first) : 0;
    _slc_purity = (_slc_nhits > 0) ? float(_slc_nnuhits) / _slc_nhits : -1.;
    _slc_completeness = (totalNuHits > 0) ? float(_slc_nnuhits) / totalNuHits : -1.;
    _slc_istruenuslice = (_slc_islc == trueNuSlice) ? 1 : 0;
    _treeSlc->Fill();
  }

  // --- PFParticle-level additions ---
  for (auto const& kv : pfpAccum) {
    int ipfpIdx = kv.first;
    PfpAccum const& acc = kv.second;
    art::Ptr<recob::PFParticle> const& pfp = pfpVec[ipfpIdx];

    _pfp_islc = acc.islc;
    _pfp_id = int(pfp->Self());
    _pfp_pdg = pfp->PdgCode();
    _pfp_isprimary = pfp->IsPrimary() ? 1 : 0;
    _pfp_nhits = acc.nhits;
    _pfp_nhits_scored = acc.nhitsScored;
    _pfp_mean_MIP     = (acc.nhitsScored > 0) ? float(acc.sumMIP / acc.nhitsScored) : -1.;
    _pfp_mean_HIP     = (acc.nhitsScored > 0) ? float(acc.sumHIP / acc.nhitsScored) : -1.;
    _pfp_mean_shower  = (acc.nhitsScored > 0) ? float(acc.sumShower / acc.nhitsScored) : -1.;
    _pfp_mean_michel  = (acc.nhitsScored > 0) ? float(acc.sumMichel / acc.nhitsScored) : -1.;
    _pfp_mean_diffuse = (acc.nhitsScored > 0) ? float(acc.sumDiffuse / acc.nhitsScored) : -1.;
    _pfp_truth_isnu_nhits = acc.nhitsIsnu;
    _pfp_truth_electron_frac = (acc.nhits > 0) ? float(acc.nhitsElectron) / acc.nhits : -1.;

    // Majority vote (by hit count) over dominant-contributor truth_pdg,
    // among hits that actually had a truth match.
    _pfp_truth_pdg = 0;
    int bestVotes = -1, totalVotes = 0;
    for (auto const& v : acc.truthPdgVotes) {
      totalVotes += v.second;
      if (v.second > bestVotes) { bestVotes = v.second; _pfp_truth_pdg = v.first; }
    }
    _pfp_nhits_truthmatched = totalVotes;
    _pfp_purity = (totalVotes > 0) ? float(bestVotes) / totalVotes : -1.;

    _treePfp->Fill();
  }
}

// --- true-semantic-label additions ---
// Recursive top-down walk over the G4 particle hierarchy, mirroring pynuml's
// walk(part, particles, depth, sl, il) (semantic half only -- instance
// labeling isn't needed for the confusion matrix). Populates fSemanticByG4Id
// for every G4 track ID in the event; hits look up their dominant truth
// contributor's entry (see the truth_g4id computation in analyze()).
void ICARUSNuGraphAnalyzer::BuildSemanticLabels(art::Event const&)
{
  fSemanticByG4Id.clear();
  auto const* pi = lar::providerFrom<cheat::ParticleInventoryService>();
  sim::ParticleList const& plist = pi->ParticleList();

  std::unordered_map<int, simb::MCParticle const*> partByTid;
  std::unordered_map<int, std::vector<int>> children;
  for (auto const& pr : plist) {
    simb::MCParticle const* p = pr.second;
    if (!p) continue;
    partByTid[p->TrackId()] = p;
    children[p->Mother()].push_back(p->TrackId());
  }

  std::function<void(int, int)> walk = [&](int trackid, int forced) {
    auto pit = partByTid.find(trackid);
    if (pit == partByTid.end()) return;
    simb::MCParticle const* part = pit->second;

    int parentPdg = 0;
    auto mit = partByTid.find(part->Mother());
    if (mit != partByTid.end()) parentPdg = mit->second->PdgCode();

    int sl, slc;
    if (forced >= 0) {
      // Sticky propagation: this subtree was already tagged (e.g. inside an
      // EM shower), so every descendant inherits the same label verbatim.
      sl = forced;
      slc = forced;
    }
    else {
      LabelResult lr = ComputeParticleLabel(part->PdgCode(), part->Process(), part->EndProcess(),
                                             part->P(), parentPdg);
      sl = lr.sl;
      slc = lr.slc;
    }

    fSemanticByG4Id[trackid] = sl;

    auto cit = children.find(trackid);
    if (cit != children.end())
      for (int child : cit->second) walk(child, slc);
  };

  for (auto const& pr : plist) {
    simb::MCParticle const* p = pr.second;
    if (p && p->Mother() == 0) walk(p->TrackId(), -1);
  }
}

// --- overlay truth-matching additions ---
// hitIdx must be this hit's position in the SAME allHits vector both
// find1FilterFromHit/find1SemanticFromHit and hitTruthAssns were built
// against -- callers must not reorder allHits between construction and use.
std::vector<Contributor> ICARUSNuGraphAnalyzer::GetContributors(
    art::Ptr<recob::Hit> const& hit, size_t hitIdx,
    cheat::BackTrackerService const* bt, detinfo::DetectorClocksData const& clockData,
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>* hitTruthAssns) const
{
  std::vector<Contributor> contributors;

  if (hitTruthAssns) {
    // Precomputed path (e.g. overlay samples where live backtracking doesn't
    // work): read the production's own Assns instead of calling
    // ChannelToTrackIDEs. ideFraction plays the same role as TrackIDE's
    // energyFrac -- fraction of this hit's charge attributed to this
    // particle. NOTE: "ideFraction" is the standard LArSoft field name for
    // this (as used the same way in AnalysisTree-style modules elsewhere in
    // the stack); verify against BackTrackerMatchingData.h at build time if
    // this doesn't compile -- that's the one part of this addition not
    // directly confirmed against this specific file's contents this session.
    std::vector<art::Ptr<simb::MCParticle>> parts;
    std::vector<anab::BackTrackerHitMatchingData const*> matchData;
    hitTruthAssns->get(hitIdx, parts, matchData);
    contributors.reserve(parts.size());
    for (size_t k = 0; k < parts.size(); ++k) {
      contributors.push_back({ parts[k]->TrackId(), matchData[k]->ideFraction });
    }
  }
  else {
    // Original, default path: live backtracking. Unchanged from before this
    // addition -- MPVMPR/pure-MC behavior is identical.
    std::vector<sim::TrackIDE> const& ides =
      bt->ChannelToTrackIDEs(clockData, hit->Channel(), hit->StartTick(), hit->EndTick());
    contributors.reserve(ides.size());
    for (auto const& tide : ides) {
      contributors.push_back({ tide.trackID, tide.energyFrac });
    }
  }
  return contributors;
}

DEFINE_ART_MODULE(ICARUSNuGraphAnalyzer)
