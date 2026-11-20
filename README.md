# 🔬 ICARUS CC1e0pi Selection Framework

> A modular CAFAna/SBNAna-based selection framework for the **BNB CC1eNp0π** channel in ICARUS, supporting both Pandora (legacy) and NuGraph2 reconstruction workflows.

---

## 📋 Table of Contents

- [Physics Signal Definition](#-physics-signal-definition)
- [What This Repo Can Produce](#-what-this-repo-can-produce)
- [Version History](#-version-history)
- [Repository Structure](#-repository-structure)
- [Environment Setup](#-environment-setup)
- [Include Order](#-include-order)
- [NuGraph vs Pandora Toggle](#-nugraph-vs-pandora-toggle)
- [Running the Code](#-running-the-code)
- [Running at Scale](#-running-at-scale)
- [Ntuple Output](#-ntuple-output)
- [Stage-1 / NuGraph2 Hit-Level Studies](#-stage-1--nugraph2-hit-level-studies)
- [Key Tunable Parameters](#-key-tunable-parameters)
- [Notes on CAF Compatibility](#-notes-on-caf-compatibility)
- [Contact](#-contact)

---

## ⚛️ Physics Signal Definition

**CC1eNp0pi** — charged-current electron neutrino interaction with:

| Particle | Requirement |
|----------|-------------|
| Primary electron | Exactly 1, KE > 200 MeV (configurable) |
| Protons | ≥ 1 visible, KE > 50 MeV, contained (configurable) |
| Pions | 0 visible (KE < 25 MeV threshold) |
| Vertex | Inside fiducial volume |

> **Note on cuts:** The initial cuts in this first version are those developed by [Riccardo Triozzi](https://github.com/rtriozzi) and collaborators for a Pandora-based NuMI νe CC selection. They serve here as a diagnostic baseline for BNB before subsequent tuning. Future versions will focus on BNB-optimised cuts.

---

## 🧭 What This Repo Can Produce

A quick map from *"I want X"* to *"run Y"*:

| You want… | Run | Lives in |
|---|---|---|
| The CC1e0pi selection ntuple (slice + per-PFP trees) | `cafe -bq` over flat CAFs | `src/CC1e0piSelection_MakeAll.C` |
| A quick efficiency sanity check on a few files | `cafe -bq` | `src/CC1e0piSelection_MakeEfficiency.C` |
| A full-sample run over thousands of CAFs | batched / parallel `cafe` | `src/run_cc1e0pi_parallel.sh`, `src/run_cc1e0pi_batches.sh` |
| Efficiency, purity, cut flow, FOM scans, Pareto front, Eν resolution | notebook on the ntuple | `notebooks/cc1e0pi_analysis.ipynb` |
| Pandora vs Pandora+NuGraph2 cut sequence straight from a flat CAF | notebook | `notebooks/caf_pandora_vs_ng2_single_file_diagnostics.ipynb` |
| How often a BNB spill has more than one neutrino (pile-up) | notebook | `notebooks/neutrino_pileup_study.ipynb` |
| A file list from grid CAF output | notebook | `notebooks/fromXMLOutput_to_FileList.ipynb` |
| NuGraph2 hit-level views, truth labels, confusion matrices | stage1 → `ICARUSNuGraphAnalyzer` → notebook | `stage1_processing/` |
| Pandora-only or NuGraph2 shower / proton ID | flip `USE_NUGRAPH` | `include/Config.h` |

---

## 🏷️ Version History

### `v1.1.0` — Scale-up + stage-1 tooling *(this release)*
- `MakeAll.C` v4: extra trees `ntuple_allsig` (every true-signal slice, no preselection) and `ntuple_truth` (one row per true signal neutrino), a `cutflow/` directory with raw counts per step, and spectra at earlier cut stages
- `NTUPLE_KEEP_PRE_FLASH` in `Config.h`: the ntuple now keeps slices that fail the flash match, so presel and FM can both be applied in Python (set to `false` for the old behaviour)
- Default mode is back to Pandora (`USE_NUGRAPH` commented out), matching this README
- Batch and parallel `cafe` runners for full-sample production
- `stage1_processing/`: FCLs, grid XMLs and notebooks for NuGraph2 hit-level studies with `ICARUSNuGraphAnalyzer`
- Example notebooks: single-file Pandora vs NG2 diagnostics, pile-up study, grid output → file list
- The v1.0 macro is tagged `macro-v3` (`git checkout macro-v3 -- src/CC1e0piSelection_MakeAll.C` to get it back)

### `v1.0.0` — Proof of Principle
- Initial framework with modular split headers
- TTree ntuple output for downstream Python analysis
- NuGraph2 toggle via `#define USE_NUGRAPH` in `Config.h`
- Supports both NuGraph2 and vanilla Pandora flat CAF inputs
- First BNB files: official CAFs produced in `v10_06_00_06p03`
- Cuts derived from Pandora-based NuMI νe CC selections as a baseline
- Per-PFP vector branches for full NuGraph vs Pandora comparison
- Efficiency, purity, and selection plots in a single `MakeAll.C` run

> Future versions will introduce BNB-optimised cuts, NuGraph-specific cut tuning, systematics handling, and cross-section machinery.

---

## 📁 Repository Structure

```
icarus-cc1e0pi-reco/
│
├── include/
│   ├── Config.h                          # Tunable constants, USE_NUGRAPH + NTUPLE_KEEP_PRE_FLASH toggles
│   ├── CC1e0piSelection_Vars_Base.h      # Geometry, basic vars, particle ID helpers
│   ├── CC1e0piSelection_Vars_Reco.h      # Mode-dependent reco vars (shower/proton)
│   ├── CC1e0piSelection_Vars_Derived.h   # Shower/proton/neutrino vars, SelectionPlots
│   ├── CC1e0piSelection_Vars_NuGraph.h   # NuGraph-only vars (compiled only if USE_NUGRAPH)
│   ├── CC1e0piSelection_Cuts.h           # All cuts, SelDef struct, cutflow steps
│   ├── CC1e0piSelection_TruthCuts.h      # Signal definition, truth cuts, interaction types
│   ├── CC1e0piSelection_Efficiency.h     # SpillMultiVar efficiency machinery + truth ntuple vars
│   └── CC1e0piSelection_Ntuple.h         # TTree branch definitions for Python analysis
│
├── src/
│   ├── CC1e0piSelection_MakeEfficiency.C # Efficiency plots + debug ntuple (fast iteration)
│   ├── CC1e0piSelection_MakeAll.C        # Full production: efficiency + purity + selection + ntuples
│   ├── run_cc1e0pi_parallel.sh           # Split a CAF list into batches, run N cafe jobs in parallel
│   └── run_cc1e0pi_batches.sh            # Same idea, one batch at a time
│
├── notebooks/                            # Analysis + example notebooks (see "What This Repo Can Produce")
├── stage1_processing/                    # NuGraph2 hit-level studies (own README inside)
│   ├── fcl/  xml/  analyzer/  notebooks/
├── CAFVersionFigures/                    # Branch-overlap figures used below
│
└── output/  batch_*/                     # ← not tracked by git
```

---

## 🛠️ Environment Setup

### 1. Start the SL7 container and set up ICARUS

```bash
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
```

### 2. Get a token

```bash
htgettoken -a htvaultprod.fnal.gov -i icarus
```

> ⚠️ **If you see `ImportError: No module named 'site'`**, your token has likely expired. Simply repeat the `htgettoken` command above.

### 3. Set up the analysis environment

Choose based on which icaruscode version you need:

#### Option A — icaruscode `v10_06_00_01p01` *(NuGraph variables)*

```bash
cd /exp/icarus/data/users/sdey2/NuGraphReco/sbnana_NG2
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source local*/setup
mrbslp
htgettoken -a htvaultprod.fnal.gov -i icarus
```

#### Option B — icaruscode `v10_06_00_04p04` *(2D deconvolution)*

Use a **fresh shell** for this, one that hasn't sourced an icaruscode dev area. Stacking two mrb areas in one shell mixes up `PATH`.

```bash
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
cd /exp/icarus/data/users/sdey2/NuGraphReco/sbnana_v04p04
source localProducts*/setup && mrbslp
htgettoken -a htvaultprod.fnal.gov -i icarus
which cafe    # should point inside sbnana_v04p04/localProducts…
```

> ⚠️ A locally built sbnana is required for v04p04 CAFs. The UPS version (`v10_01_01`) is compiled for an older kernel (`Linux64bit+3.10-2.17`) and causes heap corruption when reading v04p04 flat CAFs. The local build at `/exp/icarus/data/users/sdey2/NuGraphReco/sbnana_v04p04` compiles the same tag natively for `Linux64bit+5.14-2.17`.

> 💡 If the container start fails with a glibc/apptainer error, it can quietly wipe your bearer token. Run `htgettoken` again before reading anything from `/pnfs`.

---

## 📐 Include Order

The headers have strict dependencies — always include in this order:

```cpp
#include "Config.h"                          // no dependencies
#include "CC1e0piSelection_Vars_Base.h"       // needs Config.h
#include "CC1e0piSelection_Vars_Reco.h"       // needs Vars_Base.h
#include "CC1e0piSelection_Vars_Derived.h"    // needs Vars_Base.h + Vars_Reco.h
#include "CC1e0piSelection_Vars_NuGraph.h"    // needs Vars_Derived.h (only if USE_NUGRAPH)
#include "CC1e0piSelection_Cuts.h"            // needs all Vars headers
#include "CC1e0piSelection_TruthCuts.h"       // needs Cuts.h
#include "CC1e0piSelection_Efficiency.h"      // needs TruthCuts.h
#include "CC1e0piSelection_Ntuple.h"          // needs all of the above
```

---

## 🔀 NuGraph vs Pandora Toggle

In `Config.h`, comment or uncomment this line to switch modes:

```cpp
// #define USE_NUGRAPH   // comment out for Pandora (legacy) mode
```

| Mode | Toggle state | Shower ID | Proton ID |
|------|-------------|-----------|-----------|
| **Pandora (default)** | Commented out | Track score threshold | Chi2 PID |
| **NuGraph2** | Defined | Semantic category (`sem_cat == 2`) | Semantic category (`sem_cat == 1`) |

The `#define` controls `#ifdef USE_NUGRAPH` blocks in `Vars_Reco.h`, `Vars_NuGraph.h`, and `Ntuple.h`. Running in Pandora mode with flat CAFs that don't have NuGraph branches is fully supported — no NuGraph proxy fields are accessed.

---

## 🚀 Running the Code

Update the `listFile` path inside the `.C` file to point to your file list, then:

```bash
cd src/

# Quick test (small stats, check for errors)
cafe -bq CC1e0piSelection_MakeEfficiency.C

# Full production run (efficiency + purity + selection + ntuple)
cafe -bq CC1e0piSelection_MakeAll.C
```

**BNB file list** (official CAFs, `v10_06_00_06p03`):
```
/exp/icarus/app/users/csaia/bnb_nominalflux_caf.list
```
~28,600 files, ~5 hours as a single job. For anything this size, see [Running at Scale](#-running-at-scale) below.

**Output** goes to `output/root/`.

---

## 🏃 Running at Scale

A full sample is thousands of CAFs, so it's easier to split it into batches and run a few `cafe` jobs side by side:

```bash
cd src/
./run_cc1e0pi_parallel.sh <tag> <caf_list> [files_per_batch=50] [n_parallel=6]
```

What it does:
- splits `<caf_list>` into chunks of `files_per_batch`
- for each chunk, makes a copy of `CC1e0piSelection_MakeAll.C` with its own `listFile`, output file and function name
- runs up to `n_parallel` of them at once from `batch_<tag>/` (each with a 90-minute timeout)
- writes one ROOT file per batch to `output/root/batches_<tag>/`, renaming `.tmp` → `.root` only when a batch finishes cleanly
- refreshes your bearer token in the background while it runs

`run_cc1e0pi_batches.sh` does the same thing one batch at a time.

> 💡 **Quiet doesn't mean stuck.** Reading CAFs over xrootd can be slow, so a batch may show no output for a while. Avoid Ctrl-Z: it suspends the job rather than stopping it.

> 💡 The batch logs filter out ROOT's `Bad value` warnings so real errors are easier to spot.

> 💡 ACLiC needs to write `.so` files next to the macro. If your area is at quota, compilation fails in confusing ways, so check space first.

---

## 📊 Ntuple Output

Both macros write two TTrees under a loose preselection. By default (`NTUPLE_KEEP_PRE_FLASH = true`) that's `kNotClearCosmic && kVertexInFV`, and the flash match is left to the `pass_flash_match` / `pass_trig_flash_match` branches. Set it to `false` to apply `kFlashMatch` at fill time, as in v1.0.

| Tree | Contents |
|------|----------|
| `ntuple` | ~100 flat scalar branches per slice |
| `ntuple_pfp` | ~30 per-PFP vector branches (variable length per event) |

`MakeAll.C` also writes:

| Output | Contents |
|------|----------|
| `ntuple_allsig` | every true-signal slice, **no preselection**, with the presel / FV / FM flags and reco vertex. Useful for *"why did this signal neutrino never make it into `ntuple`?"* |
| `ntuple_truth` | one row per true CC1e0pi neutrino (E, vertex, electron KE, visible protons), the efficiency denominator |
| `cutflow/` | raw slice counts after each cumulative cut, for all slices and for signal |

> ⚠️ When comparing against an ntuple made with the v1.0 macro (`macro-v3` tag), remember that one only holds FM-passing slices. Compare like with like, e.g. efficiency after FM.

### Reading in Python

```python
import uproot
import awkward as ak
import numpy as np

f = uproot.open("output/root/CC1e0piSelection_Efficiency.root")

# flat branches → pandas DataFrame
df = f["ntuple/ntuple"].arrays(library="pd")

# per-PFP vector branches → awkward arrays
pfp = f["ntuple_pfp/ntuple_pfp"].arrays(library="ak")

# apply cuts in Python
signal_mask   = df["is_signal"] == 1
selected_mask = (df["pass_flash_match"] == 1) & (df["shower_energy"] > 0.2)

# efficiency
efficiency = (signal_mask & selected_mask).sum() / signal_mask.sum()

# threshold sweep
for threshold in [0.1, 0.15, 0.2, 0.25, 0.3]:
    sig = signal_mask & (df["true_electron_ke"] > threshold)
    eff = (sig & selected_mask).sum() / sig.sum()
    print(f"  KE > {threshold:.2f} GeV → efficiency = {eff:.3f}")
```

### Flat branch groups

| Prefix | Description |
|--------|-------------|
| `true_*` | Truth: neutrino kinematics, Q², W, particle counts, thresholds |
| `reco_vtx_*` | Reconstructed vertex position and FV flag |
| `shower_*` | Leading shower: all 3 planes of energy/dEdx/hits, NuGraph scores |
| `fm_*` | Flash match variables |
| `lead_proton_*` | Leading proton kinematics, chi2, NuGraph scores, truth matching |
| `longest_trk_*` | Muon veto inputs |
| `reco_*` | Reco neutrino energy, inelasticity, transverse momentum, e-p angle |
| `slice_*` | nu score, crumbs score, n PFPs, clear cosmic flag |
| `pass_*` | Fixed upstream cut flags (for Python-side selection) |
| `ng_*` | NuGraph event-level (only if `USE_NUGRAPH`) |
| `pfp_*` | Per-PFP: track score, NuGraph scores, hits, calorimetry, truth matching |

---

## ⚙️ Key Tunable Parameters

All thresholds live in `Config.h` — change them there and recompile:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `USE_NUGRAPH` | commented out | Pandora (off) or NuGraph2 (on) shower / proton ID |
| `NTUPLE_KEEP_PRE_FLASH` | `true` | Keep slices failing the flash match in the ntuple |
| `THRESHOLD_E` | `0.2 GeV` | Min true electron KE for signal definition |
| `VISIBILTY_THRESHOLD_P` | `0.05 GeV` | Min true proton KE |
| `VISIBILTY_THRESHOLD_PI` | `0.025 GeV` | Min true pion KE |
| `SHOWER_ENERGY_MIN` | `0.2 GeV` | Min reco shower energy cut |
| `SHOWER_DEDX_MAX` | `3.5 MeV/cm` | Max shower dEdx cut |
| `SHOWER_ANGLE_MAX` | `10 deg` | Max shower opening angle cut |
| `SHOWER_CONVGAP_MAX` | `5 cm` | Max shower conversion gap cut |
| `FM_DELTA_Z_MIN/MAX` | `0 / 100 cm` | Flash match ΔZ window |
| `FM_FLASH_T_MIN/MAX` | `-1 / 11 µs` | Flash time window |
| `PROTON_CHI2_MAX` | `100` | Max chi2 proton for proton ID |
| `MUON_CHI2_MUON_MAX` / `MUON_CHI2_PROTON_MIN` | `30` / `60` | Muon veto |
| `PION_CHI2_MIN` | `100` | Pion veto |
| `FV_MARGIN_X/Y/Z_UP/Z_DN` | `25 / 25 / 30 / 50 cm` | Fiducial volume margins |
| `CONT_MARGIN` | `5 cm` | Containment margin |

---

## 🧩 Stage-1 / NuGraph2 Hit-Level Studies

The CAF-level selection above tells you *what* gets selected. To see *why*, you need hit-level information: which hits NuGraph2 kept, what class it gave each one, and what the truth was. That lives in `stage1_processing/`:

- **`ICARUSNuGraphAnalyzer`** (truth-inclusive version) runs on stage-1 files and writes per-hit, per-slice and per-PFP trees with NuGraph2 scores and truth labels
- the FCLs and grid XMLs to run it, locally or on the grid
- notebooks for an event viewer, truth checks, recall / precision confusion matrices with event-level bootstrap errors, and an original-vs-fixed analyzer comparison

👉 Start with [`stage1_processing/README.md`](stage1_processing/README.md).

---

## 📝 Notes on CAF Compatibility

This framework was developed against flat CAFs produced in **`v10_06_00_01p01`** (the NuGraph2-enabled icaruscode version). However, flat CAF branch names and structure are highly stable across icaruscode production versions — a systematic comparison across three versions confirms this.

### v10_06_00_04p04 compatibility note
All NuGraph2 variables (`rec.slc.ng_plane.*`, `rec.slc.reco.pfp.ngscore.*`, `rec.slc.ng_filt_pass_frac`) are confirmed present in v04p04 flat CAFs. Running with `USE_NUGRAPH` defined is fully supported against v04p04 outputs, provided the local sbnana build described in Option B above is used.

### Branch overlap across icaruscode versions

A Jaccard similarity analysis of tree-qualified branch names across `v10_06_00_01p01`, `v10_06_00_04p04`, and `v10_06_00_06p03` shows **96–97.7% overlap** between versions:

![Jaccard overlap of tree-qualified branch names](CAFVersionFigures/jaccard_branch_overlap.png)

The pairwise branch count comparison confirms that the vast majority of branches are shared, with only a small number unique to each version:

![Pairwise branch overlap counts](CAFVersionFigures/pairwise_branch_counts.png)

> Analysis performed in `dev/portVersion/cafVarDiffbetweenVersions.ipynb`.

### Implications

- Variable names and structure in flat CAFs are highly stable across production versions — empirically, changes are at the **<4% level** in variable names between versions checked
- The first full BNB production run uses official CAFs from **`v10_06_00_06p03`**, and the macros compile and run correctly against those files 
- **Pandora-only analysis** (i.e. `USE_NUGRAPH` commented out) works with standard flat CAFs from any production version — no NuGraph branches are required

---

## 📬 Contact

**Sparshita Dey** — ICARUS collaboration
GitHub: [@sparshitadey](https://github.com/sparshitadey)
