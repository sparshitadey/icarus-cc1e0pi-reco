# 🧩 Stage-1 Processing: NuGraph2 Hit-Level Studies

The CAF-level selection in the main repo tells you *what* gets selected. This folder is for the step before that: running **`ICARUSNuGraphAnalyzer`** on stage-1 files to see, hit by hit, what NuGraph2 kept, what class it gave each hit, and what the truth actually was.

---

## 🗺️ The chain

```
stage0 ──► stage1 (NuGraph2 filter + Pandora) ──┬──► caf ──► CAF-level selection (main repo)
                                                └──► ICARUSNuGraphAnalyzer ──► hist_*.root ──► notebooks/ (this folder)
```

In the NuGraph2 stage-1 FCL (`stage1_run2_icarus_overlay_NuGraphReco.fcl`), the flow is:

1. Pandora runs once to make slices
2. **`NGMultiSliceCryoE/W`** scores the hits of the non-clear-cosmic slices and writes the hits that pass the filter (`ngfilteredhitsCryoE/W`)
3. Pandora runs **again from scratch** (including slicing) on those filtered hits, and that's what the CAFs are made from

The analyzer reads step 2: the scores, the slices they belong to, and the first-pass Pandora hits.

---

## 📁 What's here

| Path | What it is |
|---|---|
| `fcl/nugraph_analyser_wtruth_v2.fcl` | **The one to use.** Truth-inclusive analyzer on `NGMultiSliceCryoE`, with overlay-safe truth matching. CryoE only for now |
| `fcl/nugraph_analyser.fcl` | Predictions only, no truth. Used by the event viewer |
| `xml/BNB_NuGraph_v10_06_00_04p04_test.xml` | Grid XML for the stage1 and caf stages (BNB MC + overlay) |
| `xml/BNB_NuGraph_v10_06_00_04p04_NG2Stage1Analyser.xml` | Grid XML that runs the analyzer over existing stage-1 files |
| `analyzer/ICARUSNuGraphAnalyzer_module.cc` | **Reference copy** of the analyzer. See the note below |
| `notebooks/` | Event viewer, truth checks, confusion matrices (table further down) |

> 📌 **About the analyzer copy.** The module only builds inside icaruscode, so the copy here is for reading. The version that builds lives in [`sparshitadey/icaruscode`](https://github.com/sparshitadey/icaruscode), branch `feature/sdey_NuGraph2_TruthAnalyzer_v04p04`, commit `617e12ef`. It's a modified version of the original analyzer by Riccardo Triozzi and Giuseppe Cerati.

---

## 🔧 What the truth-inclusive analyzer adds

- **Reco fixes**: uninitialised cluster/PFP indices, an empty `NuGraphEventOutput`, and scores leaking from one hit to the next
- **Truth semantic labels**: `BuildSemanticLabels()`, ported from pynuml (`labels/standard.py` + `simple.py`), so truth uses the same class definitions as training
- **Overlay-safe truth matching**: reads the precomputed `Hit ↔ MCParticle` associations (`HitTruthAssnsLabel: "mcassociationsGausCryoE"`) instead of live backtracking
- **Extra trees**: `NuGraphEdepOutput`, `NuGraphSliceOutput`, `NuGraphPFPOutput`, alongside `NuGraphHitOutput` and `NuGraphEventOutput`

Truth sentinels in `truth_semantic`:

| Value | Meaning |
|---|---|
| `0–4` | MIP, HIP, Shower, Michel, Diffuse |
| `-1` | no truth match (on overlay samples this is most hits: they're real cosmic data) |
| `-2` | broken particle chain, should be ~0 |
| `5` | invisible particle, should be small |

> 💡 **Overlay samples and truth.** Live `BackTrackerService` returns **0% truth match** on BNB + overlay stage-1 files, with no error or warning. The `mcassociationsGausCryoE` route fixes this. On overlay, ~10% of hits being truth-matched is expected (they're the simulated neutrino), not a bug. `notebooks/bnb_overlay_truth_analysis.ipynb` explains this in its intro.

---

## ▶️ Running it

**Build** the analyzer in an icaruscode v10_06_00_04p04 mrb area with the branch above checked out.

**Locally**, on a handful of files:

```bash
lar -c fcl/nugraph_analyser_wtruth_v2.fcl -S my_stage1_files.list -T hist_nugraph_analyser_truth_mysample.root
```

**On the grid**, use `xml/BNB_NuGraph_v10_06_00_04p04_NG2Stage1Analyser.xml` with a tarball of that mrb area, and point its input list at your stage-1 files.

> ⚠️ Stage-1 files usually live on `/pnfs/.../scratch`, which gets purged. Run the analyzer before they disappear.

> ⚠️ Read the scores with `ICARUSNuGraphAnalyzer`, not with a gallery/ROOT macro that re-joins hits from separate collections. Re-joining by channel and time mixes up which score belongs to which hit, so the plots look plausible but are wrong.

---

## 📓 Notebooks

| Notebook | Reads | Shows |
|---|---|---|
| `nugraph_event_viewer.ipynb` | `nugraph_analyser.fcl` output | Event and slice overview, predicted class per hit, filter score (signal vs noise). Set `PRED_FILE` at the top |
| `truth_analysis.ipynb` | `nugraph_analyser_wtruth*.fcl` output, training sample | Truth-match checks, label distributions, confusion matrices, plus an exact-enumeration bootstrap cross-check for small samples |
| `bnb_overlay_truth_analysis.ipynb` | `nugraph_analyser_wtruth_v2.fcl` output, BNB + overlay | Same checks for overlay, recall / precision matrices, 68% event-level bootstrap intervals |
| `bnb_overlay_truth_analysis_streaming.ipynb` | analyzer output from the 5000-job BNB run | The BNB analysis scaled up, streaming the outputs instead of loading them all at once (test mode) |
| `plot_recofix_comparison.ipynb` | analyzer + RecoFix variant output | Before/after plots for the analyzer reco fixes |

Point each notebook at your own output by editing the input path near the top. The saved outputs are left in as examples.
