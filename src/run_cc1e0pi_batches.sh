#!/bin/bash
# Batch driver for CC1e0piSelection_MakeAll.C
#
# Runs the cafe macro once per 50-file batch list in BATCH_DIR, generating a
# per-batch macro copy (input list + output file substituted) so the original
# source is never touched. Skips any batch whose output already exists, so a
# crash or Ctrl-C only costs the batch in flight, not the whole run.
#
# Usage:
#   1. Get a kerberos ticket / interactive auth first (same as any cafe run).
#   2. stdbuf -oL -eL ./run_cc1e0pi_batches.sh > driver.log 2>&1 & disown
#   3. tail -f driver.log

set -u  # (not -e: one failed batch shouldn't kill the whole loop)

# Set to 1 to delete the generated per-batch macro copies (src/batches/*.C)
# at the end of a run with zero failures. Off by default -- they're tiny
# (~15-20 KB each, ~2 MB total for 95 batches) and useful as a record of
# exactly which list/output pair each batch ran against, e.g. for re-running
# a single failed batch later. Only the batch-level .root/.log outputs are
# ever worth cleaning up for space; the macros aren't.
CLEANUP_MACROS=${CLEANUP_MACROS:-0}

# Set to a positive integer to stop after processing that many batches this
# run (counting only batches actually attempted, not ones already done and
# skipped). Useful for a quick single-batch sanity check before committing
# to the full set, e.g.: MAX_BATCHES=1 ./run_cc1e0pi_batches.sh
# 0 (default) = no limit, process everything in BATCH_DIR.
MAX_BATCHES=${MAX_BATCHES:-0}

MACRO_SRC=/exp/icarus/data/users/sdey2/NuGraphReco/cc1e0pi_v3/src/CC1e0piSelection_MakeAll.C
BATCH_DIR=/exp/icarus/data/users/sdey2/NuGraphReco/cc1e0pi_v3/input/cc1e0pi_batches

# EDIT ME: path to your sbnana_v04p04 local products area (the one built
# specifically for v04p04 CAF compatibility -- see the localProducts*/setup
# you source interactively). Update this once and the environment setup
# below no longer depends on the launching shell already having it sourced.
SBNANA_AREA=/exp/icarus/data/users/sdey2/NuGraphReco/sbnana_v04p04

echo "=== Setting up CAFAna environment ==="
# setup_icarus.sh / localProducts*/setup / mrbslp reference unset variables
# internally (e.g. Options[@]) and weren't written to be nounset-safe --
# relax set -u just for this block, restore it right after.
set +u
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source "${SBNANA_AREA}"/localProducts*/setup
mrbslp
set -u

if ! command -v cafe >/dev/null 2>&1; then
    echo "ERROR: cafe not found on PATH after sourcing environment."
    echo "  Checked SBNANA_AREA=${SBNANA_AREA}"
    echo "  Confirm this path and its localProducts*/setup are correct, then retry."
    exit 1
fi
echo "cafe resolved to: $(command -v cafe)"

# Date-stamped output dir for .root outputs and logs (app is out of space,
# so these go under data instead of alongside src/).
RUN_DATE=$(date +%Y%m%d)
RUN_BASE=/exp/icarus/data/users/sdey2/NuGraphReco/cc1e0pi_v3/output/data/${RUN_DATE}
OUT_DIR="${RUN_BASE}/root"
LOG_DIR="${RUN_BASE}/logs"

# Generated per-batch macros stay NEXT TO src/ (not under the date-stamped
# output dir) -- CC1e0piSelection_MakeAll.C uses relative includes
# ("../include/Config.h" etc.) that resolve based on where cafe is run from.
# Putting the macro copies anywhere else breaks those includes.
MACRO_DIR=/exp/icarus/data/users/sdey2/NuGraphReco/cc1e0pi_v3/src/batches

mkdir -p "$OUT_DIR" "$LOG_DIR" "$MACRO_DIR"
echo "Run date: $RUN_DATE"
echo "Root/logs under: $RUN_BASE"
echo "Generated macros under: $MACRO_DIR"

# ---- sanity check the substitution targets exist exactly once (active line) ----
N_LISTFILE_ACTIVE=$(grep -cE '^\s*const std::string listFile = ".*";' "$MACRO_SRC")
N_FOUT=$(grep -cE 'TFile FOut\(".*", "RECREATE"\);' "$MACRO_SRC")
if [ "$N_LISTFILE_ACTIVE" -ne 1 ] || [ "$N_FOUT" -ne 1 ]; then
    echo "ERROR: expected exactly one active listFile line and one TFile FOut line in $MACRO_SRC"
    echo "  active listFile matches: $N_LISTFILE_ACTIVE"
    echo "  TFile FOut matches:      $N_FOUT"
    echo "Macro source may have changed -- check before running batches."
    exit 1
fi

N_TOTAL=0
N_DONE=0
N_RUN=0
N_FAIL=0

for LIST in "$BATCH_DIR"/batch_*.list; do
    N_TOTAL=$((N_TOTAL+1))
    BID=$(basename "$LIST" .list)                          # e.g. batch_000
    OUTFILE="$OUT_DIR/CC1e0piSelection_All_${BID}.root"
    LOGFILE="$LOG_DIR/${BID}.log"
    MACRO="$MACRO_DIR/CC1e0piSelection_MakeAll_${BID}.C"

    if [ -f "$OUTFILE" ]; then
        echo "[$BID] already done, skipping"
        N_DONE=$((N_DONE+1))
        continue
    fi

    if [ "$MAX_BATCHES" -gt 0 ] && [ $((N_RUN + N_FAIL)) -ge "$MAX_BATCHES" ]; then
        echo "[$BID] MAX_BATCHES=$MAX_BATCHES reached -- stopping here (remaining batches untouched)"
        break
    fi

    # Only touch the ACTIVE (non-commented) listFile line; leave the two
    # commented example lines untouched. -E for extended regex (GNU sed).
    # Also fix #include "../include/..." paths (see note above), AND rename
    # the function declaration to match the new filename -- cafe/ROOT auto-
    # invokes a function matching the file's basename (visible in the log as
    # "Processing .../CC1e0piSelection_MakeAll_batch_NNN.C++()"), so renaming
    # only the file without renaming the function inside it breaks the
    # auto-invoke and fails with "undeclared identifier".
    sed -E \
        -e "s|^(\s*)const std::string listFile = \".*\";|\1const std::string listFile = \"${LIST}\";|" \
        -e "s|TFile FOut\(\".*\", \"RECREATE\"\);|TFile FOut(\"${OUTFILE}\", \"RECREATE\");|" \
        -e 's|#include "\.\./include/|#include "../../include/|g' \
        -e "s|void CC1e0piSelection_MakeAll\(\)|void CC1e0piSelection_MakeAll_${BID}()|" \
        "$MACRO_SRC" > "$MACRO"

    # Verify the substitutions actually landed before burning time on a run
    if ! grep -qF "$LIST" "$MACRO" || ! grep -qF "$OUTFILE" "$MACRO"; then
        echo "[$BID] ERROR: listFile/output sed substitution did not take -- skipping, check $MACRO by hand"
        N_FAIL=$((N_FAIL+1))
        continue
    fi
    if grep -qE '#include "\.\./include/' "$MACRO"; then
        echo "[$BID] ERROR: include path fix did not take (still has ../include/) -- skipping, check $MACRO by hand"
        N_FAIL=$((N_FAIL+1))
        continue
    fi
    if ! grep -qF "void CC1e0piSelection_MakeAll_${BID}()" "$MACRO"; then
        echo "[$BID] ERROR: function rename did not take (won't match cafe's auto-invoke) -- skipping, check $MACRO by hand"
        N_FAIL=$((N_FAIL+1))
        continue
    fi

    NFILES=$(wc -l < "$LIST")
    echo "[$BID] $NFILES files -> $OUTFILE"
    ( cd "$(dirname "$MACRO")" && stdbuf -oL -eL cafe -bq "$MACRO" ) > "$LOGFILE" 2>&1

    if [ -f "$OUTFILE" ]; then
        echo "[$BID] done"
        N_RUN=$((N_RUN+1))
    else
        echo "[$BID] FAILED -- check $LOGFILE"
        N_FAIL=$((N_FAIL+1))
    fi
done

echo ""
echo "=== summary ==="
echo "total batches:      $N_TOTAL"
echo "already done:       $N_DONE"
echo "run this session:   $N_RUN"
echo "failed:              $N_FAIL"

if [ "$CLEANUP_MACROS" -eq 1 ]; then
    if [ "$N_FAIL" -eq 0 ]; then
        echo ""
        echo "CLEANUP_MACROS=1 and zero failures -- removing generated macros in $MACRO_DIR"
        rm -f "$MACRO_DIR"/CC1e0piSelection_MakeAll_batch_*.C
    else
        echo ""
        echo "CLEANUP_MACROS=1 but $N_FAIL batch(es) failed -- leaving all generated macros in place for debugging"
    fi
fi
