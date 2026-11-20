#!/bin/bash
# Parallel batch driver for CC1e0piSelection_MakeAll.C (cafe / CAFAna).
#
# Splits a flat-CAF file list into batches and runs N cafe jobs at once, each on its
# own copy of the macro (input list, output file and function name substituted), so the
# original macro is never touched. Batches whose output already exists are skipped, so
# rerunning the same command resumes after a crash or Ctrl-C.
#
# Run inside the SL7 container with sbnana/cafe set up and a valid token:
#   ./run_cc1e0pi_parallel.sh <tag> <file.list> [files_per_batch=50] [n_parallel=6]
# e.g.
#   nohup ./run_cc1e0pi_parallel.sh ClaraOR_1114 \
#     /exp/icarus/data/users/sdey2/NuGraphReco/ClaraModels/lists/caf_flat_ClaraOR_1114.list \
#     > driver_ClaraOR_1114.log 2>&1 &
#   grep -cE "\] ok" driver_ClaraOR_1114.log     # batches finished
#
# Outputs: $TOP/output/root/batches_<tag>/CC1e0piSelection_All_<tag>_bNNN.root
set -u

TAG=${1:?usage: $0 <tag> <file.list> [files_per_batch] [n_parallel]}
LIST=${2:?usage: $0 <tag> <file.list> [files_per_batch] [n_parallel]}
NPER=${3:-50}
NPAR=${4:-6}
TIMEOUT=${TIMEOUT:-90m}          # kill a batch that hangs (e.g. on a stuck dCache read)

TOP=${TOP:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}   # repo root (override with TOP=...)
SRC=${SRC:-$TOP/src/CC1e0piSelection_MakeAll.C}
WORK=$TOP/batch_$TAG             # same depth as src/, so the macro's ../include/ still resolves
OUT=$TOP/output/root/batches_$TAG

command -v cafe >/dev/null || { echo "cafe is not on PATH: set up sbnana first"; exit 1; }
[ -f "$SRC" ]  || { echo "macro not found: $SRC"; exit 1; }
[ -f "$LIST" ] || { echo "list not found: $LIST"; exit 1; }
mkdir -p "$WORK/lists" "$WORK/logs" "$OUT"

# Split the list once (re-running reuses the same split so resume works)
if [ -z "$(ls -A "$WORK/lists")" ]; then
    grep -v '^\s*#' "$LIST" | grep -v '^\s*$' | split -d -a 3 -l "$NPER" - "$WORK/lists/b"
fi
NB=$(ls "$WORK/lists" | wc -l)
echo "tag=$TAG  files=$(cat "$WORK"/lists/b* | wc -l)  batches=$NB  per_batch=$NPER  parallel=$NPAR"
echo "macro=$SRC"
echo "work=$WORK"
echo "out=$OUT"

# Keep the bearer token fresh for xrootd reads during a long run
( while true; do htgettoken -a htvaultprod.fnal.gov -i icarus </dev/null >/dev/null 2>&1; sleep 1800; done ) </dev/null >/dev/null 2>&1 &
TOKPID=$!
trap 'pkill -P $TOKPID 2>/dev/null; kill $TOKPID 2>/dev/null' EXIT

run_one() {
    local B=$1
    local LISTB=$WORK/lists/$B
    local OUTF=$OUT/CC1e0piSelection_All_${TAG}_${B}.root
    local FN=MakeAll_${TAG}_${B}
    local M=$WORK/$FN.C
    if [ -s "$OUTF" ]; then echo "[$B] already done"; return 0; fi

    sed -e "s#^\([[:space:]]*const std::string listFile = \)\"[^\"]*\";#\1\"$LISTB\";#" \
        -e "s#^\([[:space:]]*TFile FOut(\)\"[^\"]*\"#\1\"$OUTF.tmp\"#" \
        -e "s#void CC1e0piSelection_MakeAll()#void $FN()#" "$SRC" > "$M"
    if ! { grep -q "\"$LISTB\"" "$M" && grep -q "\"$OUTF.tmp\"" "$M" && grep -q "void $FN()" "$M"; }; then
        echo "[$B] FAILED: substitution did not match (check listFile / TFile FOut / function lines)"; return 1
    fi

    local t0=$(date +%s)
    ( cd "$WORK" && timeout "$TIMEOUT" cafe -bq "$FN.C" 2>&1 | grep -v "Bad value" > "logs/$B.log" )
    if [ -s "$OUTF.tmp" ] && grep -q "Filling" "$WORK/logs/$B.log" && ! grep -qE "^Error|segmentation|\*\*\* Break" "$WORK/logs/$B.log"; then
        mv "$OUTF.tmp" "$OUTF"
        rm -f "$WORK/${FN}_C"* "$WORK/$FN.C"          # ACLiC build products and the macro copy
        echo "[$B] ok  ($(( ($(date +%s) - t0) / 60 )) min)"
    else
        rm -f "$OUTF.tmp"
        echo "[$B] FAILED after $(( ($(date +%s) - t0) / 60 )) min -- see $WORK/logs/$B.log"
    fi
}
export -f run_one
export WORK OUT TAG SRC TIMEOUT

ls "$WORK/lists" | xargs -P "$NPAR" -I{} bash -c 'run_one {}'

NOK=$(ls "$OUT"/CC1e0piSelection_All_${TAG}_b*.root 2>/dev/null | wc -l)
echo "finished: $NOK / $NB batches have output in $OUT"
[ "$NOK" -eq "$NB" ] || echo "re-run the same command to retry the failed batches"
