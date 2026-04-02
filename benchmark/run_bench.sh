#!/bin/bash
# CoverM benchmark script
#
# Measures wall time and peak RAM for a given SIF against the CL4-AERO dataset.
# Runs four scenarios sequentially:
#   1. no BED    — mapping + basic metrics only
#   2. BED       — + --regions-bed + --regions-bed-unlabeled
#   3. BED+BG    — same + --output-bedcov (bedGraph output)
#
# Usage:
#   ./run_bench.sh --sif /path/to/coverm.sif [--label v0.11.0-myfeature]
#                  [--scenarios nobed,bed,bed_bg] [--threads 52]
#
# Output:
#   OUTDIR/bench_<label>_<scenario>.{tsv,log}
#   Summary printed to stdout

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults — edit to match your environment
# ---------------------------------------------------------------------------
DEFAULT_SINGLE="/HDD/Raw_data/WGS21_Bact_Bioaster/ONT/Preprocessed/CL4-AERO.filtered_lysed.fastq.gz"
DEFAULT_MMI="/HDD/cdemay/benchmark/wgs21/sfbio-mag-catalogue/v0.1.0_20260401-110054_benchmark/03_catalogue/GLOBAL/dad9f0975f2a2116ae49f995b21ebea8.catalogue.fasta.gz.mmi"
DEFAULT_BED="/HDD/cdemay/REPOSITORIES/sfbio-coverage/research/results/14_cross_spatial/annot_regions.bed"
DEFAULT_OUTDIR="/HDD/cdemay/REPOSITORIES/sfbio-coverage/research/results/14_cross_spatial"
DEFAULT_THREADS=52
DEFAULT_SCENARIOS="nobed,bed,bed_bg"
# ---------------------------------------------------------------------------

SIF=""
LABEL=""
SINGLE="$DEFAULT_SINGLE"
MMI="$DEFAULT_MMI"
BED="$DEFAULT_BED"
OUTDIR="$DEFAULT_OUTDIR"
THREADS="$DEFAULT_THREADS"
SCENARIOS="$DEFAULT_SCENARIOS"

usage() {
    echo "Usage: $0 --sif <sif> [--label <label>] [--single <fastq>] [--mmi <mmi>]"
    echo "          [--bed <bed>] [--outdir <dir>] [--threads N]"
    echo "          [--scenarios nobed,bed,bed_bg]"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --sif)      SIF="$2";       shift 2 ;;
        --label)    LABEL="$2";     shift 2 ;;
        --single)   SINGLE="$2";    shift 2 ;;
        --mmi)      MMI="$2";       shift 2 ;;
        --bed)      BED="$2";       shift 2 ;;
        --outdir)   OUTDIR="$2";    shift 2 ;;
        --threads)  THREADS="$2";   shift 2 ;;
        --scenarios) SCENARIOS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

[ -z "$SIF" ] && { echo "ERROR: --sif is required"; usage; }

# Auto-label from SIF filename if not provided
if [ -z "$LABEL" ]; then
    LABEL=$(basename "$SIF" .sif)
fi

BGDIR="${OUTDIR}/bedgraph_bench_${LABEL}"
mkdir -p "$OUTDIR" "$BGDIR"

COMMON_ARGS=(
    --single "$SINGLE"
    --reference "$MMI"
    --minimap2-reference-is-index
    -s '~'
    --mapper minimap2-lr-hq
    --threads "$THREADS"
    --min-covered-fraction 0
    --min-read-percent-identity 0
    --min-read-aligned-length 0
    --min-read-aligned-percent 0
    --methods relative_abundance mean trimmed_mean covered_bases covered_fraction
)

run_scenario() {
    local scenario="$1"; shift
    local tag="${LABEL}_${scenario}"
    local tsv="${OUTDIR}/bench_${tag}.tsv"
    local log="${OUTDIR}/bench_${tag}.log"

    echo "--- ${tag} ---"
    /usr/bin/time -v \
        apptainer exec "$SIF" coverm genome "${COMMON_ARGS[@]}" "$@" \
        > "$tsv" 2> "$log"
    local exit_code=$?

    local elapsed ram cpu
    elapsed=$(grep "Elapsed (wall clock)" "$log" | awk '{print $NF}')
    ram=$(grep "Maximum resident"       "$log" | awk '{print $NF}')
    cpu=$(grep "Percent of CPU"         "$log" | awk '{print $NF}')
    printf "  exit=%-3s  time=%-10s  cpu=%-8s  ram_kb=%s\n" \
        "$exit_code" "$elapsed" "$cpu" "$ram"
}

echo "========================================"
echo "  CoverM benchmark  — label: $LABEL"
echo "  SIF    : $SIF"
echo "  FASTQ  : $SINGLE"
echo "  MMI    : $MMI"
echo "  threads: $THREADS"
echo "========================================"
echo ""

IFS=',' read -ra SCENARIO_LIST <<< "$SCENARIOS"
for s in "${SCENARIO_LIST[@]}"; do
    case "$s" in
        nobed)
            run_scenario "nobed"
            ;;
        bed)
            run_scenario "bed" \
                --regions-bed "$BED" \
                --regions-bed-unlabeled
            ;;
        bed_bg)
            run_scenario "bed_bg" \
                --regions-bed "$BED" \
                --regions-bed-unlabeled \
                --output-bedcov "$BGDIR"
            ;;
        *)
            echo "WARNING: unknown scenario '$s', skipping"
            ;;
    esac
done

echo ""
echo "========================================"
echo "  Done. Logs in $OUTDIR/bench_${LABEL}_*.log"
echo "========================================"
