#!/bin/bash
#SBATCH --job-name=avg_bigwigs
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/06_avg_bigwigs_%j.out
#SBATCH --error=logs/06_avg_bigwigs_%j.err
#SBATCH --mail-type=END,FAIL

source /coh_labs/mvandenbrink/users/pkaur/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

cd "$SLURM_SUBMIT_DIR"
PROJECT_DIR="$(pwd)"

MANIFEST="${PROJECT_DIR}/0_samplesheets/sample_manifest.tsv"
BIGWIG_DIR="${PROJECT_DIR}/1_nf_output/03_peak_calling/03_bed_to_bigwig"
BW_PREFIXES=("")
OUT_DIR="${PROJECT_DIR}/8_trackplots/avg_bigwigs"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

TFS=$(      tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $3}' | grep -iv 'IgG' | sort -u)
CONDITIONS=$(tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $4}' | grep -iv 'IgG' | sort -u)

bw_path() {
    local sample_id="$1"
    local rep="${sample_id##*_}"
    local group="${sample_id%_*}"
    local bw
    for prefix in "${BW_PREFIXES[@]}"; do
        bw="${BIGWIG_DIR}/${prefix}${group}_R${rep}.bigWig"
        [[ -f "$bw" ]] && echo "$bw" && return 0
    done
    for bw in "${BIGWIG_DIR}/"*"${group}_R${rep}.bigWig"; do
        [[ -f "$bw" ]] && echo "$bw" && return 0
    done
    return 1
}

echo "=== Step 1: Checking bigWig files ==="
missing_bw=0
while IFS=$'\t' read -r sample_id _rest; do
    [[ "$sample_id" == "sample_id" ]] && continue
    [[ "$sample_id" =~ [Ii]g[Gg] ]] && continue
    if ! bw_path "${sample_id}" > /dev/null; then
        echo "  WARNING: No bigWig found for ${sample_id} in ${BIGWIG_DIR}/"
        (( missing_bw++ ))
    fi
done < "${MANIFEST}"

if [[ $missing_bw -gt 0 ]]; then
    echo "  ${missing_bw} bigWig file(s) missing. Check BIGWIG_DIR and BW_PREFIXES."
    echo "  Continuing; affected TF × condition groups will be skipped."
else
    echo "  All bigWig files found."
fi

echo ""
echo "=== Step 2: Averaging replicate bigWigs ==="

for tf in $TFS; do
    for condition in $CONDITIONS; do
        avg_bw="${OUT_DIR}/${tf}_${condition}.avg.bigWig"

        if [[ -f "$avg_bw" ]]; then
            echo "  Exists, skipping: $(basename "${avg_bw}")"
            continue
        fi

        rep_bws=()
        while IFS=$'\t' read -r sample_id _name sample_tf group; do
            [[ "$sample_id" == "sample_id" ]] && continue
            [[ "$sample_tf" != "$tf"        ]] && continue
            [[ "$group"     != "$condition" ]] && continue
            bw="$(bw_path "${sample_id}")"
            [[ -n "$bw" ]] && rep_bws+=("$bw")
        done < "${MANIFEST}"

        if [[ ${#rep_bws[@]} -eq 0 ]]; then
            echo "  WARNING: No valid bigWigs for ${tf}_${condition} — skipping"
            continue
        fi

        echo "  ${tf}_${condition}: averaging ${#rep_bws[@]} replicate(s)..."
        for bw in "${rep_bws[@]}"; do
            echo "    ${bw}"
        done

        bigwigAverage \
            --bigwigs            "${rep_bws[@]}" \
            --outFileName        "${avg_bw}" \
            --numberOfProcessors "${THREADS}"

        if [[ $? -eq 0 ]]; then
            echo "  Saved: $(basename "${avg_bw}")"
        else
            echo "  ERROR: bigwigAverage failed for ${tf}_${condition}"
        fi
    done
done

echo ""
echo "Done."
echo "  Averaged bigWigs: ${OUT_DIR}/"
ls -lh "${OUT_DIR}/"*.bigWig 2>/dev/null || echo "  (no files yet)"
