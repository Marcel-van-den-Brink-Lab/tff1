#!/bin/bash

#SBATCH --job-name=atac_global_heatmap_reps
#SBATCH --output=/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/0_broadpeak/logs/global_heatmap_reps_%j.out
#SBATCH --error=/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/0_broadpeak/logs/global_heatmap_reps_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --partition=compute


HPC_BASE="/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/0_broadpeak/results/bwa/merged_library"
BIGWIG_DIR="${HPC_BASE}/bigwig"
BLACKLIST="/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/reference/mm10-blacklist.v2.bed"

# DiffBind results CSV (on HPC — copy from local if needed)
DIFFBIND_CSV="/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/0_broadpeak/diffbind_results/diffbind_all_results.csv"

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
OUT_DIR="${PROJECT_DIR}/3_global_heatmap"
mkdir -p "${OUT_DIR}"

CONSENSUS_BED="${OUT_DIR}/consensus_peaks.bed"
MATRIX="${OUT_DIR}/global_per_replicate.mat.gz"
HEATMAP_PDF="${OUT_DIR}/global_accessibility_per_replicate_heatmap.pdf"

THREADS=${SLURM_CPUS_PER_TASK:-8}


source /coh_labs/mvandenbrink/users/pkaur/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

echo "======================================================"
echo "  ATAC-seq Global Accessibility Heatmap (per replicate)"
echo "  Date      : $(date)"
echo "  Host      : $(hostname)"
echo "  Threads   : ${THREADS}"
echo "======================================================"


GFP_BW_FILES=(  "${BIGWIG_DIR}"/GFP_REP*.mLb.clN.bigWig  )
TFF1_BW_FILES=( "${BIGWIG_DIR}"/TFF1_REP*.mLb.clN.bigWig )
ALL_BW_FILES=( "${GFP_BW_FILES[@]}" "${TFF1_BW_FILES[@]}" )

for arr_name in GFP_BW_FILES TFF1_BW_FILES; do
    eval "arr=( \"\${${arr_name}[@]}\" )"
    if [[ ${#arr[@]} -eq 0 || ! -f "${arr[0]}" ]]; then
        echo "ERROR: No files matched for ${arr_name} — check paths." >&2; exit 1
    fi
done

if [[ ! -f "${DIFFBIND_CSV}" ]]; then
    echo "ERROR: DiffBind CSV not found: ${DIFFBIND_CSV}" >&2; exit 1
fi

if [[ ! -f "${BLACKLIST}" ]]; then
    echo "  WARNING: Blacklist not found — continuing without it."
    BLACKLIST_FLAG=""
else
    BLACKLIST_FLAG="--blackListFileName ${BLACKLIST}"
fi

# Build sample labels by stripping path and suffix
SAMPLE_LABELS=()
for bw in "${ALL_BW_FILES[@]}"; do
    SAMPLE_LABELS+=( "$(basename "${bw}" .mLb.clN.bigWig)" )
done

# Build one colormap per column: Blues x3 for GFP, Reds x3 for TFF1
COLOR_MAPS=()
for _ in "${GFP_BW_FILES[@]}";  do COLOR_MAPS+=( "Blues" ); done
for _ in "${TFF1_BW_FILES[@]}"; do COLOR_MAPS+=( "Reds"  ); done

echo "  GFP  BigWigs : ${#GFP_BW_FILES[@]}  (${GFP_BW_FILES[*]##*/})"
echo "  TFF1 BigWigs : ${#TFF1_BW_FILES[@]}  (${TFF1_BW_FILES[*]##*/})"
echo "  Columns      : ${SAMPLE_LABELS[*]}"
echo "  Color maps   : ${COLOR_MAPS[*]}"


echo ""
echo "------------------------------------------------------"
echo "  Step 1: Build consensus peaks BED  ($(date '+%H:%M:%S'))"
echo "------------------------------------------------------"

if [[ -f "${CONSENSUS_BED}" ]]; then
    echo "  Reusing existing: ${CONSENSUS_BED}"
else
    awk -F',' '
        NR==1 { next }
        {
            chr = $1; gsub(/"/, "", chr)
            start = $2
            end   = $3
            print chr "\t" start "\t" end
        }
    ' "${DIFFBIND_CSV}" \
      | sort -k1,1 -k2,2n \
      > "${CONSENSUS_BED}"
fi

N_PEAKS=$(wc -l < "${CONSENSUS_BED}")
echo "  Consensus peaks: ${N_PEAKS} regions → ${CONSENSUS_BED}"


echo ""
echo "------------------------------------------------------"
echo "  Step 2: computeMatrix (all replicates)  ($(date '+%H:%M:%S'))"
echo "------------------------------------------------------"

computeMatrix reference-point \
    --referencePoint center \
    --scoreFileName "${ALL_BW_FILES[@]}" \
    --regionsFileName "${CONSENSUS_BED}" \
    --outFileName "${MATRIX}" \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 2000 \
    --binSize 10 \
    --samplesLabel "${SAMPLE_LABELS[@]}" \
    --skipZeros \
    ${BLACKLIST_FLAG} \
    --numberOfProcessors "${THREADS}"

echo "  Matrix written: ${MATRIX}"


echo ""
echo "------------------------------------------------------"
echo "  Step 3: Compute color scale  ($(date '+%H:%M:%S'))"
echo "------------------------------------------------------"

ZMAX=$(python3 - "${MATRIX}" <<'PYEOF'
import gzip, sys
import numpy as np

vals = []
with gzip.open(sys.argv[1], "rt") as fh:
    fh.readline()
    for line in fh:
        for v in line.strip().split("\t")[6:]:
            try:
                x = float(v)
                if np.isfinite(x) and x > 0:
                    vals.append(x)
            except ValueError:
                pass

zmax = np.percentile(vals, 99) if vals else 5.0
print(f"{zmax:.4f}")
PYEOF
)

echo "  Color scale: zMin=0  zMax=${ZMAX}"


echo ""
echo "------------------------------------------------------"
echo "  Step 4: plotHeatmap  ($(date '+%H:%M:%S'))"
echo "------------------------------------------------------"

plotHeatmap \
    --matrixFile "${MATRIX}" \
    --outFileName "${HEATMAP_PDF}" \
    --colorMap "${COLOR_MAPS[@]}" \
    --missingDataColor "#FFFFFF" \
    --zMin 0 \
    --zMax "${ZMAX}" \
    --whatToShow "plot and heatmap" \
    --heatmapHeight 16 \
    --heatmapWidth 5 \
    --xAxisLabel "" \
    --refPointLabel "Peak center" \
    --regionsLabel "Consensus peaks" \
    --samplesLabel "${SAMPLE_LABELS[@]}" \
    --plotTitle "Global Chromatin Accessibility" \
    --sortUsing mean \
    --sortRegions descend \
    --legendLocation none \
    --dpi 300

echo "  Heatmap written: ${HEATMAP_PDF}"

echo ""
echo "======================================================"
echo "  Done: $(date)"
echo "  Output: ${HEATMAP_PDF}"
echo "  Peaks shown: ${N_PEAKS}"
echo "  Columns: ${SAMPLE_LABELS[*]}"
echo "  Color scale: zMin=0  zMax=${ZMAX}"
echo "======================================================"

exit 0
