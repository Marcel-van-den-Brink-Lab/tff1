#!/bin/bash

#SBATCH --job-name=atac_global_heatmap
#SBATCH --output=/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/0_broadpeak/logs/global_heatmap_%j.out
#SBATCH --error=/coh_labs/mvandenbrink/users/pkaur/6_tff1/2_bulk_atac/atacseq_pipeline/0_broadpeak/logs/global_heatmap_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
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
MATRIX="${OUT_DIR}/global.mat.gz"
HEATMAP_PDF="${OUT_DIR}/global_accessibility_heatmap.pdf"

GFP_BW="${BIGWIG_DIR}/GFP_REP1.mLb.clN.bigWig"
TFF1_BW="${BIGWIG_DIR}/TFF1_REP1.mLb.clN.bigWig"

THREADS=${SLURM_CPUS_PER_TASK:-8}


source /coh_labs/mvandenbrink/users/pkaur/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

echo "======================================================"
echo "  ATAC-seq Global Accessibility Heatmap (side-by-side)"
echo "  Date      : $(date)"
echo "  Host      : $(hostname)"
echo "  Threads   : ${THREADS}"
echo "======================================================"

for f in "${GFP_BW}" "${TFF1_BW}" "${DIFFBIND_CSV}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: File not found: ${f}" >&2; exit 1
    fi
done

if [[ ! -f "${BLACKLIST}" ]]; then
    echo "  WARNING: Blacklist not found — continuing without it."
    BLACKLIST_FLAG=""
else
    BLACKLIST_FLAG="--blackListFileName ${BLACKLIST}"
fi

echo ""
echo "------------------------------------------------------"
echo "  Step 1: Build consensus peaks BED  ($(date '+%H:%M:%S'))"
echo "------------------------------------------------------"

awk -F',' '
    NR==1 { next }          # skip header
    {
        chr = $1; gsub(/"/, "", chr)   # remove surrounding quotes
        start = $2
        end   = $3
        print chr "\t" start "\t" end
    }
' "${DIFFBIND_CSV}" \
  | sort -k1,1 -k2,2n \
  > "${CONSENSUS_BED}"

N_PEAKS=$(wc -l < "${CONSENSUS_BED}")
echo "  Consensus peaks: ${N_PEAKS} regions → ${CONSENSUS_BED}"

echo ""
echo "------------------------------------------------------"
echo "  Step 2: computeMatrix (GFP + TFF1)  ($(date '+%H:%M:%S'))"
echo "------------------------------------------------------"

computeMatrix reference-point \
    --referencePoint center \
    --scoreFileName "${GFP_BW}" "${TFF1_BW}" \
    --regionsFileName "${CONSENSUS_BED}" \
    --outFileName "${MATRIX}" \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 2000 \
    --binSize 10 \
    --samplesLabel "GFP" "TFF1" \
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
    fh.readline()                               # skip JSON header
    for line in fh:
        for v in line.strip().split("\t")[6:]:  # first 6 cols = region metadata
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
    --colorMap Blues Reds \
    --missingDataColor "#FFFFFF" \
    --zMin 0 \
    --zMax "${ZMAX}" \
    --whatToShow "plot and heatmap" \
    --heatmapHeight 15 \
    --heatmapWidth 5 \
    --xAxisLabel "" \
    --refPointLabel "Peak center" \
    --regionsLabel "Consensus peaks" \
    --samplesLabel "GFP" "TFF1" \
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
echo "  Color scale: zMin=0  zMax=${ZMAX}"
echo "======================================================"

exit 0
