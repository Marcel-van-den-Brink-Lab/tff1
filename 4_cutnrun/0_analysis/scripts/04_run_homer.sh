#!/bin/bash
#SBATCH --job-name=homer
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04_homer_%j.out
#SBATCH --error=logs/04_homer_%j.err
#SBATCH --mail-type=END,FAIL

source /coh_labs/mvandenbrink/users/pkaur/miniconda3/etc/profile.d/conda.sh
conda activate cutnrun

cd "$SLURM_SUBMIT_DIR"
PROJECT_DIR="$(pwd)"

GENOME="mm10"
PEAK_DIR="${PROJECT_DIR}/1_nf_output/03_peak_calling/05_consensus_peaks"
BAM_DIR="${PROJECT_DIR}/2_aligned"
TAG_DIR="${PROJECT_DIR}/6_homer/tag_dirs"
MOTIF_DIR="${PROJECT_DIR}/6_homer/motifs"
ANNOT_DIR="${PROJECT_DIR}/6_homer/annotation"
MOTIF_GENE_DIR="${PROJECT_DIR}/6_homer/motif_gene_links"

echo "=== Building tag directories ==="
for bam in "${BAM_DIR}"/*.final.bam; do
    sample=$(basename "${bam}" .final.bam)
    if [ ! -d "${TAG_DIR}/${sample}" ]; then
        makeTagDirectory "${TAG_DIR}/${sample}" \
            "${bam}" \
            -genome "${GENOME}" \
            -checkGC
        echo "  Tag dir created: ${sample}"
    else
        echo "  Tag dir exists, skipping: ${sample}"
    fi
done

echo ""
echo "=== Running motif enrichment ==="
for bed in "${PEAK_DIR}"/*peak_counts.bed; do
    sample=$(basename "${bed}" .bed)
    out="${MOTIF_DIR}/${sample}"
    if [ ! -d "${out}" ]; then
        mkdir -p "${out}"
        findMotifsGenome.pl \
            "${bed}" \
            "${GENOME}" \
            "${out}" \
            -size 200 \
            -mask \
            -p "${SLURM_CPUS_PER_TASK}"
        echo "  Motif analysis done: ${sample}"
    else
        echo "  Motif dir exists, skipping: ${sample}"
    fi
done

echo ""
echo "=== Annotating peaks ==="
for bed in "${PEAK_DIR}"/*peak_counts.bed; do
    sample=$(basename "${bed}" .bed)
    if [ ! -f "${ANNOT_DIR}/${sample}_annotated.txt" ]; then
        annotatePeaks.pl \
            "${bed}" \
            "${GENOME}" \
            -annStats "${ANNOT_DIR}/${sample}_annStats.txt" \
            > "${ANNOT_DIR}/${sample}_annotated.txt"
        echo "  Annotation done: ${sample}"
    else
        echo "  Annotation exists, skipping: ${sample}"
    fi
done

echo ""
echo "=== Linking motifs to genes ==="
mkdir -p "${MOTIF_GENE_DIR}"

for bed in "${PEAK_DIR}"/*peak_counts.bed; do
    sample=$(basename "${bed}" .bed)

    known_dir="${MOTIF_DIR}/${sample}/knownResults"
    combined_known="${known_dir}/combined_known.motif"
    if [ ! -f "${MOTIF_GENE_DIR}/${sample}_known_motif_genes.txt" ]; then
        if [ -d "${known_dir}" ] && compgen -G "${known_dir}/known*.motif" > /dev/null 2>&1; then
            cat "${known_dir}"/known*.motif > "${combined_known}"
            annotatePeaks.pl \
                "${bed}" \
                "${GENOME}" \
                -m "${combined_known}" \
                -nmotifs \
                > "${MOTIF_GENE_DIR}/${sample}_known_motif_genes.txt"
            echo "  Known motif-gene links done: ${sample}"
        else
            echo "  WARNING: No known*.motif files in ${known_dir}, skipping: ${sample}"
        fi
    else
        echo "  Known motif-gene links exist, skipping: ${sample}"
    fi

    denovo_motif="${MOTIF_DIR}/${sample}/homerMotifs.all.motifs"
    if [ ! -f "${MOTIF_GENE_DIR}/${sample}_denovo_motif_genes.txt" ]; then
        if [ -f "${denovo_motif}" ]; then
            annotatePeaks.pl \
                "${bed}" \
                "${GENOME}" \
                -m "${denovo_motif}" \
                -nmotifs \
                > "${MOTIF_GENE_DIR}/${sample}_denovo_motif_genes.txt"
            echo "  De novo motif-gene links done: ${sample}"
        else
            echo "  WARNING: De novo motif file not found, skipping: ${sample}"
        fi
    else
        echo "  De novo motif-gene links exist, skipping: ${sample}"
    fi
done

echo ""
echo "HOMER analysis complete."
echo "  Tag directories : ${TAG_DIR}/"
echo "  Motif results   : ${MOTIF_DIR}/"
echo "  Annotations     : ${ANNOT_DIR}/"
echo "  Motif-gene links: ${MOTIF_GENE_DIR}/"
