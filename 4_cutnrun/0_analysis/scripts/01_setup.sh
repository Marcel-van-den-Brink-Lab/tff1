#!/bin/bash
#SBATCH --job-name=setup
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/01_setup_%j.out
#SBATCH --error=logs/01_setup_%j.err

module load Nextflow

cd "$SLURM_SUBMIT_DIR"

GENOME="mm10"

IGENOMES_BASE="igenomes"
IGENOMES_PATH="${IGENOMES_BASE}/Mus_musculus/UCSC/${GENOME}"

if [ -d "${IGENOMES_PATH}/Sequence/WholeGenomeFasta" ]; then
    echo "iGenomes mm10 already exists, skipping."
else
    mkdir -p "${IGENOMES_PATH}"
    TARBALL="Mus_musculus_UCSC_${GENOME}.tar.gz"
    wget -c -O "${IGENOMES_BASE}/${TARBALL}" \
        "https://s3.amazonaws.com/igenomes.illumina.com/Mus_musculus/UCSC/${GENOME}/${TARBALL}"
    tar -xzf "${IGENOMES_BASE}/${TARBALL}" -C "${IGENOMES_BASE}"
    rm -f "${IGENOMES_BASE}/${TARBALL}"
    echo "iGenomes mm10 download complete."
fi

BLACKLIST_DIR="references/blacklist"
BLACKLIST_FILE="${BLACKLIST_DIR}/mm10-blacklist.v2.bed.gz"

if [ -f "${BLACKLIST_FILE}" ]; then
    echo "Blacklist already exists, skipping."
else
    mkdir -p "${BLACKLIST_DIR}"
    wget -c -O "${BLACKLIST_FILE}" \
        "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz"
    echo "Blacklist download complete."
fi
