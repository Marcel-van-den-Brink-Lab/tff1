#!/bin/bash
#SBATCH --job-name=cutandrun_nf
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/02_run_pipeline_%j.out
#SBATCH --error=logs/02_run_pipeline_%j.err
#SBATCH --mail-type=END,FAIL

module load Nextflow
module load singularity

cd "$SLURM_SUBMIT_DIR"

nextflow run nf-core/cutandrun \
    -r 3.2.2 \
    -profile singularity \
    -c nextflow.config \
    -params-file params.yaml \
    -resume