#!/bin/bash
# =============================================================================
# run_atacseq.slurm — SLURM submission script for nf-core/atacseq v2.1.2
# Project: atacseq_pipeline | https://nf-co.re/atacseq/2.1.2/
#
# Resources requested here are for the Nextflow HEAD PROCESS ONLY (the
# orchestrator). Individual pipeline tasks (alignment, peak-calling, etc.)
# are submitted as separate SLURM jobs by Nextflow, using resources defined
# in nextflow.config (process labels). The 40 CPU / 128 GB cap in
# nextflow.config and params.yml applies to those child jobs.
#
# Before submitting:
#   1. Fill in --partition and --account below.
#   2. Adjust module names to match your HPC (run: module avail nextflow).
#   3. Verify params.yml (input samplesheet, outdir, genome, read_length).
#   4. Optionally validate: python3 scripts/validate_samplesheet.py samplesheet.csv
# =============================================================================

#SBATCH --job-name=broad_atacseq
#SBATCH --output=logs/atacseq_%j.out
#SBATCH --error=logs/atacseq_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --partition=compute  

mkdir -p logs

module purge
module load Nextflow    # <-- adjust to your HPC module name, e.g. nextflow/23.10.1
module load singularity # <-- adjust to your HPC module name, e.g. singularity/3.9.0


export NXF_SINGULARITY_CACHEDIR="${HOME}/.singularity_cache/nf-core"
mkdir -p "${NXF_SINGULARITY_CACHEDIR}"

# JVM heap for the Nextflow orchestrator (NOT for pipeline tasks).
export NXF_OPTS="-Xms1g -Xmx4g"


export SINGULARITY_CACHEDIR="${HOME}/.singularity_cache/singularity"  
export SINGULARITY_TMPDIR="${TMPDIR:-/tmp}"                           

mkdir -p "${SINGULARITY_CACHEDIR}"
mkdir -p "${SINGULARITY_TMPDIR}"


echo "======================================================"
echo "  nf-core/atacseq v2.1.2"
echo "  Date      : $(date)"
echo "  Host      : $(hostname)"
echo "  SLURM job : ${SLURM_JOB_ID}"
echo "  Work dir  : $(pwd)"
echo "  NXF_SINGULARITY_CACHEDIR : ${NXF_SINGULARITY_CACHEDIR}"
echo "  SINGULARITY_TMPDIR       : ${SINGULARITY_TMPDIR}"
echo "======================================================"


nextflow run nf-core/atacseq \
    -r 2.1.2 \
    -profile singularity \
    -c nextflow.config \
    -params-file params.yml \
    -resume

PIPELINE_EXIT=$?

echo ""
echo "======================================================"
echo "  Run finished at : $(date)"
echo "  Exit code       : ${PIPELINE_EXIT}"
if [[ ${PIPELINE_EXIT} -eq 0 ]]; then
    echo "  Status : SUCCESS"
    echo "  Results: see 'outdir' in params.yml"
    echo "  Start with: results/multiqc/multiqc_report.html"
else
    echo "  Status : FAILED"
    echo "  Check the log above and results/pipeline_info/execution_trace.txt"
    echo "  Re-submit with -resume to retry failed tasks only."
fi
echo "======================================================"

exit ${PIPELINE_EXIT}
