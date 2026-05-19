#!/bin/bash
#SBATCH --job-name=deeptools_heatmaps
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_deeptools_%j.out
#SBATCH --error=logs/05_deeptools_%j.err
#SBATCH --mail-type=END,FAIL

source /coh_labs/mvandenbrink/users/pkaur/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

module load BedTools

cd "$SLURM_SUBMIT_DIR"
PROJECT_DIR="$(pwd)"

MANIFEST="${PROJECT_DIR}/0_samplesheets/sample_manifest.tsv"
BIGWIG_DIR="${PROJECT_DIR}/1_nf_output/03_peak_calling/03_bed_to_bigwig"
BW_PREFIXES=()
PEAK_DIR="${PROJECT_DIR}/3_peaks"
PEAK_STRINGENCY="relaxed"
PEAK_SUFFIX=".seacr.peaks.${PEAK_STRINGENCY}.bed"
UPSTREAM=2000
DOWNSTREAM=2000
REGION_MODE="tss"
GENCODE_GZ="${PROJECT_DIR}/7_deeptools/gencode.vM10.annotation.gtf.gz"
GTF="${PROJECT_DIR}/7_deeptools/gencode.vM10.annotation.gtf"
COLORMAP="RdBu_r"
OUT_DIR="${PROJECT_DIR}/7_deeptools"
CONSENSUS_DIR="${OUT_DIR}/consensus_peaks"
MATRIX_DIR="${OUT_DIR}/matrices"
HEATMAP_DIR="${OUT_DIR}/heatmaps"
TSS_BED="${OUT_DIR}/tss.mm10.bed"
CONSENSUS_BW_DIR="${OUT_DIR}/consensus_bigwigs"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${CONSENSUS_DIR}" "${CONSENSUS_BW_DIR}" "${MATRIX_DIR}" "${HEATMAP_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

TFS=$(      tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $3}' | sort -u)
CONDITIONS=$(tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $4}' | sort -u)

case "$REGION_MODE" in
    peaks|tss|both) ;;
    *) echo "ERROR: REGION_MODE must be 'peaks', 'tss', or 'both' (got '${REGION_MODE}')"; exit 1 ;;
esac

if [[ "$REGION_MODE" == "tss" || "$REGION_MODE" == "both" ]]; then
    if [[ ! -f "$GTF" && -f "$GENCODE_GZ" ]]; then
        echo "Decompressing: $(basename "${GENCODE_GZ}")"
        gunzip -k "${GENCODE_GZ}"
    fi
    if [[ ! -f "$GTF" ]]; then
        echo "ERROR: GTF not found: ${GTF}"
        echo "       Place the GENCODE GTF (or .gz) in: $(dirname "${GTF}")/"
        exit 1
    fi
fi

active_modes=()
[[ "$REGION_MODE" == "peaks" || "$REGION_MODE" == "both" ]] && active_modes+=("peaks")
[[ "$REGION_MODE" == "tss"   || "$REGION_MODE" == "both" ]] && active_modes+=("tss")

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

peak_path() {
    local sample_name="$1"
    local sample_id="$2"
    local rep="${sample_id##*_}"
    local pf

    pf="${PEAK_DIR}/${sample_name}${PEAK_SUFFIX}"
    [[ -f "$pf" ]] && echo "$pf" && return 0

    for pf in "${PEAK_DIR}/${sample_name}"*"REP${rep}"*"${PEAK_SUFFIX}"; do
        [[ -f "$pf" ]] && echo "$pf" && return 0
    done

    for pf in "${PEAK_DIR}/${sample_name}"*"_R${rep}"*"${PEAK_SUFFIX}" \
              "${PEAK_DIR}/${sample_name}"*"_${rep}"*"${PEAK_SUFFIX}"; do
        [[ -f "$pf" ]] && echo "$pf" && return 0
    done

    local -a matches=( "${PEAK_DIR}/${sample_name}"*"${PEAK_SUFFIX}" )
    if [[ ${#matches[@]} -eq 1 && -f "${matches[0]}" ]]; then
        echo "${matches[0]}" && return 0
    fi

    pf="${PEAK_DIR}/${sample_id}${PEAK_SUFFIX}"
    [[ -f "$pf" ]] && echo "$pf" && return 0

    return 1
}

if [[ "$REGION_MODE" == "tss" || "$REGION_MODE" == "both" ]]; then
    echo "=== Step 0: Generating TSS BED from GTF ==="
    if [[ -f "$TSS_BED" ]]; then
        echo "  Exists, skipping: ${TSS_BED}"
    else
        echo "  Parsing: ${GTF}"
        awk 'BEGIN{OFS="\t"} !/^#/ && $3=="transcript" {
            name = ""
            if (match($0, /gene_name "([^"]+)"/, m))   name = m[1]
            else if (match($0, /gene_id "([^"]+)"/, m)) name = m[1]
            else name = "unknown"

            if ($7 == "+") print $1, $4-1, $4, name, 0, $7
            else           print $1, $5-1, $5, name, 0, $7
        }' "${GTF}" \
            | sort -k1,1 -k2,2n \
            | uniq \
            > "${TSS_BED}"

        n_tss=$(wc -l < "${TSS_BED}")
        echo "  Created: ${TSS_BED} (${n_tss} TSS entries)"
        if [[ "$n_tss" -eq 0 ]]; then
            echo "  ERROR: TSS BED is empty — check that the GTF has 'transcript' features."
            exit 1
        fi
    fi
    echo ""
fi

echo "=== Step 1: Checking bigWig files ==="
missing_bw=0
while IFS=$'\t' read -r sample_id _rest; do
    [[ "$sample_id" == "sample_id" ]] && continue
    if ! bw_path "${sample_id}" > /dev/null; then
        echo "  WARNING: No bigWig found for ${sample_id} in ${BIGWIG_DIR}/"
        (( missing_bw++ ))
    fi
done < "${MANIFEST}"

if [[ $missing_bw -gt 0 ]]; then
    echo "  ${missing_bw} bigWig file(s) missing — check BIGWIG_DIR and BW_PREFIXES."
    echo "  Continuing; affected samples will be skipped."
else
    echo "  All bigWig files found."
fi

if [[ "$REGION_MODE" == "peaks" || "$REGION_MODE" == "both" ]]; then
    echo ""
    echo "=== Step 1b: Checking peak BED files ==="
    missing_peaks=0
    found_peaks=0
    while IFS=$'\t' read -r sample_id sample_name _rest; do
        [[ "$sample_id" == "sample_id" ]] && continue
        pf="$(peak_path "${sample_name}" "${sample_id}")"
        if [[ -z "$pf" ]]; then
            echo "  WARNING: No peak file for ${sample_id} (sample_name=${sample_name}) — searched ${PEAK_DIR}/${sample_name}*${PEAK_SUFFIX}"
            (( missing_peaks++ ))
        else
            (( found_peaks++ ))
        fi
    done < "${MANIFEST}"

    if [[ $found_peaks -eq 0 ]]; then
        echo ""
        echo "  ERROR: No peak BED files found in ${PEAK_DIR}/"
        echo "         Check that PEAK_DIR and PEAK_SUFFIX are correct, and that"
        echo "         03_organize_outputs.sh has been run."
        if [[ "$REGION_MODE" == "both" ]]; then
            echo "         Falling back to REGION_MODE=tss for this run."
            active_modes=("tss")
        else
            echo "         Cannot continue without peak files. Exiting."
            exit 1
        fi
    elif [[ $missing_peaks -gt 0 ]]; then
        echo "  ${missing_peaks} peak file(s) missing, ${found_peaks} found."
        echo "  Continuing; samples with missing peaks will be skipped."
    else
        echo "  All peak BED files found."
    fi
fi

echo ""
echo "=== Step 2: Building condition-level consensus peaks ==="

if [[ "$REGION_MODE" == "tss" ]]; then
    echo "  REGION_MODE=tss — skipping peak consensus generation."
fi

for tf in $TFS; do
    for condition in $CONDITIONS; do
        [[ "$REGION_MODE" == "tss" ]] && break 2
        consensus_bed="${CONSENSUS_DIR}/${tf}_${condition}.consensus.bed"

        if [[ -f "$consensus_bed" ]]; then
            echo "  Exists, skipping: $(basename "${consensus_bed}")"
            continue
        fi

        peak_files=()
        while IFS=$'\t' read -r sample_id sample_name sample_tf group; do
            [[ "$sample_id" == "sample_id" ]] && continue
            [[ "$sample_tf" != "$tf"        ]] && continue
            [[ "$group"     != "$condition" ]] && continue
            pf="$(peak_path "${sample_name}" "${sample_id}")"
            [[ -n "$pf" ]] && peak_files+=("$pf")
        done < "${MANIFEST}"

        if [[ ${#peak_files[@]} -eq 0 ]]; then
            echo "  WARNING: No peak files found for ${tf}_${condition} — skipping"
            continue
        fi

        cat "${peak_files[@]}" \
            | cut -f1-3 \
            | bedtools sort -i - \
            | bedtools merge -i - \
            > "$consensus_bed"
        echo "  Created: $(basename "${consensus_bed}") ($(wc -l < "${consensus_bed}") regions)"
    done
done

echo ""
echo "=== Step 2b: Averaging replicate bigWigs per TF × condition ==="

for tf in $TFS; do
    for condition in $CONDITIONS; do
        avg_bw="${CONSENSUS_BW_DIR}/${tf}_${condition}.avg.bigWig"

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
            [[ -f "$bw" ]] && rep_bws+=("$bw")
        done < "${MANIFEST}"

        if [[ ${#rep_bws[@]} -eq 0 ]]; then
            echo "  WARNING: No valid bigWigs for ${tf}_${condition} — skipping"
            continue
        fi

        echo "  Averaging ${#rep_bws[@]} bigWig(s) → $(basename "${avg_bw}")"
        bigwigAverage \
            --bigwigs            "${rep_bws[@]}" \
            --outFileName        "$avg_bw" \
            --numberOfProcessors "${THREADS}"
    done
done

echo ""
echo "=== Step 3a: Per-condition heatmaps ==="

for mode in "${active_modes[@]}"; do
    echo ""
    echo "  --- mode: ${mode} ---"
    for tf in $TFS; do
        for condition in $CONDITIONS; do
            tag="${tf}_${condition}.${mode}"
            matrix="${MATRIX_DIR}/${tag}.replicate.mat.gz"
            heatmap_pdf="${HEATMAP_DIR}/${tag}.replicate_heatmap.pdf"

            bw_files=()
            sample_labels=()
            rep_peak_files=()

            while IFS=$'\t' read -r sample_id sample_name sample_tf group; do
                [[ "$sample_id" == "sample_id" ]] && continue
                [[ "$sample_tf" != "$tf"        ]] && continue
                [[ "$group"     != "$condition" ]] && continue
                bw="$(bw_path "${sample_id}")"
                pf="$(peak_path "${sample_name}" "${sample_id}")"
                if [[ -f "$bw" ]]; then
                    bw_files+=("$bw")
                    sample_labels+=("$sample_id")
                fi
                [[ -n "$pf" ]] && rep_peak_files+=("$pf")
            done < "${MANIFEST}"

            if [[ ${#bw_files[@]} -eq 0 ]]; then
                echo "  WARNING: No bigWigs for ${tag} — skipping"
                continue
            fi

            if [[ "$mode" == "peaks" ]]; then
                regions_bed="${MATRIX_DIR}/${tf}_${condition}.combined_peaks.bed"
                if [[ ${#rep_peak_files[@]} -gt 0 && ! -f "$regions_bed" ]]; then
                    cat "${rep_peak_files[@]}" \
                        | cut -f1-3 \
                        | bedtools sort -i - \
                        | bedtools merge -i - \
                        > "$regions_bed"
                fi
                if [[ ! -f "$regions_bed" ]]; then
                    echo "  WARNING: No peak regions for ${tag} — skipping"
                    continue
                fi
                ref_point="center"
                title_suffix="per-replicate peaks"
            else
                regions_bed="$TSS_BED"
                ref_point="TSS"
                title_suffix="TSS ±${UPSTREAM}bp"
            fi

            if [[ ! -f "$matrix" ]]; then
                echo "  Computing matrix: ${tag}"
                computeMatrix reference-point \
                    --referencePoint          "$ref_point" \
                    --beforeRegionStartLength "${UPSTREAM}" \
                    --afterRegionStartLength  "${DOWNSTREAM}" \
                    --regionsFileName         "${regions_bed}" \
                    --scoreFileName           "${bw_files[@]}" \
                    --samplesLabel            "${sample_labels[@]}" \
                    --outFileName             "$matrix" \
                    --numberOfProcessors      "${THREADS}" \
                    --skipZeros
            fi

            if [[ ! -f "$heatmap_pdf" && -f "$matrix" ]]; then
                echo "  Plotting: ${tag}"
                plotHeatmap \
                    --matrixFile    "$matrix" \
                    --outFileName   "$heatmap_pdf" \
                    --plotTitle     "${tf} — ${condition} (${title_suffix})" \
                    --colorMap      "${COLORMAP}" \
                    --whatToShow    "plot, heatmap and colorbar" \
                    --heatmapHeight 15 \
                    --heatmapWidth  4
            fi
        done
    done
done

echo ""
echo "=== Step 3b: Cross-condition consensus heatmaps (averaged bigWigs) ==="

for mode in "${active_modes[@]}"; do
    echo ""
    echo "  --- mode: ${mode} ---"
    for tf in $TFS; do
        tag="${tf}.cross_condition_consensus.${mode}"
        matrix="${MATRIX_DIR}/${tag}.mat.gz"
        heatmap_pdf="${HEATMAP_DIR}/${tag}_heatmap.pdf"

        if [[ "$mode" == "peaks" ]]; then
            union_bed="${MATRIX_DIR}/${tf}.union_consensus.bed"
            if [[ ! -f "$union_bed" ]]; then
                cond_beds=()
                for condition in $CONDITIONS; do
                    cb="${CONSENSUS_DIR}/${tf}_${condition}.consensus.bed"
                    [[ -f "$cb" ]] && cond_beds+=("$cb")
                done
                if [[ ${#cond_beds[@]} -eq 0 ]]; then
                    echo "  WARNING: No consensus BEDs for ${tf} — skipping"
                    continue
                fi
                cat "${cond_beds[@]}" \
                    | cut -f1-3 \
                    | bedtools sort -i - \
                    | bedtools merge -i - \
                    > "$union_bed"
                echo "  Union consensus for ${tf}: $(wc -l < "${union_bed}") regions"
            fi
            regions_bed="$union_bed"
            ref_point="center"
            title_suffix="consensus peaks"
        else
            regions_bed="$TSS_BED"
            ref_point="TSS"
            title_suffix="TSS ±${UPSTREAM}bp"
        fi

        avg_bws=()
        avg_labels=()
        for condition in $CONDITIONS; do
            avg_bw="${CONSENSUS_BW_DIR}/${tf}_${condition}.avg.bigWig"
            if [[ -f "$avg_bw" ]]; then
                avg_bws+=("$avg_bw")
                avg_labels+=("${tf} ${condition}")
            fi
        done

        if [[ ${#avg_bws[@]} -lt 2 ]]; then
            echo "  Skipping ${tf}: need at least 2 conditions with averaged bigWigs (found ${#avg_bws[@]})"
            continue
        fi

        if [[ ! -f "$matrix" ]]; then
            echo "  Computing matrix: ${tag}"
            echo "    regions : ${regions_bed} ($(wc -l < "${regions_bed}") lines)"
            echo "    bigWigs : ${avg_bws[*]}"
            echo "    labels  : ${avg_labels[*]}"
            computeMatrix reference-point \
                --referencePoint          "$ref_point" \
                --beforeRegionStartLength "${UPSTREAM}" \
                --afterRegionStartLength  "${DOWNSTREAM}" \
                --regionsFileName         "${regions_bed}" \
                --scoreFileName           "${avg_bws[@]}" \
                --samplesLabel            "${avg_labels[@]}" \
                --outFileName             "$matrix" \
                --numberOfProcessors      "${THREADS}" \
                --skipZeros
            if [[ $? -ne 0 || ! -f "$matrix" ]]; then
                echo "  ERROR: computeMatrix failed for ${tag} — skipping plot"
                continue
            fi
        fi

        if [[ ! -f "$heatmap_pdf" && -f "$matrix" ]]; then
            echo "  Plotting: ${tag}"
            plotHeatmap \
                --matrixFile    "$matrix" \
                --outFileName   "$heatmap_pdf" \
                --plotTitle     "${tf} — all conditions (${title_suffix})" \
                --colorMap      "${COLORMAP}" \
                --whatToShow    "plot, heatmap and colorbar" \
                --heatmapHeight 15 \
                --heatmapWidth  4
            if [[ $? -ne 0 ]]; then
                echo "  ERROR: plotHeatmap failed for ${tag}"
            fi
        fi
    done
done

echo ""
echo "Done."
echo "  Consensus peaks : ${CONSENSUS_DIR}/"
echo "  Matrices        : ${MATRIX_DIR}/"
echo "  Heatmaps        : ${HEATMAP_DIR}/"
