# TFF1

This repository contains five independent analysis projects focused on regulatory T cells (Tregs) and thymic epithelial cells (TECs) in mouse models. Projects span single-cell RNA-seq (scRNA-seq), single-cell TCR-seq (scTCR-seq), bulk ATAC-seq, and CUT&RUN data. Each project is self-contained with its own notebooks, analysis steps, and conda environments.

---

## Repository Structure

| Folder | Project |
|--------|---------|
| [`0_pan_tissue_cd4s/`](0_pan_tissue_cd4s/) | Pan-tissue CD4 T cell analysis, focused on Tregs |
| [`1_wt_ko_tecs/`](1_wt_ko_tecs/) | WT vs TFF1 knockout thymic epithelial cell (TEC) analysis |
| [`2_treg_sltbi/`](2_treg_sltbi/) | Treg analysis following sublethal total body irradiation (SLTBI) |
| [`3_bulk_atacseq/`](3_bulk_atacseq/) | Bulk ATAC-seq differential accessibility and motif enrichment (TEC organoids: Tff1 OE vs GFP EV) |
| [`4_cutnrun/`](4_cutnrun/) | CUT&RUN differential binding, motif analysis, and visualization (TEC organoids: Tff1 OE vs GFP EV) |

---

## Projects

### 1. Pan-Tissue CD4 T Cell Analysis

**Dataset:** Mouse scRNA-seq of CD4 T cells isolated from multiple tissues (bone marrow, lung, skin, spleen, and thymus). CD4 T cells are filtered down to regulatory T cells (Tregs) for downstream analysis.

#### Python Notebooks (`1_pyzone/0_notebooks/0_scRNA/`)

| Step | Notebook | Description |
|------|----------|-------------|
| QC | [`0_qc/0_qc.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/0_qc/0_qc.ipynb) | Quality control filtering of raw scRNA-seq data |
| Annotation | [`1_annotation/0_bm.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/0_bm.ipynb) | Cell type annotation — bone marrow |
| Annotation | [`1_annotation/2_lung.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/2_lung.ipynb) | Cell type annotation — lung |
| Annotation | [`1_annotation/3_skin.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/3_skin.ipynb) | Cell type annotation — skin |
| Annotation | [`1_annotation/4_spleen.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/4_spleen.ipynb) | Cell type annotation — spleen |
| Annotation | [`1_annotation/5_thymus.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/5_thymus.ipynb) | Cell type annotation — thymus |
| Annotation | [`1_annotation/6_combined_tissue.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/6_combined_tissue.ipynb) | Combined cross-tissue annotation |
| Annotation | [`1_annotation/7_tregs.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/1_annotation/7_tregs.ipynb) | Treg-specific annotation and subclustering |
| DEG | [`2_deg/0_tissue_tregs.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/2_deg/0_tissue_tregs.ipynb) | Differentially expressed gene (DEG) analysis across tissues |
| Figures | [`3_figures/0_treg_fig1.ipynb`](0_pan_tissue_cd4s/1_pyzone/0_notebooks/0_scRNA/3_figures/0_treg_fig1.ipynb) | Figure generation for Treg results |

#### R Scripts (`2_rzone/`)

Volcano plots of DEG results, one per tissue:

| Script | Tissue |
|--------|--------|
| [`0_thymus_tregs_volcano_plot.R`](0_pan_tissue_cd4s/2_rzone/0_thymus_tregs_volcano_plot.R) | Thymus |
| [`1_spleen_tregs_volcano_plot.R`](0_pan_tissue_cd4s/2_rzone/1_spleen_tregs_volcano_plot.R) | Spleen |
| [`2_bm_tregs_volcano_plot.R`](0_pan_tissue_cd4s/2_rzone/2_bm_tregs_volcano_plot.R) | Bone marrow |
| [`3_liver_tregs_volcano_plot.R`](0_pan_tissue_cd4s/2_rzone/3_liver_tregs_volcano_plot.R) | Liver |
| [`4_skin_tregs_volcano_plot.R`](0_pan_tissue_cd4s/2_rzone/4_skin_tregs_volcano_plot.R) | Skin |

> **Note:** R package dependencies for these scripts are not yet captured in a conda environment file.

#### Environment
```bash
conda env create -f 0_pan_tissue_cd4s/envs/scrna.yaml
conda activate scrna
```

---

### 2. WT vs TFF1 Knockout Thymic Epithelial Cell Analysis

**Dataset:** Mouse scRNA-seq comparing thymic epithelial cells (TECs) from wild-type (WT) and TFF1 knockout (KO) animals. The analysis characterizes how loss of TFF1 affects TEC composition and transcriptional state, including stromal populations.

#### Python Notebooks (`1_pyzone/0_notebooks/0_scRNA/`)

| Step | Notebook | Description |
|------|----------|-------------|
| QC | [`0_qc/0_qc.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/0_qc/0_qc.ipynb) | Quality control filtering of raw scRNA-seq data |
| QC | [`0_qc/1_broad_annotations.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/0_qc/1_broad_annotations.ipynb) | Broad cell type annotation to identify major populations |
| Annotation | [`1_annotation/0_ec.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/1_annotation/0_ec.ipynb) | Endothelial cell (EC) subclustering and annotation |
| Annotation | [`1_annotation/1_fb.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/1_annotation/1_fb.ipynb) | Fibroblast (FB) subclustering and annotation |
| Annotation | [`1_annotation/2_tec.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/1_annotation/2_tec.ipynb) | TEC subclustering and annotation |
| Annotation | [`1_annotation/3_combine_cell_subsets.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/1_annotation/3_combine_cell_subsets.ipynb) | Integration of all annotated cell subsets |
| DEG | [`2_deg/0_tec_degs.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/2_deg/0_tec_degs.ipynb) | Differentially expressed gene analysis between WT and KO TECs |
| MiloR | [`3_miloR/tecs_miloR.ipynb`](1_wt_ko_tecs/1_pyzone/0_notebooks/0_scRNA/3_miloR/tecs_miloR.ipynb) | Differential abundance analysis using MiloR |

#### Environment
```bash
conda env create -f 1_wt_ko_tecs/envs/scrna.yaml
conda activate scrna
```

---

### 3. Treg Analysis Post Sublethal Total Body Irradiation (SLTBI)

**Dataset:** Mouse Tregs profiled following sublethal total body irradiation (SLTBI). This project includes matched scRNA-seq and scTCR-seq data, enabling simultaneous transcriptional and clonotype-level analysis of Tregs.

#### Analysis order

> Run demultiplexing (R) first, then proceed with the Python scRNA-seq and scTCR-seq notebooks.

#### R Notebook — Demultiplexing (`1_rzone/0_notebooks/`)

| Notebook | Description |
|----------|-------------|
| [`0_demultiplex.ipynb`](2_treg_sltbi/1_rzone/0_notebooks/0_demultiplex.ipynb) | Sample demultiplexing using Seurat/HTODemux |

```bash
conda env create -f 2_treg_sltbi/envs/seurat.yaml
conda activate seurat
```

#### Python Notebooks — scRNA-seq (`0_pyzone/0_scRNA/0_notebooks/`)

| Step | Notebook | Description |
|------|----------|-------------|
| QC | [`1_qc.ipynb`](2_treg_sltbi/0_pyzone/0_scRNA/0_notebooks/1_qc.ipynb) | Quality control filtering of scRNA-seq data |
| Analysis | [`2_treg_analysis.ipynb`](2_treg_sltbi/0_pyzone/0_scRNA/0_notebooks/2_treg_analysis.ipynb) | Treg clustering, annotation, and transcriptional analysis |
| DEG | [`3_degs.ipynb`](2_treg_sltbi/0_pyzone/0_scRNA/0_notebooks/3_degs.ipynb) | Differentially expressed gene analysis |

```bash
conda env create -f 2_treg_sltbi/envs/scrna.yaml
conda activate scrna
```

#### Python Notebook — scTCR-seq (`0_pyzone/1_scTCR/0_notebooks/`)

| Notebook | Description |
|----------|-------------|
| [`0_tcr.ipynb`](2_treg_sltbi/0_pyzone/1_scTCR/0_notebooks/0_tcr.ipynb) | T cell receptor clonotype analysis using matched scTCR-seq data |

```bash
conda env create -f 2_treg_sltbi/envs/sctcr.yaml
conda activate sctcr
```

---

### 4. Bulk ATAC-seq Differential Accessibility Analysis

**Dataset:** Bulk ATAC-seq from mouse thymic epithelial cell (TEC) organoids comparing Tff1 overexpression (OE) to GFP empty vector (EV) control. Raw reads are processed with the nf-core/atacseq v2.1.2 Nextflow pipeline, followed by differential accessibility analysis, heatmap visualization, and transcription factor motif enrichment.

#### Pipeline Scripts (`0_analysis/nextflow_scripts/`)

| File | Description |
|------|-------------|
| [`run_atacseq.sh`](3_bulk_atacseq/0_analysis/nextflow_scripts/run_atacseq.sh) | SLURM submission script for nf-core/atacseq v2.1.2 |
| [`params.yml`](3_bulk_atacseq/0_analysis/nextflow_scripts/params.yml) | Pipeline parameters |
| [`nextflow.config`](3_bulk_atacseq/0_analysis/nextflow_scripts/nextflow.config) | Nextflow configuration |
| [`samplesheet.csv`](3_bulk_atacseq/0_analysis/nextflow_scripts/samplesheet.csv) | Input sample sheet |
| [`run_global_heatmap.sh`](3_bulk_atacseq/0_analysis/nextflow_scripts/run_global_heatmap.sh) | Global signal heatmap generation |
| [`run_global_heatmap_per_replicate.sh`](3_bulk_atacseq/0_analysis/nextflow_scripts/run_global_heatmap_per_replicate.sh) | Per-replicate heatmap generation |

#### Analysis Notebooks (`0_analysis/`)

| Step | Notebook | Language | Description |
|------|----------|----------|-------------|
| DiffBind | [`0_diffbind/diffbind_atac.ipynb`](3_bulk_atacseq/0_analysis/0_diffbind/diffbind_atac.ipynb) | R | Differential accessibility analysis using DiffBind |
| Heatmaps | [`1_stacked_heatmaps/stacked_heatmap.ipynb`](3_bulk_atacseq/0_analysis/1_stacked_heatmaps/stacked_heatmap.ipynb) | Python | Stacked ATAC-seq signal heatmaps |
| Heatmaps | [`1_stacked_heatmaps/replicate_heatmaps.ipynb`](3_bulk_atacseq/0_analysis/1_stacked_heatmaps/replicate_heatmaps.ipynb) | Python | Per-replicate ATAC-seq heatmaps |
| chromVAR | [`2_chromvar/chromvar_motif_enrichment.ipynb`](3_bulk_atacseq/0_analysis/2_chromvar/chromvar_motif_enrichment.ipynb) | R | Transcription factor motif enrichment using chromVAR |

#### Reference

- [`reference/mm10-blacklist.v2.bed`](3_bulk_atacseq/reference/mm10-blacklist.v2.bed) — ENCODE mm10 blacklist regions

#### Environment
```bash
conda env create -f 3_bulk_atacseq/envs/bulkatac.yaml
conda activate bulkatac
```

---

### 5. CUT&RUN Differential Binding Analysis

**Dataset:** CUT&RUN profiling of transcription factors (FosB, JunB) and histone mark H3K4me3 in mouse thymic epithelial cell (TEC) organoids comparing Tff1 overexpression (OE) to GFP empty vector (EV) control. Raw reads are processed with the nf-core/cutandrun Nextflow pipeline with SEACR peak calling, followed by DiffBind differential binding analysis, HOMER motif enrichment, deepTools heatmaps, and track-level visualization.

#### Pipeline Scripts (`0_analysis/scripts/`)

| File | Description |
|------|-------------|
| [`01_setup.sh`](4_cutnrun/0_analysis/scripts/01_setup.sh) | Environment and directory setup |
| [`02_run_pipeline.sh`](4_cutnrun/0_analysis/scripts/02_run_pipeline.sh) | SLURM submission script for nf-core/cutandrun |
| [`04_run_homer.sh`](4_cutnrun/0_analysis/scripts/04_run_homer.sh) | HOMER known motif analysis |
| [`05_deeptools_heatmaps.sh`](4_cutnrun/0_analysis/scripts/05_deeptools_heatmaps.sh) | deepTools signal heatmap generation |
| [`06_avg_bigwigs.sh`](4_cutnrun/0_analysis/scripts/06_avg_bigwigs.sh) | Average bigWig track generation across replicates |
| [`params.yaml`](4_cutnrun/0_analysis/scripts/params.yaml) | Pipeline parameters |

#### Analysis Notebooks (`0_analysis/`)

| Step | Notebook | Language | Description |
|------|----------|----------|-------------|
| DiffBind | [`0_diffbind/fosb_diffbind_seacr.ipynb`](4_cutnrun/0_analysis/0_diffbind/fosb_diffbind_seacr.ipynb) | R | Differential binding analysis — FosB |
| DiffBind | [`0_diffbind/junb_diffbind_seacr.ipynb`](4_cutnrun/0_analysis/0_diffbind/junb_diffbind_seacr.ipynb) | R | Differential binding analysis — JunB |
| DiffBind | [`0_diffbind/h3k4me3_diffbind_seacr.ipynb`](4_cutnrun/0_analysis/0_diffbind/h3k4me3_diffbind_seacr.ipynb) | R | Differential binding analysis — H3K4me3 |
| HOMER | [`1_homer/known_motif_analysis.ipynb`](4_cutnrun/0_analysis/1_homer/known_motif_analysis.ipynb) | Python | Known transcription factor motif enrichment in H3K4me3 peaks |
| Heatmaps | [`2_heatmaps/stacked_heatmap.ipynb`](4_cutnrun/0_analysis/2_heatmaps/stacked_heatmap.ipynb) | Python | Stacked CUT&RUN signal heatmaps |
| Trackplots | [`3_trackplots/trackplots.ipynb`](4_cutnrun/0_analysis/3_trackplots/trackplots.ipynb) | Python | Genome browser-style track visualization |

#### Reference

- [`references/blacklist/mm10-blacklist.v2.bed.gz`](4_cutnrun/0_analysis/references/blacklist/mm10-blacklist.v2.bed.gz) — ENCODE mm10 blacklist regions (gzipped)

#### Environment
```bash
conda env create -f 4_cutnrun/envs/cutnrun.yaml
conda activate cutnrun
```

---

## Getting Started

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- [Jupyter](https://jupyter.org/) (included in the conda environments)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/Marcel-van-den-Brink-Lab/tff1.git
   cd tff1
   ```

2. Create and activate the appropriate environment for the project you want to run (see environment instructions under each project above).

3. Launch Jupyter and open the notebooks in the order listed for each project. 

> Data files are not included in this repository. Please check GEO for input data.
