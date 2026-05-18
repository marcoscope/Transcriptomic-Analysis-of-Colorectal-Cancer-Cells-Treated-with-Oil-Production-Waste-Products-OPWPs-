# Transcriptomic-Analysis-of-Colorectal-Cancer-Cells-Treated-with-Oil-Production-Waste-Products-OPWPs-
An end-to-end RNA-seq pipeline replicating transcriptomic changes in colorectal cancer cells treated with olive mill waste extracts (OPWPs).
# End-to-End RNA-seq Pipeline: Transcriptomic Profiling of Colorectal Cancer Cells Under OPWP Treatment

## 📌 Project Overview
This project presents a complete, high-throughput RNA-sequencing (RNA-seq) data analysis pipeline executed within our bioinformatics laboratory. The study investigates the cellular and molecular responses of colorectal cancer cell lines (HCT116 and LoVo) treated with **Oil Production Waste Products (OPWPs)** and its primary polyphenolic component, **Hydroxytyrosol**.

The pipeline encompasses raw read quality control, adapter trimming, splice-aware genome alignment, transcript quantification, and advanced differential expression and functional enrichment using R. Our lab-bench pipeline successfully uncovers significant mitochondrial alterations and transcriptomic shifts induced by the treatment.

---

## 🛠️ Laboratory Workflow & Pipeline Steps

### Phase 1: Galaxy-Based Upstream Processing
Our upstream analysis was fully structured and processed using standard tools integrated within the Galaxy interface:

1. **Raw Data Quality Control (`FastQC`)**
   * Performed initial quality checks on raw paired-end reads to assess base quality scores, GC content, and adapter contamination.
2. **Adapter & Quality Trimming (`Cutadapt`)**
   * Processed reads to remove Illumina TruSeq universal adapters. 
   * Strict parameters were enforced, resulting in a conservative ~5bp clean trimming (reducing read length from 151bp to ~146bp) to preserve maximum mapping integrity.
3. **Post-Trimming Validation (`MultiQC`)**
   * Aggregated all individual FastQC metrics into a single report to verify successful adapter removal and a sharp increase in overall base qualities.
4. **Splice-Aware Alignment (`STAR Aligner`)**
   * Mapped the clean reads against the Human Reference Genome (`GRCh38`). STAR was chosen to handle splicing junctions accurately across the exons.
5. **Transcript Quantification (`featureCounts`)**
   * Assigned mapped reads to genomic features to generate a clean, raw count matrix. The optimization achieved highly efficient assigned read rates (>60-70%), creating a reliable matrix for downstream R processing.

### Phase 2: R-Based Downstream Analysis & Enrichment
Once the raw count matrix was generated, the analysis shifted to an R environment for statistical processing and advanced data visualization:

6. **Differential Gene Expression (DGE) with `DESeq2`**
   * Normalization of raw counts and calculation of Log2 Fold Changes to identify significantly altered genes under treatment conditions.
7. **Gene Set Enrichment Analysis (GSEA)**
   * Performed functional annotation based on the MSigDB Hallmark gene sets to extract core biological pathways significantly affected.

---

## 📊 Key Laboratory Findings & Visualizations

Our pipeline successfully captured strong transcriptomic signals showing that OPWP treatment induces major metabolic and survival stress in cancer cells:

### 1. Upregulation of Oxidative Phosphorylation
GSEA analysis revealed a highly significant upregulation of the **Oxidative Phosphorylation** pathway ($NES = 2.05$, $FDR = 1.56 \times 10^{-4}$), demonstrating a major stimulation of mitochondrial respiratory chain functions.

### 2. Global Functional Enrichment (Top Pathways)
Using specialized `ggplot2` visualizations, we captured the top 12 pathways for both Upregulated and Downregulated genes, shedding light on extensive **DNA Damage Repair** responses (such as P53 Signaling and Mismatch Repair) activated by the treatment.

---

## 📁 Repository Structure
* `📁 data/` : Contains the structured raw count matrix generated from featureCounts.
* `📁 scripts/` : Contains the complete, reproducible R scripts used for plotting and analysis.
* `📁 results/` : Contains high-resolution plots (Dot plots and Enrichment plots).

