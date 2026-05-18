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

### Phase 2: R-Based Downstream Analysis
Once the raw count matrix was generated, the analysis shifted to an R environment for statistical processing and advanced data visualization:

6. **Differential Gene Expression (DGE) with `DESeq2`**
   * Normalization of raw counts and calculation of Log2 Fold Changes to identify significantly altered genes under treatment conditions.
     <img width="1041" height="547" alt="Volcano plot" src="https://github.com/user-attachments/assets/76cb9bb5-e4a9-4a0e-ac52-cbd86d1fd10c" />
  **Data Visualization via Hierarchical Clustering `pheatmap`**
* We performed hierarchical clustering on the top 10 differentially expressed genes (DEGs) to visually confirm clear and consistent Z-score expression shifts between the OPWP-treated group and the control samples.
  <img width="1152" height="783" alt="Heatmap" src="https://github.com/user-attachments/assets/406fc654-0239-46d2-8533-39fa5a9dfe62" />
7. **Gene Set Enrichment Analysis (GSEA)**
   * Performed functional enrichment analysis using GSEA based on the MSigDB Hallmark gene sets. This approach allowed us to profile global biological pathways and metabolic shifts comprehensively without relying on arbitrary p-value cutoffs, ensuring a robust representation of the cellular changes induced by the treatment. *(Source code available in the scripts folder)*

---

## 📊 Laboratory Visualizations (Generated using R)

All downstream graphical plots and biological profiling visualizations were programmatically generated within the R environment utilizing `ggplot2` and specialized enrichment plotting packages.

### 1. Core Pathway Analysis (GSEA Enrichment Plot)
Gene Set Enrichment Analysis (GSEA) revealed a highly significant upregulation of the **Oxidative Phosphorylation** pathway ($NES = 2.05655$, $FDR = 1.56 \times 10^{-4}$), demonstrating a major stimulation of mitochondrial respiratory chain functions under OPWP treatment.

<details>
<summary><b>👁️ View GSEA Enrichment Plot</b></summary>
<img width="2377" height="2557" alt="لا" src="https://github.com/user-attachments/assets/00db32ec-7e33-4918-be45-a9ae6eea94f0" />


##### **Visual Result:**
> **Note:** Replace `results/gsea_oxidative.png` with your actual uploaded GSEA plot path in the repository.
![GSEA Hallmark Oxidative Phosphorylation](results/gsea_oxidative.png)
</details>

---

### 2. Global Pathway Enrichment (Top 12 Pathways via `ggplot2`)
Advanced horizontal `ggplot2` dot plots were generated to visualize the top 12 Upregulated and Downregulated hallmark pathways based on their Normalized Enrichment Scores (NES) and FDR q-values.

#### 🟢 Top 12 Upregulated Pathways
<details>
<summary><b>📊 View Upregulated Dot Plot</b></summary>
<img width="1041" height="547" alt="Rplot04" src="https://github.com/user-attachments/assets/882d93ea-6d3c-495f-84d0-2a6d049222a7" />
##### **Visual Result:**
> **Note:** Replace `results/upregulated_plot.png` with your actual image path.
![Top 12 Upregulated Pathways](results/upregulated_plot.png)
</details>


#### 🔵 Top 12 Downregulated Pathways
<details>
<summary><b>📊 View Downregulated Dot Plot</b></summary>
<img width="1041" height="547" alt="Rplot03" src="https://github.com/user-attachments/assets/c7711964-5b16-40a2-bf55-d6dc57947446" />
##### **Visual Result:**
> **Note:** Replace `results/downregulated_plot.png` with your actual image path.
![Top 12 Downregulated Pathways](results/downregulated_plot.png)
</details>


---

## 📁 Repository Structure
* `📁 data/` : Contains the structured raw count matrix generated from featureCounts.
* `📁 scripts/` : Contains the complete, reproducible R scripts used for plotting and analysis.
* `📁 results/` : Contains high-resolution plots (Dot plots and Enrichment plots).

