# **Bisulfite Sequencing Analysis: Differential Methylation and Enrichment Studies**

## **Overview**

This repository contains a **Bisulfite Sequencing (BS-seq) analysis pipeline** designed to study **differential DNA methylation** in **breast carcinoma pre- and post-chemotherapy** samples. The workflow covers:  
\- **Methylation calling and identification of Differentially Methylated Loci (DMLs) and Regions (DMRs)**  
\- **Gene pathway enrichment analysis for DML/DMR-associated genes**  
\- **Transcription Factor (TF) enrichment analysis to identify regulatory factors involved in methylation changes**

## **Dataset Information**

* **Dataset:** Xmal-RRBS breast carcinoma bisulfite sequencing  
* **NCBI GEO Accession:** `GSE202402`  
* **Reference Genome:** *Homo sapiens* (hg38)  
* **Input Data:** Raw sequencing data from **pre- and post-chemotherapy tumor samples**


## **Workflow Overview**

The analysis is divided into **three parts**:

### **1. Jupyter Notebook: Data Processing & Methylation Calling (`BisulfiteSeq_Preprocessing.ipynb`)**

**Preprocessing Steps**

* Download bisulfite sequencing data (FASTQ)  
* Perform **Quality Control (QC)** with **FastQC**  
* Trim low-quality reads using **Trim Galore**

**Read Alignment & Methylation Calling**

* Align reads to **hg38 genome** using **Bismark**  
* Remove duplicate reads with **Samtools & Picard**  
* Perform **methylation calling** and generate **DSS input files** for DML/DMR analysis

**Output:**

* Filtered **CpG methylation matrices**  
* `dss_input.csv` for further analysis in R


### **2. R Script: DML/DMR Analysis with Pathway Enrichment (`BisulfiteSequencing_DML_DMR_Pathway.R`)**

**DML & DMR Identification**

* Load methylation data from `dss_input.csv`  
* Perform **Differentially Methylated Loci (DML) analysis** using **DSS package**  
* Identify **Differentially Methylated Regions (DMRs)**

**Functional Annotation & Enrichment**

* Map **DML/DMR sites to genes** using **ChIPseeker**  
* Perform **GO (Gene Ontology) enrichment** for **biological processes (BP) & molecular functions (MF)**  
* Perform **KEGG pathway enrichment** to identify **cancer-associated pathways**

**Output:**

* `significant_DMLs.csv` (DML sites)  
* `significant_DMRs.csv` (DMR regions)  
* `GO_Enrichment_BP_DMRs.csv`, `KEGG_Enrichment_DMRs.csv` (Functional enrichment results)  
* **Visualizations:** Volcano plots, Heatmaps


### **3. R Script: Transcription Factor (TF) Enrichment Analysis (`BisulfiteSeq_DML_DMR_TF.R`)**

**TF Enrichment Workflow**

* Load **annotated DML/DMR sites**  
* Retrieve **known TF motifs from JASPAR2022**  
* Scan for **TF motifs in DMLs & DMRs** using **motifmatchr**  
* Perform **background correction & statistical testing**

**Output:**

* `TF_enrichment_rawscore.csv` (Significant TF motifs based on raw scores)  
* `TF_enrichment_FDR.csv` (TF motifs corrected for multiple testing)  
* **TF enrichment visualizations (bar plots)**


## **Repository Structure**

`Bisulfite_Sequencing_Analysis/`

`│── data/                                  # Raw and processed methylation data`
`│── results/                               # Output files, figures, and reports`
`│── reference/                			 # Reference genome and annotations`
`│── BisulfiteSeq_Preprocessing.ipynb    # Jupyter notebook for preprocessing & methylation calling`
`│── BisulfiteSequencing_DML_DMR_Pathway.R  # DML/DMR pathway analysis`
`│── BisulfiteSeq_DML_DMR_TF.R  # TF enrichment analysis`
`│── README.md               # Project documentation`
`│── requirements.txt        # Dependencies and installation requirements`


## **Running the Workflow**

### **1. Run Jupyter Notebook (Preprocessing & Methylation Calling)**

To execute the **data acquisition, QC, alignment, and methylation calling steps**, open the notebook:

```bash
jupyter notebook bisulfite_sequencing_workflow.ipynb
```

### **2. Run DML/DMR Analysis in R**

To perform **DML/DMR identification and pathway enrichment**, execute:

``` r
source("BisulfiteSequencing_DML_DMR_Pathway_Enrichment.R")
```

### **3. Run TF Enrichment Analysis in R**

To scan **TF motifs in DML/DMRs and identify enriched TFs**, execute:

``` r
source("BisulfiteSequencing_DML_DMR_TF_Enrichment.R")
```

---

## **Dependencies & Installation**

### **Conda Environment Setup**

``` bash
conda create --name bisulfite_analysis_env python=3.9
conda activate bisulfite_analysis_env
pip install -r requirements.txt
```

### **Bioinformatics Tools (Conda)**

```bash
conda install -y -c bioconda fastqc trim-galore bismark samtools picard
```

### **R Package Installation**

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DSS", "bsseq", "ggplot2", "pheatmap", "clusterProfiler", "ChIPseeker",

                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "JASPAR2022", 

                       "TFBSTools", "motifmatchr", "regioneR", "Biostrings"))
```

## **Contributors**

This project is part of an effort to explore **epigenetic changes in cancer** using **bisulfite sequencing** and **computational genomics**.

**Feel free to contribute or collaborate\!**


## **Acknowledgments**

* This project utilizes **publicly available bisulfite sequencing datasets** from **NCBI GEO (GSE202402)**.

