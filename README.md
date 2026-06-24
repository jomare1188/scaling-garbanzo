# RNAseq analysis of *Fusarium verticillioides*

## Overview

This repository accompanies the study of *Fusarium verticillioides* and its interaction with a virus


## Palette of colours

"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF"



## RNAseq Workflow Description

| sample      | fastq_1 | fastq_2 | group | strandedness | rep |
|------------|---------|---------|--------|--------------|-----|
| Fv_IAC_VF_1 | /dados01/samuel/israel/BCL/1_Marcio_L1-ds.6e0b458ce2384e8cb21357ec3add4b18/1_Marcio_S25_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/1_Marcio_L1-ds.6e0b458ce2384e8cb21357ec3add4b18/1_Marcio_S25_L001_R2_001.fastq.gz | Fv_IAC_VF | auto | 1 |
| Fv_IAC_VF_2 | /dados01/samuel/israel/BCL/2_Marcio_L1-ds.a5021cd52e334b1783ff3660db27a95b/2_Marcio_S26_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/2_Marcio_L1-ds.a5021cd52e334b1783ff3660db27a95b/2_Marcio_S26_L001_R2_001.fastq.gz | Fv_IAC_VF | auto | 2 |
| Fv_IAC_VF_3 | /dados01/samuel/israel/BCL/3_Marcio_L2-ds.aba760a0b0c64830ba843995dc4ca4d8/3_Marcio_S55_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/3_Marcio_L2-ds.aba760a0b0c64830ba843995dc4ca4d8/3_Marcio_S55_L002_R2_001.fastq.gz | Fv_IAC_VF | auto | 3 |
| Fv_IAC_VF_4 | /dados01/samuel/israel/BCL/4_Marcio_L2-ds.226c81a392dd4e33b68b4f663e8fed38/4_Marcio_S56_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/4_Marcio_L2-ds.226c81a392dd4e33b68b4f663e8fed38/4_Marcio_S56_L002_R2_001.fastq.gz | Fv_IAC_VF | auto | 4 |
| Fv_IAC_VI_1 | /dados01/samuel/israel/BCL/5_Marcio_L1-ds.15178dbbc2574ba8b4897a507ee4ddd7/5_Marcio_S27_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/5_Marcio_L1-ds.15178dbbc2574ba8b4897a507ee4ddd7/5_Marcio_S27_L001_R2_001.fastq.gz | Fv_IAC_VI | auto | 1 |
| Fv_IAC_VI_2 | /dados01/samuel/israel/BCL/6_Marcio_L1-ds.10e33f2626ce4884b5b266f3026a9d90/6_Marcio_S28_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/6_Marcio_L1-ds.10e33f2626ce4884b5b266f3026a9d90/6_Marcio_S28_L001_R2_001.fastq.gz | Fv_IAC_VI | auto | 2 |
| Fv_IAC_VI_3 | /dados01/samuel/israel/BCL/7_Marcio_L2-ds.f194cc72af0f4b14aa2ddc55e6c36c28/7_Marcio_S57_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/7_Marcio_L2-ds.f194cc72af0f4b14aa2ddc55e6c36c28/7_Marcio_S57_L002_R2_001.fastq.gz | Fv_IAC_VI | auto | 3 |
| Fv_IAC_VI_4 | /dados01/samuel/israel/BCL/8_Marcio_L2-ds.2677da67026447f7b907796488a58de7/8_Marcio_S58_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/8_Marcio_L2-ds.2677da67026447f7b907796488a58de7/8_Marcio_S58_L002_R2_001.fastq.gz | Fv_IAC_VI | auto | 4 |


### 1. **References**

- `Genome Assembly`: `/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.fa.gz`
- `Proteins`: `/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_protein.faa`
- `GTF`: `/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz`

### 2. **Protein Annotation**

We used `emapper-2.1.3` from `EggNOG v5.0` to get annotations (included GO annotations) for the proteins of the genome based on orthology relationships. 
- Results: `/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/eggnog_anot/eggnog_anot.emapper.annotations`

### 3. **RNAseq processing**

We used a `Nextflow v25.04.7` pipeline `rnaseq (v3.12.0)` from nf-core (https://nf-co.re/rnaseq/3.12.0) to preprocces, align and quantify RNAseq data

We used the default method from `rnaseq (v3.12.0)` which uses `STAR` aligner and `Salmon` to quantify transcript abundance.

Full report of preprocess and aligment can be found in 
`/dados02/jorge/israel_rnaseq/rnaseq/run01/multiqc/star_salmon/multiqc_report.html`

generated reads: 449.2M
aligned reads to reference: 272.2M

### 4. **Exploratory Analysis**

for deseq2 we used vst (variance stabilization transformation) and to filter low count genes:
 - min_count <- 5
 - min_samples <- 5
 - lfc_threshold <- 1
 - padj_threshold <- 0.05

- Principal component analysis: We load the quantification data produced by Salmon into DESeq2 and used the transformed counts.

   - PCA DESeq2
![PCA_DESeq2](rnaseq/run01/star_salmon/deseq2_qc/all_pca.png)

  - Initial genes: 16290
  - Genes after filtering: 9330
  - Genes removed: 6960 (42.7%)


### 5. **Differential Expression Analysis (DEA)**

We conducted a differential expression analysis (DEA) on this contrasts filtering for p-adj (FDR) < 0.05 and |log 2 fold change| > 1

Fv_IAC_VI_vs_Fv_IAC_VF

for example downregulated genes in Fv_IAC_VI_vs_Fv_IAC_VF mean genes differentially less expressed in Fv_IAC_VI than in Fv_IAC_VF and consecuentely upregulated genes means more expressed in Fv_IAC_VF than in Fv_IAC_VI



- results DESeq2: `/dados02/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc`


### Summary of DEG Counts (|LFC| > 1, padj < 0.05)

| Contrast | Description | Upregulated | Downregulated | Total_DEGs |
|-----------|-------------|------------:|--------------:|-----------:|
| Fv_IAC_VI_vs_Fv_IAC_VF | Fv_IAC_VI vs Fv_IAC_VF | 931 | 1336 | 2267 |



### 6. **Functional Enrichment Analysis**

To get insights about the function and the processes that are represented by the sets of up-regulated and down-regulated genes we carried out over representation analysis (ORA) for gene ontology terms (GO)

- GO: We used topGO R package (v2.58.0), p-value < 0.05 and corrected for multiple testing using BH (FDR) procedure

- Up genes
![up_enrichment](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_upregulated.png)
![up_table](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_upregulated.csv)

- Down genes
![down_enrichment](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_downregulated.png)
![down_table](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_downregulated.csv)


See erichment results for all contrastas in the server

- `/dados02/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment`









