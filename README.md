# RNAseq analysis of *Fusarium verticillioides*

## Overview

This repository accompanies the study of *Fusarium verticillioides* and its interaction with a mycovirus.

Two strains of the same fungal background were compared in axenic culture. The DESeq2 `group`
factor has exactly two levels:

| `group` level | Mycovirus |
|---|---|
| `FvIAC-VF` | absent |
| `FvIAC-HT` | present |



**Contrast direction (important).** All differential expression is reported for the contrast

```r
results(dds, contrast = c("group", "FvIAC-HT", "FvIAC-VF"))   # numerator, denominator
```

Therefore:

- **Upregulated** = expressed **more** in `FvIAC-HT` than in `FvIAC-VF`.
- **Downregulated** = expressed **less** in `FvIAC-HT` than in `FvIAC-VF`.
- Positive log2FoldChange / positive NES = higher in HT strain.


## Palette of colours

`#E64B35FF` for `FvIAC-HT`

`#4DBBD5FF` for `FvIAC-VF`

Other colours that can be used: `#00A087FF`, `#3C5488FF`, `#F39B7FFF`

## RNAseq Workflow Description

| sample | fastq_1 | fastq_2 | group | strandedness | rep |
|---|---|---|---|---|---|
| Fv_IAC_VF_1 | /dados01/samuel/israel/BCL/1_Marcio_L1-ds.6e0b458ce2384e8cb21357ec3add4b18/1_Marcio_S25_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/1_Marcio_L1-ds.6e0b458ce2384e8cb21357ec3add4b18/1_Marcio_S25_L001_R2_001.fastq.gz | FvIAC-VF | auto | 1 |
| Fv_IAC_VF_2 | /dados01/samuel/israel/BCL/2_Marcio_L1-ds.a5021cd52e334b1783ff3660db27a95b/2_Marcio_S26_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/2_Marcio_L1-ds.a5021cd52e334b1783ff3660db27a95b/2_Marcio_S26_L001_R2_001.fastq.gz | FvIAC-VF | auto | 2 |
| Fv_IAC_VF_3 | /dados01/samuel/israel/BCL/3_Marcio_L2-ds.aba760a0b0c64830ba843995dc4ca4d8/3_Marcio_S55_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/3_Marcio_L2-ds.aba760a0b0c64830ba843995dc4ca4d8/3_Marcio_S55_L002_R2_001.fastq.gz | FvIAC-VF | auto | 3 |
| Fv_IAC_VF_4 | /dados01/samuel/israel/BCL/4_Marcio_L2-ds.226c81a392dd4e33b68b4f663e8fed38/4_Marcio_S56_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/4_Marcio_L2-ds.226c81a392dd4e33b68b4f663e8fed38/4_Marcio_S56_L002_R2_001.fastq.gz | FvIAC-VF | auto | 4 |
| Fv_IAC_VI_1 | /dados01/samuel/israel/BCL/5_Marcio_L1-ds.15178dbbc2574ba8b4897a507ee4ddd7/5_Marcio_S27_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/5_Marcio_L1-ds.15178dbbc2574ba8b4897a507ee4ddd7/5_Marcio_S27_L001_R2_001.fastq.gz | FvIAC-HT | auto | 1 |
| Fv_IAC_VI_2 | /dados01/samuel/israel/BCL/6_Marcio_L1-ds.10e33f2626ce4884b5b266f3026a9d90/6_Marcio_S28_L001_R1_001.fastq.gz | /dados01/samuel/israel/BCL/6_Marcio_L1-ds.10e33f2626ce4884b5b266f3026a9d90/6_Marcio_S28_L001_R2_001.fastq.gz | FvIAC-HT | auto | 2 |
| Fv_IAC_VI_3 | /dados01/samuel/israel/BCL/7_Marcio_L2-ds.f194cc72af0f4b14aa2ddc55e6c36c28/7_Marcio_S57_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/7_Marcio_L2-ds.f194cc72af0f4b14aa2ddc55e6c36c28/7_Marcio_S57_L002_R2_001.fastq.gz | FvIAC-HT | auto | 3 |
| Fv_IAC_VI_4 | /dados01/samuel/israel/BCL/8_Marcio_L2-ds.2677da67026447f7b907796488a58de7/8_Marcio_S58_L002_R1_001.fastq.gz | /dados01/samuel/israel/BCL/8_Marcio_L2-ds.2677da67026447f7b907796488a58de7/8_Marcio_S58_L002_R2_001.fastq.gz | FvIAC-HT | auto | 4 |

**Design:** 2 groups × **4 biological replicates** = 8 libraries, paired-end

### 1. **References**

- `Genome Assembly`: `/dados04/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.fa.gz`
- `Proteins`: `/dados04/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_protein.faa`
- `GTF`: `/dados04/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz`

### 2. **Protein Annotation**

We used `emapper-2.1.3` from `EggNOG v5.0` to obtain annotations (including GO terms) for the proteins
of the genome, based on orthology relationships.

- Results: `/dados04/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/eggnog_anot/eggnog_anot.emapper.annotations`

### 3. **RNAseq processing**

We used a `Nextflow v25.04.7` pipeline, `rnaseq (v3.12.0)` from nf-core
(https://nf-co.re/rnaseq/3.12.0), to preprocess, align and quantify the RNAseq data.

We used the default method of `rnaseq (v3.12.0)`, which uses the `STAR` aligner and `Salmon` to
quantify transcript abundance.

The full preprocessing and alignment report can be found at
`/dados04/jorge/israel_rnaseq/rnaseq/run01/multiqc/star_salmon/multiqc_report.html`

- Generated reads: 449.2M
- Reads aligned to reference: 272.2M (60.6%)

### 4. **Exploratory Analysis**

Transcript-level quantifications produced by Salmon were imported into DESeq2 and aggregated to gene
level. Counts were transformed with `vst` (variance stabilizing transformation) for exploratory
visualisation.

Filtering and significance parameters:

| Parameter | Value |
|---|---|
| `min_count` | 5 |
| `min_samples` | 5 |
| `lfc_threshold` | 1 |
| `padj_threshold` | 0.05 |

A gene was retained if it had ≥ 5 counts in ≥ 5 of the 8 libraries.

- Initial genes: 16,290
- Genes after filtering: 9,330
- Genes removed: 6,960 (42.7%)

- Principal component analysis (on vst-transformed counts):

![PCA_DESeq2](rnaseq/run01/star_salmon/deseq2_qc/all_pca.png)

### 5. **Differential Expression Analysis (DEA)**

We conducted a differential expression analysis (DEA) on the contrast below, filtering for
p-adj (FDR) < 0.05 and |log2 fold change| > 1.

**Contrast:** `FvIAC_HT_vs_FvIAC_VF` — numerator `FvIAC-HT`, denominator `FvIAC-VF`.

Interpretation of the contrast:

- **Upregulated** genes in `FvIAC_HT_vs_FvIAC_VF` are expressed **more** in `FvIAC-HT`
  (virus-harbouring) than in `FvIAC-VF` (virus-free).
- **Downregulated** genes in `FvIAC_HT_vs_FvIAC_VF` are expressed **less** in `FvIAC-HT`
  than in `FvIAC-VF`.

- DESeq2 results: `/dados04/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc`

#### Summary of DEG counts (|LFC| > 1, padj < 0.05)

| Contrast | Description | Upregulated | Downregulated | Total DEGs |
|---|---|---:|---:|---:|
| FvIAC_HT_vs_FvIAC_VF | `FvIAC-HT` vs `FvIAC-VF` | 931 | 1,336 | 2,267 |

### 6. **Functional Enrichment Analysis (ORA)**

To gain insight into the functions and processes represented by the sets of up- and down-regulated
genes, we carried out over-representation analysis (ORA) of Gene Ontology terms.

- **GO:** `topGO` R package (v2.58.0), p-value < 0.05, corrected for multiple testing with the
  Benjamini–Hochberg (FDR) procedure. Background: all 9,330 genes retained after filtering.

- Up-regulated genes:

![up_enrichment](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_upregulated.png)

  Table: `rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_upregulated.csv`

- Down-regulated genes:

![down_enrichment](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_downregulated.png)

  Table: `rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment/GO_BP_downregulated.csv`

Full enrichment results on the server:
`/dados04/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/GO_enrichment`

### 7. **Gene Set Enrichment Analysis (GSEA)**

#### 7.1 Rationale

The ORA in Section 6 operates on the *thresholded* DEG lists (|LFC| > 1, padj < 0.05), which discards
~76% of the tested genes and treats a gene with LFC = 0.98 identically to one with LFC = 0.01. In an
exploratory design with few biological replicates, a biologically meaningful response frequently
appears as a **coordinated, modest shift of many genes in a pathway** rather than as a handful of
large-fold-change outliers — a pattern that ORA is structurally unable to detect.

We therefore complemented the ORA with a hypothesis-driven Gene Set Enrichment Analysis, which uses
the **complete ranking of all tested genes** and applies no significance cutoff. Gene sets were not
taken from a generic database: they were **constructed a priori** to represent the specific phenotypic
hypotheses raised by the accompanying experiments (H1–H5 below), so that each hypothesis maps onto an
explicit, testable transcriptional signature.

#### 7.2 Hypotheses

Phenotypic hypotheses arising from the accompanying experiments:

| # | Hypothesis | Tested by GSEA |
|---|---|---|
| **H1** | The mycovirus alters the **colour** and morphology of fungal colonies. | Colour only |
| **H2** | The mycovirus interferes with the fungus's capacity for **spore production**. | Yes |
| **H3** | The mycovirus modifies the fungus's ability to use different **carbon sources**. | Yes |
| **H4** | The mycovirus modulates the **metabolic pathways** of the fungal host. | Yes |
| **H5** | The mycovirus influences the profile of **volatile organic compounds** emitted by the fungus. | Yes |
| **H6** | The mycovirus affects the fungus's ability to cause **disease in the plant host**. | No |
| **H7** | The mycovirus alters the **tritrophic interaction** among fungus, plant and insect via chemical signals. | No |

**H6 and H7 were not tested.** Both require readouts absent from an axenic, fungus-only
transcriptome — an infected plant host and an insect, respectively. Expression of virulence-associated
genes in culture is not evidence of altered *in planta* disease, and a three-way interaction cannot be
evidenced from a single-organism dataset. Both are deferred to dedicated experimental designs
(*in planta* infection RNAseq; multi-organism assays).

**The morphology component of H1 was likewise excluded**, as colony morphology has no compact,
well-defined transcriptional gene set; only the pigmentation component was retained.

#### 7.3 Gene set construction

Gene sets were built from the eggNOG-mapper v2.1.3 annotation (Section 2) using the `Preferred_name`,
`Description`, `PFAMs`, `EC`, `KEGG_ko`, `KEGG_Pathway`, `GOs` and `CAZy` fields.

**Identifier resolution.** eggNOG annotates protein accessions (`XP_*`), whereas DESeq2 reports gene
identifiers (`FVEG_*`). Proteins were mapped to genes via the `protein_id` / `gene_id` attributes of
the CDS features in the reference GTF. Genes with multiple protein isoforms were collapsed to a single
record by taking the union of their annotation terms.

**Statistical universe.** The background for all tests was restricted to genes that were both
(i) retained after low-count filtering **and** assigned a non-NA adjusted p-value by DESeq2 — i.e.
genes that were actually *eligible* to be called differentially expressed — and (ii) present in the
eggNOG annotation. Genes removed by DESeq2 independent filtering were excluded; including them would
inflate the apparent enrichment of every set.

#### 7.4 Gene set sizes and annotation coverage

**Curated sets (n = 27).**


For each set: the hypothesis it serves, its size in the 8,742-gene universe, and the
exact matching rules used against the eggNOG columns. A gene is included if it matches
**any** criterion (logical OR). Field tags: `Name` = `Preferred_name` (exact, case-insensitive);
`Desc~` = `Description` regex; `PFAM` = `PFAMs`; `EC` = `EC` prefix; `CAZy` = `CAZy` family/class;
`KEGG` = `KEGG_Pathway`; `GO` = `GOs`.

| Gene set | Hyp. | n | Status | Evidence field(s) | Exact match terms / patterns |
|---|---|--:|---|---|---|
| `H1_pigment_bikaverin` | H1 | 0 | not testable | Name; Desc~ | Name: BIK1, BIK2, BIK3, BIK4, BIK5, BIK6, bik1, PKS4 · Desc~: `bikaverin` |
| `H1_pigment_fusarubin` | H1 | 0 | not testable | Name; Desc~ | Name: PGL1, fsr1–fsr6, PKS3 · Desc~: `fusarubin\|naphthoquinone` |
| `H1_pigment_carotenoid` | H1 | 4 | tested | Name; Desc~; EC | Name: carRA, carB, carO, carX, carT, CRTYB · Desc~: `phytoene\|carotene\|carotenoid\|lycopene` · EC: 1.3.5.5, 5.5.1.19 |
| `H1_pigment_melanin_DHN` | H1 | 14 | tested | Name; Desc~; PFAM | Name: PKS12, ARP1, ARP2, ABR1, ABR2, BRN1, BRN2, SCD1, THR1, AYG1 · Desc~: `scytalone\|trihydroxynaphthalene\|tetrahydroxynaphthalene\|melanin\|laccase` · PFAM: Cu-oxidase, Cu-oxidase_2, Cu-oxidase_3 |
| `H1_pigment_PKS_backbone` | H1 | 18 | tested | PFAM; Desc~ | PFAM: ketoacyl-synt, Ketoacyl-synt_C, PKS_AT, PKS_DE, PKS_ER, Acyl_transf_1 · Desc~: `polyketide synthase` |
| `H2_conidiation_core` | H2 | 5 | tested | Name; Desc~ | Name: brlA, abaA, wetA, vosA, medA, stuA, flbA–flbE, fluG, rodA, dewA, con6, con10, chsD · Desc~: `conidiation\|conidiophore\|conidium` |
| `H2_velvet_LaeA` | H2 | 6 | tested | Name; PFAM; Desc~ | Name: veA, velB, velC, vosA, laeA, LAE1, VEL1, VEL2, VEL3 · PFAM: Velvet · Desc~: `velvet` |
| `H2_sporulation_GO` | H2 | 117 | tested | GO | GO: GO:0030435, GO:0043934, GO:0043938, GO:0075307, GO:0048315, GO:1903666, GO:0034293 |
| `H2_MAPK_signaling` | H2 | 8 | tested | Name; Desc~ | Name: FUS3, KSS1, STE7, STE11, STE12, HOG1, PBS2, SSK2, SLT2, BCK1, MKK1, MST11, MST7, GPMK1, FMK1 · Desc~: `mitogen-activated protein kinase` |
| `H2_cAMP_PKA` | H2 | 11 | tested | Name; Desc~ | Name: CYR1, PKA1, PKA2, TPK1, TPK2, BCY1, GPA1, GPA2, GPA3, PDE1, PDE2, RAS1, RAS2 · Desc~: `adenylate cyclase\|cAMP-dependent protein kinase` |
| `H2_light_response` | H2 | 5 | tested | Name; Desc~ | Name: wc-1, wc-2, WCO1, WC1, WC2, vvd, phy1, phy2, cryD, envoy · Desc~: `blue light\|photoreceptor\|white collar\|opsin\|cryptochrome` |
| `H3_CAZy_GH_all` | H3 | 69 | tested | CAZy | CAZy class: GH (any GH family) |
| `H3_CAZy_GT_all` | H3 | 59 | tested | CAZy | CAZy class: GT (any GT family) |
| `H3_CAZy_PL_all` | H3 | 1 |  not testable | CAZy | CAZy class: PL (any PL family) |
| `H3_CAZy_CE_all` | H3 | 2 | not testable | CAZy | CAZy class: CE (any CE family) |
| `H3_CAZy_AA_all` | H3 | 4 | tested | CAZy | CAZy class: AA (any AA family) |
| `H3_sugar_transporters` | H3 | 272 | tested | PFAM; Desc~ | PFAM: Sugar_tr, MFS_1, MFS_2 · Desc~: `sugar transporter\|hexose transporter\|monosaccharide transport` |
| `H3_carbon_catabolite_repr` | H3 | 8 | tested | Name; Desc~ | Name: creA, cre1, creB, creC, creD, snf1, SNF1, hxk1, HXK1, glk1, mig1 · Desc~: `carbon catabolite` |
| `H3_glycolysis_TCA` | H3 | 204 | tested | KEGG | KEGG: ko00010, ko00020, ko00030, ko00051, ko00052, ko00500, ko00620 |
| `H3_alt_carbon_pathways` | H3 | 193 | tested | KEGG | KEGG: ko00040, ko00053, ko00520, ko00561, ko00630, ko00640, ko00650 |
| `H4_secondary_metab_backbone` | H4 | 69 | tested | PFAM; Desc~ | PFAM: ketoacyl-synt, Condensation, AMP-binding, Terpene_synth, Terpene_synth_C, DMATS · Desc~: `polyketide synthase\|nonribosomal peptide synthetase\|terpene synthase\|prenyltransferase` |
| `H4_cytochrome_P450` | H4 | 69 | tested | PFAM; Desc~ | PFAM: p450 · Desc~: `cytochrome P450` |
| `H4_oxidative_stress` | H4 | 27 | tested | Name; Desc~ | Name: CAT1, CAT2, CAT3, SOD1, SOD2, CTT1, TRX1, TRX2, GPX1, GPX2, AP1, YAP1, SKN7 · Desc~: `catalase\|superoxide dismutase\|glutathione peroxidase\|thioredoxin\|peroxiredoxin` |
| `H4_glutathione_detox` | H4 | 28 | tested | KEGG | KEGG: ko00480 (glutathione metabolism) |
| `H4_pentose_phosphate` | H4 | 24 | tested | KEGG | KEGG: ko00030 (pentose phosphate pathway) |
| `H4_sulfur_methionine` | H4 | 67 | tested | KEGG | KEGG: ko00920 (sulfur), ko00270 (cysteine/methionine) |
| `H4_RNAi_antiviral` | H4 | 5 | tested | Name; PFAM; Desc~ | Name: dcl1, dcl2, DCL1, DCL2, ago1, ago2, AGO1, AGO2, qde1, qde2, qde3, rdrp1, sad1 · PFAM: Dicer_dimer, PAZ, Piwi, RNase_III · Desc~: `argonaute\|dicer\|RNA-dependent RNA polymerase` |
| `H5_terpene_synthases` | H5 | 21 | tested | PFAM; Desc~ | PFAM: Terpene_synth, Terpene_synth_C, TRI5, Prenyltrans, polyprenyl_synt · Desc~: `terpene\|sesquiterpene\|trichodiene\|geranyl\|farnesyl` |
| `H5_alcohol_aldehyde_DH` | H5 | 144 | tested | PFAM; EC; Desc~ | PFAM: ADH_N, ADH_zinc_N, Aldedh, ADH_N_2 · EC: 1.1.1.1, 1.2.1.3, 1.1.1.2 · Desc~: `alcohol dehydrogenase\|aldehyde dehydrogenase` |
| `H5_Ehrlich_fusel` | H5 | 4 | tested | Name; Desc~ | Name: ARO8, ARO9, ARO10, PDC1, PDC5, PDC6, ADH1, ADH2, SFA1, BAT1, BAT2 · Desc~: `pyruvate decarboxylase\|aromatic aminotransferase\|branched-chain aminotransferase` |
| `H5_fatty_acid_oxylipin` | H5 | 5 | tested | Name; PFAM; Desc~ | Name: ppo1, ppo2, ppo3, lox1, LOX, ppoA, ppoB, ppoC · PFAM: Lipoxygenase, An_peroxidase · Desc~: `lipoxygenase\|linoleate diol synthase\|oxylipin\|fatty acid desaturase` |
| `H5_esterases_lipases` | H5 | 163 | tested | PFAM; Desc~ | PFAM: Abhydrolase_1, Abhydrolase_3, Lipase_3, COesterase · Desc~: `carboxylesterase\|lipase` |
| `H5_methyltransferases_SAM` | H5 | 70 | tested | PFAM; Desc~ | PFAM: Methyltransf_2, Methyltransf_11, Methyltransf_12, Methyltransf_31 · Desc~: `O-methyltransferase` |


Set sizes were computed **within the statistical universe** (8,742 genes) before any testing. Sets
with fewer than `minSize = 3` genes could not be tested and are reported as **not testable** rather
than as negative results.

**Summary**

| Hypothesis | Sets defined | Tested | Not testable |
|---|---:|---:|---:|
| H1 — pigmentation | 5 | 3 | 2 |
| H2 — sporulation | 6 | 6 | 0 |
| H3 — carbon sources | 9 | 7 | 2 |
| H4 — metabolism | 7 | 7 | 0 |
| H5 — volatiles | 6 | 6 | 0 |
| **Total** | **33** | **29** | **4** |

Full table: `geneset_sizes.csv`

![gsea](rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/hypothesis_genesets/fgsea_hypothesis_sets.png)

#### 7.5 Statistical testing

**GSEA.** Enrichment was computed with the `fgsea` R package (`fgseaMultilevel`). All genes in the
universe were ranked by the **DESeq2 Wald statistic** (`stat` = log2FoldChange / lfcSE). The Wald
statistic was preferred over raw log2FC because, with few replicates, log2FC alone is dominated by
low-expression genes with large standard errors.

**Statistical universe.** The background for all tests was restricted to genes that were both
(i) assigned complete DESeq2 statistics (9,174 genes; Section 5) and (ii) present in the eggNOG
annotation, giving a final universe of **8,742 genes**. Within this universe, 900 of the 931
up-regulated and 1,269 of the 1,336 down-regulated DEGs are annotated and were used for the
complementary ORA.

| Stage | Genes |
|---|---:|
| Genes in the annotation (GTF) | 16,290 |
| Retained after low-count filtering | 9,330 |
| With complete DESeq2 statistics | 9,174 |
| **Statistical universe (tested & eggNOG-annotated)** | **8,742** |
| — of which up-regulated | 900 |
| — of which down-regulated | 1,269 |

**Curated sets.** 33 sets were defined a priori across the five hypotheses:



| Output | Content |
|---|---|
| `fgsea_all_sets.csv` | NES, p, padj, set size, leading edge — all sets |
| `fisher_all_sets.csv` | Fisher's exact test, all sets × direction (up / down / any) |
| `fisher_significant.csv` | Fisher hits at padj < 0.05 |
| `DEG_hits_per_hypothesis_set.csv` | Gene-level table: LFC, lfcSE, stat, padj, eggNOG annotation |
| `geneset_sizes.csv` | Annotation-coverage diagnostic |
| `hypothesis_summary.csv` | H1–H5 overview (best set, NES, padj per hypothesis) |
| `fgsea_hypothesis_sets.png` | NES per set, faceted by hypothesis; **point size = number of genes in the set** |

Output directory:
`/dados04/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc/FvIAC_HT_vs_FvIAC_VF/hypothesis_genesets`

Software: R v4.5.3, `fgsea` v1.36.2, `dplyr_1.2.0`, `tidyr_1.3.2`, `readr_2.2.0`, `ggplot2_4.0.2`


