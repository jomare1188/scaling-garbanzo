# ============================================================================
# HYPOTHESIS-DRIVEN GENE-SET ANALYSIS
# Fusarium verticillioides +/- mycovirus  (FvIAC_HT vs FvIAC_VF)
#
# Goal: move beyond generic GO-BP ORA and test the specific hypotheses
#       coming from the phenotypic experiments, using the richer columns of
#       the eggNOG-mapper output (Preferred_name, Description, PFAMs, EC,
#       KEGG_ko, KEGG_Pathway, CAZy, COG_category).
#
# Tests performed per gene set:
#   (a) Fisher's exact test  -> set enriched among UP / DOWN DEGs ?
#   (b) fgsea (optional)     -> set coordinately shifted along ranked LFC ?
#   (c) per-gene hit table   -> which genes, which LFC, which padj
#
# ID resolution chain (same logic as top_GO_dev2.r):
#   eggNOG query  XP_018741972.1 --(GTF CDS attrs)--> FVEG_*
#   DEG gene_id   XM_* or FVEG_* --(tx2gene)-------> FVEG_*
#
# Required: dplyr, readr, tidyr, ggplot2
# Optional: fgsea  (BiocManager::install("fgsea"))
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})

# ============================================================================
# PARAMETERS -- EDIT THIS BLOCK
# ============================================================================

annotation_file <- "/dados04/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/eggnog_anot/eggnog_anot.emapper.annotations"
gtf_file        <- "/dados04/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz"
tx2gene_file    <- "/dados04/jorge/israel_rnaseq/rnaseq/run01/star_salmon/salmon.merged.tx2gene.tsv"

results_dir <- "/dados04/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc"
contrast    <- "FvIAC_HT_vs_FvIAC_VF"

# FULL DESeq2 results table (ALL tested genes -> defines the statistical
# universe and enables fgsea).
# Columns: gene_id, baseMean_*, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
full_res_file <- file.path(results_dir, contrast,
                           paste0("all_genes_", contrast, ".csv"))

# Ranking metric for fgsea: "stat" (Wald statistic, recommended -- it accounts
# for the standard error, which matters a lot with n=3) or "log2FoldChange".
rank_metric <- "stat"

# DEG files (already produced by your pipeline)
deg_up_file   <- file.path(results_dir, contrast, paste0("DEGs_upregulated_",   contrast, ".csv"))
deg_down_file <- file.path(results_dir, contrast, paste0("DEGs_downregulated_", contrast, ".csv"))

out_dir <- file.path(results_dir, contrast, "hypothesis_genesets")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_set_size <- 3      # gene sets smaller than this are reported but not tested
fdr_cut      <- 0.05
run_fgsea    <- TRUE   # set FALSE to skip

# ============================================================================
# STEP 1 -- PARSE eggNOG ANNOTATIONS
# ============================================================================

cat("[1] Parsing eggNOG annotations...\n")

raw_lines <- readLines(annotation_file)
hdr_idx   <- grep("^#query", raw_lines)[1]
if (is.na(hdr_idx)) stop("Could not find the '#query' header line in the annotation file.")

body_lines <- raw_lines[hdr_idx:length(raw_lines)]
body_lines <- body_lines[!grepl("^##", body_lines)]   # drop '## emapper' comment lines
body_lines <- body_lines[nzchar(body_lines)]

egg <- read.delim(
  text = paste(body_lines, collapse = "\n"),
  sep = "\t", header = TRUE, quote = "", comment.char = "",
  stringsAsFactors = FALSE, check.names = FALSE
)
colnames(egg)[1] <- "query"          # '#query' -> 'query'

# turn eggNOG's "-" placeholder into NA
egg[egg == "-"] <- NA

cat(sprintf("    annotated proteins: %d\n", nrow(egg)))

# ============================================================================
# STEP 2 -- XP_* -> FVEG_* MAP FROM GTF
# ============================================================================

cat("[2] Building XP_ -> FVEG_ map from GTF...\n")

gtf <- read.table(gtf_file, sep = "\t", header = FALSE, comment.char = "#",
                  quote = "", stringsAsFactors = FALSE, fill = TRUE)
gtf <- gtf[gtf[[3]] == "CDS", , drop = FALSE]
attrs <- gtf[[9]]

get_attr <- function(a, key) {
  pat <- paste0(key, ' "([^"]+)"')
  ifelse(grepl(pat, a, perl = TRUE), sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", a, perl = TRUE), NA)
}

prot2gene <- data.frame(
  protein_id = get_attr(attrs, "protein_id"),
  gene_id    = get_attr(attrs, "gene_id"),
  stringsAsFactors = FALSE
) %>% filter(!is.na(protein_id), !is.na(gene_id)) %>% distinct()

xp2fveg <- setNames(prot2gene$gene_id, prot2gene$protein_id)

egg$gene <- unname(xp2fveg[egg$query])
n_unmapped <- sum(is.na(egg$gene))
cat(sprintf("    mapped: %d   unmapped: %d\n", sum(!is.na(egg$gene)), n_unmapped))
egg <- egg %>% filter(!is.na(gene))

# collapse to one row per gene (a gene can have several protein isoforms):
# keep the best-scoring annotation, but concatenate the multi-value fields
collapse_uniq <- function(x) {
  v <- unique(unlist(strsplit(na.omit(x), ",")))
  if (length(v) == 0) NA_character_ else paste(v, collapse = ",")
}

ann <- egg %>%
  group_by(gene) %>%
  summarise(
    Preferred_name = collapse_uniq(Preferred_name),
    Description    = paste(unique(na.omit(Description)), collapse = " | "),
    PFAMs          = collapse_uniq(PFAMs),
    EC             = collapse_uniq(EC),
    KEGG_ko        = collapse_uniq(KEGG_ko),
    KEGG_Pathway   = collapse_uniq(KEGG_Pathway),
    CAZy           = collapse_uniq(CAZy),
    COG_category   = collapse_uniq(COG_category),
    GOs            = collapse_uniq(GOs),
    .groups = "drop"
  ) %>%
  mutate(Description = ifelse(Description == "", NA, Description))

cat(sprintf("    unique annotated genes: %d\n", nrow(ann)))

# ============================================================================
# STEP 3 -- LOAD DEGs AND FULL RESULTS, RESOLVE IDs
# ============================================================================

cat("[3] Loading DESeq2 output...\n")

# DEG / results files already use FVEG_* gene_ids, but keep the tx2gene
# translation as a safety net in case any XM_* transcript ids slip through.
if (!is.null(tx2gene_file) && file.exists(tx2gene_file)) {
  tx2gene <- read.table(tx2gene_file, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "")
  tx_map <- setNames(c(tx2gene$gene_id, tx2gene$gene_id),
                     c(tx2gene$transcript_id, tx2gene$gene_id))
} else {
  tx_map <- character(0)
}

to_fveg <- function(ids) {
  ids <- as.character(ids)
  if (length(tx_map) == 0) return(ids)
  out <- unname(tx_map[ids])
  ifelse(is.na(out), ids, out)
}

up   <- read_csv(deg_up_file,   show_col_types = FALSE) %>% mutate(gene = to_fveg(gene_id))
down <- read_csv(deg_down_file, show_col_types = FALSE) %>% mutate(gene = to_fveg(gene_id))

up_genes   <- unique(up$gene)
down_genes <- unique(down$gene)
deg_genes  <- union(up_genes, down_genes)
cat(sprintf("    DEG files: %d up, %d down\n", length(up_genes), length(down_genes)))

if (!is.null(full_res_file) && file.exists(full_res_file)) {
  full_res <- read_csv(full_res_file, show_col_types = FALSE) %>%
    mutate(gene = to_fveg(gene_id))

  # Universe = genes that DESeq2 actually tested AND could have called DE.
  # padj = NA means the gene was removed by independent filtering, so it was
  # never eligible to be a DEG -> excluding it keeps Fisher's test honest.
  tested <- full_res %>% filter(!is.na(padj)) %>% pull(gene) %>% unique()
  universe <- intersect(tested, ann$gene)

  cat(sprintf("    all_genes rows: %d | with padj (testable): %d | testable & annotated: %d\n",
              nrow(full_res), length(tested), length(universe)))
} else {
  full_res <- NULL
  universe <- ann$gene
  cat("    WARNING: no full results file -> universe = all annotated genes.\n")
  cat("             fgsea will be skipped.\n")
  run_fgsea <- FALSE
}

up_genes   <- intersect(up_genes,   universe)
down_genes <- intersect(down_genes, universe)
deg_genes  <- intersect(deg_genes,  universe)
cat(sprintf("    UP in universe: %d | DOWN in universe: %d\n",
            length(up_genes), length(down_genes)))

ann_u <- ann %>% filter(gene %in% universe)

# ============================================================================
# STEP 4 -- GENE-SET BUILDERS
# ============================================================================
# Helpers that query the eggNOG columns. All are case-insensitive and
# operate only on the universe.

by_name <- function(names_vec) {
  pat <- paste0("(^|,)(", paste(names_vec, collapse = "|"), ")($|,)")
  ann_u$gene[grepl(pat, ann_u$Preferred_name, ignore.case = TRUE)]
}
by_desc <- function(regex) {
  ann_u$gene[grepl(regex, ann_u$Description, ignore.case = TRUE, perl = TRUE)]
}
by_pfam <- function(pfam_vec) {
  pat <- paste0("(^|,)(", paste(pfam_vec, collapse = "|"), ")($|,)")
  ann_u$gene[grepl(pat, ann_u$PFAMs, ignore.case = TRUE)]
}
by_ec <- function(ec_prefix_vec) {
  pat <- paste0("(^|,)(", paste(gsub("\\.", "\\\\.", ec_prefix_vec), collapse = "|"), ")")
  ann_u$gene[grepl(pat, ann_u$EC)]
}
by_cazy_class <- function(cls) {   # "GH", "GT", "PL", "CE", "AA", "CBM"
  pat <- paste0("(^|,)", cls, "[0-9]")
  ann_u$gene[grepl(pat, ann_u$CAZy)]
}
by_cazy_family <- function(fam_vec) {  # e.g. c("GH28","PL1")
  pat <- paste0("(^|,)(", paste(fam_vec, collapse = "|"), ")($|_|,)")
  ann_u$gene[grepl(pat, ann_u$CAZy)]
}

# ---------------------------------------------------------------------------
# H1 -- PIGMENTATION / COLONY MORPHOLOGY
# ---------------------------------------------------------------------------
sets <- list()

sets[["H1_pigment_bikaverin"]]      <- unique(c(by_name(c("BIK1","BIK2","BIK3","BIK4","BIK5","BIK6","bik1","PKS4")),
                                                by_desc("bikaverin")))
sets[["H1_pigment_fusarubin"]]      <- unique(c(by_name(c("PGL1","fsr1","fsr2","fsr3","fsr4","fsr5","fsr6","PKS3")),
                                                by_desc("fusarubin|naphthoquinone")))
sets[["H1_pigment_carotenoid"]]     <- unique(c(by_name(c("carRA","carB","carO","carX","carT","CRTYB")),
                                                by_desc("phytoene|carotene|carotenoid|lycopene"),
                                                by_ec(c("1.3.5.5","5.5.1.19"))))
sets[["H1_pigment_melanin_DHN"]]    <- unique(c(by_name(c("PKS12","ARP1","ARP2","ABR1","ABR2","BRN1","BRN2","SCD1","THR1","AYG1")),
                                                by_desc("scytalone|trihydroxynaphthalene|tetrahydroxynaphthalene|melanin|laccase"),
                                                by_pfam(c("Cu-oxidase","Cu-oxidase_2","Cu-oxidase_3"))))
sets[["H1_pigment_PKS_backbone"]]   <- unique(c(by_pfam(c("ketoacyl-synt","Ketoacyl-synt_C","PKS_AT","PKS_DE","PKS_ER","Acyl_transf_1")),
                                                by_desc("polyketide synthase")))

# ---------------------------------------------------------------------------
# H2 -- SPORULATION / CONIDIATION
# ---------------------------------------------------------------------------
sets[["H2_conidiation_core"]]       <- unique(c(by_name(c("brlA","abaA","wetA","vosA","medA","stuA","flbA","flbB","flbC","flbD","flbE","fluG","rodA","dewA","con6","con10","chsD")),
                                                by_desc("conidiation|conidiophore|conidium")))
sets[["H2_velvet_LaeA"]]            <- unique(c(by_name(c("veA","velB","velC","vosA","laeA","LAE1","VEL1","VEL2","VEL3")),
                                                by_pfam(c("Velvet")),
                                                by_desc("velvet")))
sets[["H2_sporulation_GO"]]         <- ann_u$gene[grepl("GO:0030435|GO:0043934|GO:0043938|GO:0075307|GO:0048315|GO:1903666|GO:0034293", ann_u$GOs)]
sets[["H2_MAPK_signaling"]]         <- unique(c(by_name(c("FUS3","KSS1","STE7","STE11","STE12","HOG1","PBS2","SSK2","SLT2","BCK1","MKK1","MST11","MST7","GPMK1","FMK1")),
                                                by_desc("mitogen-activated protein kinase")))
sets[["H2_cAMP_PKA"]]               <- unique(c(by_name(c("CYR1","PKA1","PKA2","TPK1","TPK2","BCY1","GPA1","GPA2","GPA3","PDE1","PDE2","RAS1","RAS2")),
                                                by_desc("adenylate cyclase|cAMP-dependent protein kinase")))
sets[["H2_light_response"]]         <- unique(c(by_name(c("wc-1","wc-2","WCO1","WC1","WC2","vvd","phy1","phy2","cryD","envoy")),
                                                by_desc("blue light|photoreceptor|white collar|opsin|cryptochrome")))

# ---------------------------------------------------------------------------
# H3 -- CARBON SOURCE UTILIZATION
# ---------------------------------------------------------------------------
sets[["H3_CAZy_GH_all"]]            <- by_cazy_class("GH")
sets[["H3_CAZy_GT_all"]]            <- by_cazy_class("GT")
sets[["H3_CAZy_PL_all"]]            <- by_cazy_class("PL")
sets[["H3_CAZy_CE_all"]]            <- by_cazy_class("CE")
sets[["H3_CAZy_AA_all"]]            <- by_cazy_class("AA")
sets[["H3_sugar_transporters"]]     <- unique(c(by_pfam(c("Sugar_tr","MFS_1","MFS_2")),
                                                by_desc("sugar transporter|hexose transporter|monosaccharide transport")))
sets[["H3_carbon_catabolite_repr"]] <- unique(c(by_name(c("creA","cre1","creB","creC","creD","snf1","SNF1","hxk1","HXK1","glk1","mig1")),
                                                by_desc("carbon catabolite")))
sets[["H3_glycolysis_TCA"]]         <- ann_u$gene[grepl("ko00010|ko00020|ko00030|ko00051|ko00052|ko00500|ko00620", ann_u$KEGG_Pathway)]
sets[["H3_alt_carbon_pathways"]]    <- ann_u$gene[grepl("ko00040|ko00053|ko00520|ko00561|ko00630|ko00640|ko00650", ann_u$KEGG_Pathway)]

# ---------------------------------------------------------------------------
# H4 -- HOST METABOLIC PATHWAYS  (handled mostly by the KEGG loop below)
# ---------------------------------------------------------------------------
sets[["H4_secondary_metab_backbone"]] <- unique(c(by_pfam(c("ketoacyl-synt","Condensation","AMP-binding","Terpene_synth","Terpene_synth_C","DMATS")),
                                                  by_desc("polyketide synthase|nonribosomal peptide synthetase|terpene synthase|prenyltransferase")))
sets[["H4_cytochrome_P450"]]          <- unique(c(by_pfam(c("p450")), by_desc("cytochrome P450")))
sets[["H4_oxidative_stress"]]         <- unique(c(by_name(c("CAT1","CAT2","CAT3","SOD1","SOD2","CTT1","TRX1","TRX2","GPX1","GPX2","AP1","YAP1","SKN7")),
                                                  by_desc("catalase|superoxide dismutase|glutathione peroxidase|thioredoxin|peroxiredoxin")))
sets[["H4_glutathione_detox"]]        <- ann_u$gene[grepl("ko00480", ann_u$KEGG_Pathway)]
sets[["H4_pentose_phosphate"]]        <- ann_u$gene[grepl("ko00030", ann_u$KEGG_Pathway)]
sets[["H4_sulfur_methionine"]]        <- ann_u$gene[grepl("ko00920|ko00270", ann_u$KEGG_Pathway)]
sets[["H4_RNAi_antiviral"]]           <- unique(c(by_name(c("dcl1","dcl2","DCL1","DCL2","ago1","ago2","AGO1","AGO2","qde1","qde2","qde3","rdrp1","sad1")),
                                                  by_pfam(c("Dicer_dimer","PAZ","Piwi","RNase_III")),
                                                  by_desc("argonaute|dicer|RNA-dependent RNA polymerase")))

# ---------------------------------------------------------------------------
# H5 -- VOLATILE ORGANIC COMPOUNDS (proxy sets)
# ---------------------------------------------------------------------------
sets[["H5_terpene_synthases"]]      <- unique(c(by_pfam(c("Terpene_synth","Terpene_synth_C","TRI5","Prenyltrans","polyprenyl_synt")),
                                                by_desc("terpene|sesquiterpene|trichodiene|geranyl|farnesyl")))
sets[["H5_alcohol_aldehyde_DH"]]    <- unique(c(by_pfam(c("ADH_N","ADH_zinc_N","Aldedh","ADH_N_2")),
                                                by_ec(c("1.1.1.1","1.2.1.3","1.1.1.2")),
                                                by_desc("alcohol dehydrogenase|aldehyde dehydrogenase")))
sets[["H5_Ehrlich_fusel"]]          <- unique(c(by_name(c("ARO8","ARO9","ARO10","PDC1","PDC5","PDC6","ADH1","ADH2","SFA1","BAT1","BAT2")),
                                                by_desc("pyruvate decarboxylase|aromatic aminotransferase|branched-chain aminotransferase")))
sets[["H5_fatty_acid_oxylipin"]]    <- unique(c(by_name(c("ppo1","ppo2","ppo3","lox1","LOX","ppoA","ppoB","ppoC")),
                                                by_pfam(c("Lipoxygenase","An_peroxidase")),
                                                by_desc("lipoxygenase|linoleate diol synthase|oxylipin|fatty acid desaturase")))
sets[["H5_esterases_lipases"]]      <- unique(c(by_pfam(c("Abhydrolase_1","Abhydrolase_3","Lipase_3","COesterase")),
                                                by_desc("carboxylesterase|lipase")))
sets[["H5_methyltransferases_SAM"]] <- unique(c(by_pfam(c("Methyltransf_2","Methyltransf_11","Methyltransf_12","Methyltransf_31")),
                                                by_desc("O-methyltransferase")))

# H6 (virulence) and H7 (tritrophic signalling) are out of scope for this
# analysis -- see the feasibility ranking. Add them back here if that changes.

# ============================================================================
# STEP 5 -- KEGG PATHWAY SETS (data-driven, for H4/H3)
# ============================================================================

kegg_long <- ann_u %>%
  filter(!is.na(KEGG_Pathway)) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  filter(grepl("^ko[0-9]{5}$", KEGG_Pathway)) %>%
  distinct(gene, KEGG_Pathway)

kegg_sets <- split(kegg_long$gene, kegg_long$KEGG_Pathway)
kegg_sets <- kegg_sets[lengths(kegg_sets) >= min_set_size]
cat(sprintf("[4] KEGG pathway sets (>= %d genes): %d\n", min_set_size, length(kegg_sets)))

# CAZy family sets (data-driven, for H3/H6)
cazy_long <- ann_u %>%
  filter(!is.na(CAZy)) %>%
  separate_rows(CAZy, sep = ",") %>%
  mutate(CAZy = sub("_.*$", "", CAZy)) %>%
  filter(grepl("^(GH|GT|PL|CE|AA|CBM)[0-9]+$", CAZy)) %>%
  distinct(gene, CAZy)

cazy_sets <- split(cazy_long$gene, cazy_long$CAZy)
cazy_sets <- cazy_sets[lengths(cazy_sets) >= min_set_size]
cat(sprintf("    CAZy family sets (>= %d genes): %d\n", min_set_size, length(cazy_sets)))

# ============================================================================
# STEP 6 -- FISHER'S EXACT TEST
# ============================================================================

fisher_set <- function(set_genes, deg, universe, set_name, direction) {
  set_genes <- intersect(unique(set_genes), universe)
  k  <- length(intersect(set_genes, deg))     # in set & DEG
  m  <- length(set_genes)                     # in set
  n  <- length(universe) - m                  # not in set
  dd <- length(deg)                           # DEGs
  if (m == 0) return(NULL)
  mat <- matrix(c(k, m - k, dd - k, n - (dd - k)), nrow = 2)
  ft  <- fisher.test(mat, alternative = "greater")
  data.frame(
    set        = set_name,
    direction  = direction,
    set_size   = m,
    n_DEG      = dd,
    observed   = k,
    expected   = round(m * dd / length(universe), 2),
    fold_enr   = round(ifelse(m * dd == 0, NA, k / (m * dd / length(universe))), 2),
    odds_ratio = round(unname(ft$estimate), 2),
    pvalue     = ft$p.value,
    genes      = paste(sort(intersect(set_genes, deg)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

run_all <- function(set_list, label) {
  res <- bind_rows(lapply(names(set_list), function(s) {
    bind_rows(
      fisher_set(set_list[[s]], up_genes,   universe, s, "up"),
      fisher_set(set_list[[s]], down_genes, universe, s, "down"),
      fisher_set(set_list[[s]], deg_genes,  universe, s, "any")
    )
  }))
  if (nrow(res) == 0) return(res)
  res %>%
    group_by(direction) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    ungroup() %>%
    mutate(collection = label) %>%
    arrange(pvalue) %>%
    select(collection, set, direction, set_size, observed, expected,
           fold_enr, odds_ratio, pvalue, padj, genes)
}

cat("[5] Fisher's exact tests...\n")
res_hyp  <- run_all(sets[lengths(sets) >= min_set_size], "curated_hypothesis")
res_kegg <- run_all(kegg_sets, "KEGG_pathway")
res_cazy <- run_all(cazy_sets, "CAZy_family")

res_all <- bind_rows(res_hyp, res_kegg, res_cazy)
write_csv(res_all, file.path(out_dir, "fisher_all_sets.csv"))
write_csv(res_all %>% filter(padj < fdr_cut),
          file.path(out_dir, "fisher_significant.csv"))

cat("\n--- Curated hypothesis sets (top hits) ---\n")
print(as.data.frame(res_hyp %>% filter(direction != "any") %>%
                      select(set, direction, set_size, observed, expected, fold_enr, pvalue, padj) %>%
                      head(25)), row.names = FALSE)

# empty / tiny sets -> useful diagnostic (annotation gaps, not biology)
empty_sets <- data.frame(
  set = names(sets),
  size_in_universe = sapply(sets, function(s) length(intersect(s, universe)))
) %>% arrange(size_in_universe)
write_csv(empty_sets, file.path(out_dir, "geneset_sizes.csv"))
cat("\nSets with < min_set_size genes in the universe (check annotation coverage):\n")
print(empty_sets %>% filter(size_in_universe < min_set_size), row.names = FALSE)

# ============================================================================
# STEP 7 -- PER-GENE HIT TABLE
# ============================================================================

cat("\n[6] Building per-gene hit tables...\n")

set_membership <- bind_rows(lapply(names(sets), function(s) {
  g <- intersect(sets[[s]], universe)
  if (length(g) == 0) return(NULL)
  data.frame(gene = g, set = s, stringsAsFactors = FALSE)
}))

deg_stats <- bind_rows(
  data.frame(gene = up_genes,   direction = "up",   stringsAsFactors = FALSE),
  data.frame(gene = down_genes, direction = "down", stringsAsFactors = FALSE)
)

# pull LFC / padj / baseMean from the all_genes table (authoritative source)
if (!is.null(full_res)) {
  deg_stats <- deg_stats %>%
    left_join(full_res %>%
                select(gene, baseMean, log2FoldChange, lfcSE, stat, padj) %>%
                distinct(gene, .keep_all = TRUE),
              by = "gene")
}

hits <- set_membership %>%
  inner_join(deg_stats, by = "gene") %>%
  left_join(ann_u %>% select(gene, Preferred_name, Description, PFAMs, EC, CAZy, KEGG_Pathway),
            by = "gene") %>%
  arrange(set, desc(abs(log2FoldChange)))

write_csv(hits, file.path(out_dir, "DEG_hits_per_hypothesis_set.csv"))
cat(sprintf("    %d gene-set/DEG hits written.\n", nrow(hits)))

# ============================================================================
# STEP 8 -- fgsea ON RANKED LFC (optional, uses ALL tested genes)
# ============================================================================

if (run_fgsea && requireNamespace("fgsea", quietly = TRUE)) {
  cat(sprintf("[7] Running fgsea on genes ranked by '%s'...\n", rank_metric))
  library(fgsea)

  if (!rank_metric %in% colnames(full_res))
    stop("rank_metric '", rank_metric, "' not found in the all_genes table.")

  rank_df <- full_res %>%
    filter(gene %in% universe, !is.na(.data[[rank_metric]])) %>%
    group_by(gene) %>%
    summarise(metric = mean(.data[[rank_metric]]), .groups = "drop")

  ranks <- sort(setNames(rank_df$metric, rank_df$gene), decreasing = TRUE)
  cat(sprintf("    ranked genes: %d  (metric: %s)\n", length(ranks), rank_metric))

  all_sets <- c(sets[lengths(sets) >= min_set_size], kegg_sets, cazy_sets)
  all_sets <- lapply(all_sets, function(s) intersect(s, names(ranks)))
  all_sets <- all_sets[lengths(all_sets) >= min_set_size]

  fg <- fgsea(pathways = all_sets, stats = ranks, minSize = min_set_size,
              maxSize = 1000, nPermSimple = 10000)
  fg <- fg %>% as.data.frame() %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    arrange(pval)
  write_csv(fg, file.path(out_dir, "fgsea_all_sets.csv"))
  fgsea_res <- fg   # kept for the hypothesis-level summary below

  cat("\n--- fgsea (padj < 0.05) ---\n")
  print(fg %>% filter(padj < fdr_cut) %>%
          select(pathway, size, NES, pval, padj) %>% head(30), row.names = FALSE)

  # ---- plot: all curated sets, point size = number of genes in the set ----
  fg_cur <- fg %>%
    filter(pathway %in% names(sets)) %>%
    mutate(
      hypothesis = sub("_.*$", "", pathway),
      label      = paste0(pathway, "  (n=", size, ")"),
      sig        = ifelse(padj < fdr_cut, paste0("padj < ", fdr_cut),
                          ifelse(padj < 0.25, "padj < 0.25", "n.s."))
    )

  if (nrow(fg_cur) > 0) {
    p <- ggplot(fg_cur, aes(x = NES, y = reorder(label, NES))) +
      geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
      geom_segment(aes(x = 0, xend = NES, yend = reorder(label, NES)),
                   colour = "grey70", linewidth = 0.4) +
      geom_point(aes(size = size, colour = sig)) +
      scale_size_continuous(name = "genes in set", range = c(2.5, 9),
                            breaks = scales::breaks_pretty(4)) +
      scale_colour_manual(
        name   = "significance",
        values = setNames(c("#c0392b", "#e0a63a", "grey65"),
                          c(paste0("padj < ", fdr_cut), "padj < 0.25", "n.s.")),
        drop   = FALSE
      ) +
      facet_grid(hypothesis ~ ., scales = "free_y", space = "free_y") +
      labs(
        x        = "Normalized Enrichment Score (NES)",
        y        = NULL,
        title    = paste0("Hypothesis gene sets - ", contrast),
        subtitle = "NES > 0: set shifted UP in FvIAC_HT relative to FvIAC_VF\npoint size = number of genes in the set"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        strip.background   = element_rect(fill = "grey92", colour = NA),
        strip.text.y       = element_text(face = "bold"),
        plot.subtitle      = element_text(size = 10, colour = "grey30")
      )

    ggsave(file.path(out_dir, "fgsea_hypothesis_sets.png"), p,
           width = 26, height = 22, units = "cm", dpi = 300)
    cat("    plot written: fgsea_hypothesis_sets.png\n")
  }
} else if (run_fgsea) {
  cat("[7] fgsea not installed - skipping. BiocManager::install('fgsea')\n")
}

# ============================================================================
# STEP 9 -- HYPOTHESIS-LEVEL SUMMARY
# ============================================================================

cat("\n[8] Hypothesis-level summary...\n")

summary_h <- res_all %>%
  filter(collection == "curated_hypothesis", direction != "any") %>%
  mutate(hypothesis = sub("_.*$", "", set)) %>%
  group_by(hypothesis) %>%
  summarise(
    sets_tested      = n_distinct(set),
    sets_sig_fisher  = n_distinct(set[padj < fdr_cut]),
    best_set_fisher  = set[which.min(pvalue)][1],
    best_dir         = direction[which.min(pvalue)][1],
    best_padj_fisher = min(padj, na.rm = TRUE),
    total_DEG_hits   = sum(observed),
    .groups = "drop"
  )

# add the fgsea view (the statistic I'd trust more with n=3)
if (exists("fgsea_res")) {
  fg_h <- fgsea_res %>%
    filter(pathway %in% names(sets)) %>%
    mutate(hypothesis = sub("_.*$", "", pathway)) %>%
    group_by(hypothesis) %>%
    summarise(
      sets_sig_fgsea  = sum(padj < fdr_cut, na.rm = TRUE),
      best_set_fgsea  = pathway[which.min(pval)][1],
      best_NES        = round(NES[which.min(pval)][1], 2),
      best_padj_fgsea = min(padj, na.rm = TRUE),
      .groups = "drop"
    )
  summary_h <- summary_h %>% left_join(fg_h, by = "hypothesis")
}

summary_h <- summary_h %>%
  arrange(pmin(best_padj_fisher,
               if ("best_padj_fgsea" %in% names(.)) best_padj_fgsea else Inf,
               na.rm = TRUE))

write_csv(summary_h, file.path(out_dir, "hypothesis_summary.csv"))
print(as.data.frame(summary_h), row.names = FALSE)

cat("\n=================================================================\n")
cat("DONE. Output in:\n  ", out_dir, "\n")
cat("  fisher_all_sets.csv            all sets x direction\n")
cat("  fisher_significant.csv         padj < ", fdr_cut, "\n")
cat("  DEG_hits_per_hypothesis_set.csv  gene-level table (LFC, padj, annot)\n")
cat("  geneset_sizes.csv              annotation-coverage diagnostic\n")
cat("  fgsea_all_sets.csv             GSEA on ranked LFC\n")
cat("  hypothesis_summary.csv         H1..H5 overview\n")
cat("  fgsea_hypothesis_sets.png      NES per set, point size = set size\n")
cat("=================================================================\n")
