library("ggplot2")
library("dplyr")
library("DESeq2")
library("tximport")
library("tidyverse")
library("wesanderson")
library("viridis")
library("ggsci")
library("showtext")
showtext_auto()

# NOTE: for best-quality text rendering in the SVG outputs, install "svglite"
# (install.packages("svglite")) — ggsave() will use it automatically if present.


# ============================================================================
# PARAMETERS - EDIT THIS SECTION TO CONFIGURE YOUR ANALYSIS
# ============================================================================

colors <- (pal_npg("nrc")(5))

base_dir <- "/dados02/jorge/israel_rnaseq/rnaseq"
output_dir <- "/dados02/jorge/israel_rnaseq/rnaseq/run01/star_salmon/deseq2_qc"
# Design formula
design_formula <- ~ group

# Gene filtering thresholds
min_count <- 5
min_samples <- 5

# DEG thresholds
lfc_threshold <- 1
padj_threshold <- 0.05

# PCA parameters
pca_top_genes <- 100000000

# Define all pairwise contrasts between the 4 groups
# Format: numerator vs denominator (positive LFC = higher in numerator)
groups <- c("FvIAC-VF", "FvIAC-HT")

contrasts_to_run <- do.call(c, lapply(seq_along(groups), function(i) {
  lapply(seq_along(groups), function(j) {
    if (j >= i) return(NULL)
    num <- groups[i]
    den <- groups[j]
    # Build a safe filename: replace hyphens and underscores
    safe_num <- gsub("[-]", "_", num)
    safe_den <- gsub("[-]", "_", den)
    list(
      numerator   = num,
      denominator = den,
      name        = paste0(safe_num, "_vs_", safe_den),
      description = paste0(num, " vs ", den)
    )
  })
}))
contrasts_to_run <- Filter(Negate(is.null), contrasts_to_run)

cat(sprintf("Contrasts to run (%d total):\n", length(contrasts_to_run)))
for (ct in contrasts_to_run) cat(sprintf("  - %s\n", ct$description))

# ----------------------------------------------------------------------------
# Helper: save a ggplot object as PNG + PDF + SVG in one call
# ----------------------------------------------------------------------------
save_plot_multi <- function(plot, basepath, width = 19.5, height = 15, dpi = 320) {
  for (ext in c("png", "pdf", "svg")) {
    fname <- paste0(basepath, ".", ext)
    ggsave(plot, filename = fname, units = "cm", width = width, height = height, dpi = dpi)
    cat(sprintf("  ✓ saved: %s\n", fname))
  }
}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=================================================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS PIPELINE\n")
cat("=================================================================\n\n")

# ----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# ----------------------------------------------------------------------------

cat("[1/6] Loading metadata and count data...\n")

# Read metadata - filter for infected samples
metadata <- read.table(
  file.path(base_dir, "/samplesheet.csv"),
  header = TRUE, 
  sep = ","
)
# here filter by genotype
#tmp <- metadata %>% filter(genotype == "35S::AtCKX2")
#metadata <- tmp

# Display sample distribution
cat("\nSample distribution:\n")
sample_table <- table(metadata$group)
print(sample_table)
cat("\nTotal infected samples:", nrow(metadata), "\n")

# Load Salmon quantification files
sample_files <- file.path(
  base_dir, 
  "/run01/star_salmon",
  metadata$sample, 
  "quant.sf"
)
names(sample_files) <- metadata$sample

# Load transcript-to-gene mapping
tx2gene <- read.table(
  file.path(base_dir, "/run01/star_salmon/salmon.merged.tx2gene.tsv"),
  header = TRUE
)

# Import count data
count_data <- tximport(
  files = sample_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE,
  ignoreAfterBar = TRUE
)



coldata <- metadata

rownames(coldata) <- coldata$sample

# Verify alignment
stopifnot(all(colnames(count_data$counts) == rownames(coldata)))

cat("✓ Data loaded successfully\n")
cat(sprintf("  Samples: %d\n", ncol(count_data$counts)))
cat(sprintf("  Genes: %d\n", nrow(count_data$counts)))

# ----------------------------------------------------------------------------
# 2. CREATE DESEQ2 OBJECT AND FILTER GENES
# ----------------------------------------------------------------------------

cat("\n[2/6] Creating DESeq2 object and filtering genes...\n")

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(
  txi = count_data,
  colData = coldata,
  design = design_formula
)

cat(sprintf("  Initial genes: %d\n", nrow(dds)))

# Filter low-count genes
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds_filtered <- dds[keep, ]

cat(sprintf("  Genes after filtering: %d\n", nrow(dds_filtered)))
cat(sprintf("  Genes removed: %d (%.1f%%)\n", 
            sum(!keep), 100 * sum(!keep) / nrow(dds)))

# ----------------------------------------------------------------------------
# 3. RUN DESEQ2
# ----------------------------------------------------------------------------

cat("\n[3/6] Running DESeq2 analysis...\n")
dds_filtered <- DESeq(dds_filtered, parallel = TRUE)
cat("✓ DESeq2 analysis complete\n")

# ----------------------------------------------------------------------------
# 4. GENERATE QC PLOTS
# ----------------------------------------------------------------------------

cat("\n[4/6] Generating quality control plots...\n")

# VST transformation
vst <- varianceStabilizingTransformation(dds_filtered)

# PCA plot
#colors <- wes_palette("Darjeeling1", 5, type = "discrete")
#colors <- viridis(4)
pca_data <- plotPCA(
  vst, 
  intgroup = c("group"),
  returnData = TRUE, 
  ntop = pca_top_genes
)
percentVar <- round(100 * attr(pca_data, "percentVar"))


p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  scale_colour_manual(values = colors) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 22),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

pca_basepath <- file.path(output_dir, "all_pca")
save_plot_multi(p, pca_basepath, width = 19.5, height = 15)
cat(sprintf("✓ PCA plots saved (png/pdf/svg): %s.*\n", pca_basepath))

# ----------------------------------------------------------------------------
# 5. RUN DIFFERENTIAL EXPRESSION CONTRASTS
# ----------------------------------------------------------------------------

cat("\n[5/6] Running differential expression contrasts...\n")
cat(sprintf("  Number of contrasts: %d\n\n", length(contrasts_to_run)))

# Function to build a minimalistic volcano plot for one contrast.
# Uses the same color palette as the PCA plot (colors[1] = up, colors[2] = down).
make_volcano_plot <- function(res_df, contrast_info, lfc_thresh, padj_thresh, colors) {

  volcano_df <- res_df %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig_category = case_when(
        padj < padj_thresh & log2FoldChange >  lfc_thresh ~ "Up",
        padj < padj_thresh & log2FoldChange < -lfc_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )

  # Cap Inf values (padj == 0) so points aren't dropped from the plot
  finite_max <- max(volcano_df$neg_log10_padj[is.finite(volcano_df$neg_log10_padj)], na.rm = TRUE)
  volcano_df$neg_log10_padj[is.infinite(volcano_df$neg_log10_padj)] <- finite_max * 1.05

  n_up   <- sum(volcano_df$sig_category == "Up")
  n_down <- sum(volcano_df$sig_category == "Down")

  volcano_colors <- c("Up" = colors[1], "Down" = colors[2], "NS" = "grey80")

  x_range <- range(volcano_df$log2FoldChange, na.rm = TRUE)
  y_max   <- max(volcano_df$neg_log10_padj, na.rm = TRUE)

  v <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = sig_category)) +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_color_manual(values = volcano_colors, guide = "none") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    annotate("text", x = x_range[2], y = y_max, hjust = 1, vjust = 1,
             label = paste0("Up: ", n_up), color = colors[1], size = 4.5,
             family = "Times New Roman") +
    annotate("text", x = x_range[1], y = y_max, hjust = 0, vjust = 1,
             label = paste0("Down: ", n_down), color = colors[2], size = 4.5,
             family = "Times New Roman") +
    labs(
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ "adjusted p-value"),
      title = contrast_info$description
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Times New Roman", size = 22),
      plot.title = element_text(size = 16, hjust = 0.5),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )

  return(v)
}

# Function to run a single contrast
run_contrast <- function(contrast_info, dds, output_dir, lfc_thresh, padj_thresh, colors) {

  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("Processing: %s\n", contrast_info$name))
  cat(sprintf("Description: %s\n", contrast_info$description))
  cat(sprintf("Samples in analysis: %d\n", ncol(dds)))

  # Extract results directly — design is ~ group, no refitting needed
  dea_contrast <- results(
    dds,
    contrast = c("group", contrast_info$numerator, contrast_info$denominator),
    alpha    = padj_thresh,
    parallel = TRUE
  )

  # Calculate per-condition base means
  numerator_samples   <- colData(dds)[["group"]] == contrast_info$numerator
  denominator_samples <- colData(dds)[["group"]] == contrast_info$denominator

  baseMean_numerator   <- rowMeans(counts(dds, normalized = TRUE)[, numerator_samples, drop = FALSE])
  baseMean_denominator <- rowMeans(counts(dds, normalized = TRUE)[, denominator_samples, drop = FALSE])

  # Assemble results data frame
  res_df <- data.frame(
    gene_id              = rownames(dea_contrast),
    baseMean_numerator   = baseMean_numerator,
    baseMean_denominator = baseMean_denominator,
    as.data.frame(dea_contrast)
  ) %>%
    filter(complete.cases(.))

  colnames(res_df)[2:3] <- c(
    paste0("baseMean_", gsub("[-]", "_", contrast_info$numerator)),
    paste0("baseMean_", gsub("[-]", "_", contrast_info$denominator))
  )

  # Filter DEGs
  deg_all  <- res_df %>% filter(abs(log2FoldChange) > lfc_thresh & padj < padj_thresh) %>%
                arrange(desc(abs(log2FoldChange)))
  deg_up   <- deg_all %>% filter(log2FoldChange >  lfc_thresh)
  deg_down <- deg_all %>% filter(log2FoldChange < -lfc_thresh)

  cat(sprintf("\nResults:\n"))
  cat(sprintf("  Total genes tested:        %d\n",   nrow(res_df)))
  cat(sprintf("  Upregulated in %-20s %d\n",   paste0(contrast_info$numerator, ":"), nrow(deg_up)))
  cat(sprintf("  Downregulated in %-19s %d\n", paste0(contrast_info$numerator, ":"), nrow(deg_down)))
  cat(sprintf("  Total DEGs:                %d\n",   nrow(deg_all)))

  # Save outputs
  contrast_dir <- file.path(output_dir, contrast_info$name)
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

  write.csv(res_df,   file.path(contrast_dir, paste0("all_genes_",        contrast_info$name, ".csv")), row.names = FALSE)
  write.csv(deg_all,  file.path(contrast_dir, paste0("DEGs_",             contrast_info$name, ".csv")), row.names = FALSE)
  write.csv(deg_up,   file.path(contrast_dir, paste0("DEGs_upregulated_", contrast_info$name, ".csv")), row.names = FALSE)
  write.csv(deg_down, file.path(contrast_dir, paste0("DEGs_downregulated_",contrast_info$name, ".csv")), row.names = FALSE)

  cat(sprintf("✓ Results saved to: %s/\n", contrast_dir))

  # Volcano plot (minimalistic, same palette as PCA), saved as png/pdf/svg
  v <- make_volcano_plot(res_df, contrast_info, lfc_thresh, padj_thresh, colors)
  volcano_basepath <- file.path(contrast_dir, paste0("volcano_", contrast_info$name))
  save_plot_multi(v, volcano_basepath, width = 16, height = 14)
  cat(sprintf("✓ Volcano plots saved (png/pdf/svg): %s.*\n", volcano_basepath))

  return(list(
    contrast = contrast_info,
    results  = res_df,
    degs     = deg_all,
    up       = deg_up,
    down     = deg_down,
    volcano  = v
  ))
}
# Run all contrasts
contrast_results <- lapply(contrasts_to_run, function(contrast) {
  tryCatch({
    run_contrast(contrast, dds_filtered, output_dir, lfc_threshold, padj_threshold, colors)
  }, error = function(e) {
    cat(sprintf("✗ ERROR processing %s: %s\n", contrast$name, e$message))
    return(NULL)
  })
})

# ----------------------------------------------------------------------------
# 6. SAVE DESEQ2 OBJECTS AND SUMMARY
# ----------------------------------------------------------------------------

cat("\n[6/6] Saving DESeq2 objects and creating summary...\n")

save(dds_filtered, vst, coldata, contrast_results,
     file = file.path(output_dir, "deseq2_objects.RData"))

cat(sprintf("✓ DESeq2 objects saved\n"))

# Create summary table
summary_df <- data.frame(
  Contrast = character(),
  Description = character(),
  Upregulated = integer(),
  Downregulated = integer(),
  Total_DEGs = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(contrast_results)) {
  if (!is.null(contrast_results[[i]])) {
    summary_df <- rbind(summary_df, data.frame(
      Contrast = contrast_results[[i]]$contrast$name,
      Description = contrast_results[[i]]$contrast$description,
      Upregulated = nrow(contrast_results[[i]]$up),
      Downregulated = nrow(contrast_results[[i]]$down),
      Total_DEGs = nrow(contrast_results[[i]]$degs)
    ))
  }
}

write.csv(summary_df, 
          file.path(output_dir, "DEG_summary.csv"), 
          row.names = FALSE)

# ----------------------------------------------------------------------------
# FINAL SUMMARY
# ----------------------------------------------------------------------------

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=================================================================\n\n")
cat(sprintf("Output directory: %s\n\n", output_dir))

cat("Summary of DEG counts (|LFC|>1, padj<0.05):\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
print(summary_df, row.names = FALSE)

cat("\n=================================================================\n")
