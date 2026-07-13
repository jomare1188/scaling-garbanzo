# =============================================================================
# KEGG and GO Enrichment Analysis Pipeline - eggNOG Format
# =============================================================================
# Description: Performs KEGG enrichment analysis and creates GO-KEGG 
#              interaction networks for differentially expressed genes
#              Modified to work with eggNOG annotation files
# =============================================================================

# Load required libraries
library(readr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)
library(svglite)

# =============================================================================
# 1. DATA LOADING FUNCTIONS
# =============================================================================

#' Read and process eggNOG annotations
#'
#' @param annot_file Path to eggNOG annotation file
#' @return List containing KO vector and gene-KO mapping dataframe
read_eggnog_annotations <- function(annot_file) {
  
  # Read eggNOG file (skip comment lines starting with #)
  annot_df <- read_tsv(annot_file, comment = "#", show_col_types = FALSE)

  # Remove .t(number) suffix from query column (e.g., g1.t1 -> g1)
  annot_df <- annot_df %>% mutate(query = sub("\\.t\\d+$", "", query)) 

  # Extract gene-KO mapping
  kos_genes <- annot_df %>%
    select(gene = query, KO = KEGG_ko) %>%
    filter(!is.na(KO) & KO != "" & KO != "-") %>%
    # Split multiple KOs (comma-separated in eggNOG)
    separate_rows(KO, sep = ",") %>%
    # Remove 'ko:' prefix if present
    mutate(KO = sub("^ko:", "", KO)) %>%
    mutate(KO = trimws(KO)) %>%
    filter(KO != "" & KO != "-") %>%
    distinct()
  
  kos <- unique(kos_genes$KO)
  
  message("  Loaded ", nrow(annot_df), " gene annotations")
  message("  Found ", nrow(kos_genes), " gene-KO mappings")
  message("  Unique KOs: ", length(kos))
  
  return(list(
    kos = kos,
    kos_genes = kos_genes  # This is the gene-to-KO mapping needed later
  ))
}


#' Read differentially expressed genes
#'
#' @param up_file Path to up-regulated genes file
#' @param down_file Path to down-regulated genes file
#' @return List with up and down regulated gene vectors
read_deg_files <- function(up_file, down_file) {
  
  up_genes <- read.table(up_file, header = T, sep = ",")$gene_id
 # up_genes <- sub("\\.v2\\.1$", "", up_genes)
  
  down_genes <- read.table(down_file, header = T, sep = ",")$gene_id
 # down_genes <- sub("\\.v2\\.1$", "", down_genes)
  
  return(list(
    up = up_genes,
    down = down_genes
  ))
}


#' Read GO annotations from eggNOG file
#'
#' @param annot_file Path to eggNOG annotation file
#' @param go_results_file Path to GO enrichment results
#' @return List with GO annotations and results
read_go_data <- function(annot_file, go_results_file) {
  
  # Read eggNOG file
  annot_df <- read_tsv(annot_file, comment = "#", show_col_types = FALSE)
  
  annot_df <- annot_df %>% mutate(query = sub("\\.t\\d+$", "", query))

  # Extract gene-GO mapping
  all_go <- annot_df %>%
    select(gene = query, GO = GOs) %>%
    filter(!is.na(GO) & GO != "" & GO != "-") %>%
    # Split multiple GO terms (comma-separated in eggNOG)
    separate_rows(GO, sep = ",") %>%
    mutate(GO = trimws(GO)) %>%
    filter(GO != "" & GO != "-") %>%
    distinct()
  
  message("  Found ", nrow(all_go), " gene-GO mappings")
  
  # Read GO enrichment results
  go_results <- read.table(go_results_file, header = TRUE, sep = ",")
  
  return(list(
    annotations = all_go,
    results = go_results
  ))
}


# =============================================================================
# 2. ENRICHMENT ANALYSIS FUNCTIONS
# =============================================================================

#' Perform KEGG enrichment analysis
#'
#' @param gene_kos Vector of KOs for genes of interest
#' @param universe_kos Vector of all KOs (background)
#' @param pvalue_cutoff P-value cutoff (default: 0.05)
#' @param qvalue_cutoff Q-value cutoff (default: 0.05)
#' @return enrichResult object from clusterProfiler
perform_kegg_enrichment <- function(gene_kos, universe_kos, 
                                   pvalue_cutoff = 0.05, 
                                   qvalue_cutoff = 0.05) {
  
  result <- enrichKEGG(
    gene = gene_kos,
    universe = universe_kos,
    organism = "ko",
    pAdjustMethod = "fdr",
    keyType = "kegg",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )
  
  return(result)
}


# =============================================================================
# 3. NETWORK PREPARATION FUNCTIONS
# =============================================================================

#' Prepare KEGG edge list for network
#'
#' @param kegg_result enrichResult object from KEGG enrichment
#' @param kos_genes Dataframe with gene-KO mapping
#' @param dea_data DEA results with log2FoldChange
#' @param padj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @return Dataframe with gene-term edges and log2FC
prepare_kegg_edges <- function(kegg_result, kos_genes, dea_data, 
                               padj_cutoff = 0.05) {
  
  # Extract significant KEGG terms
  kegg_df <- kegg_result@result %>%
    filter(p.adjust < padj_cutoff) %>%
    select(Description, geneID)
  
  # Create edge list
  edge_list <- kegg_df %>%
    mutate(geneID = strsplit(geneID, "/")) %>%
    unnest(geneID) %>%
    mutate(geneID = trimws(geneID)) %>%
    left_join(kos_genes, by = c("geneID" = "KO")) %>%
    select(gene, Description, geneID)
  
  # Add log2 fold change
  edge_list_with_lfc <- edge_list %>%
    left_join(dea_data %>% select(gene, log2FoldChange), by = "gene") %>%
    filter(!is.na(log2FoldChange))
  
  return(edge_list_with_lfc)
}


#' Prepare GO edge list for network
#'
#' @param go_results GO enrichment results
#' @param go_annotations All gene-GO annotations
#' @param dea_data DEA results with log2FoldChange
#' @param top_n Number of top GO terms to include (default: 20)
#' @return Dataframe with gene-GO term edges
prepare_go_edges <- function(go_results, go_annotations, dea_data, 
                             top_n = 20) {
  
  # Process GO results
  go_genes <- go_results %>%
    left_join(go_annotations, by = c("GO.ID" = "GO")) %>%
    group_by(GO.ID, Term, Annotated, Significant, Expected, Classic, p.adj) %>%
    summarise(gene = paste(unique(gene), collapse = "/"), .groups = "drop") %>%
    arrange(p.adj) %>%
    head(top_n)
  
  # Create edge list
  edge_list <- go_genes %>%
    separate_rows(gene, sep = "/") %>%
    select(gene, Term)
  
  # Add log2 fold change
  edge_list_with_lfc <- edge_list %>%
    left_join(dea_data %>% select(gene, log2FoldChange), by = "gene") %>%
    filter(!is.na(log2FoldChange))
  
  return(edge_list_with_lfc)
}


#' Combine GO and KEGG edge lists
#'
#' @param go_edges GO edge list
#' @param kegg_edges KEGG edge list
#' @param dea_data DEA results with log2FoldChange
#' @return Combined edge list with gene-term connections
combine_edge_lists <- function(go_edges, kegg_edges, dea_data) {
  
  # Format and label edges
  go_formatted <- go_edges %>%
    select(gene, Term) %>%
    mutate(Class = "GO") %>%
    rename(Description = Term)
  
  kegg_formatted <- kegg_edges %>%
    select(gene, Description) %>%
    mutate(Class = "KEGG")
  
  # Combine
  combined <- rbind(go_formatted, kegg_formatted) %>%
    filter(!is.na(gene) & gene != "NA") %>%
    left_join(dea_data %>% select(gene, log2FoldChange), by = "gene")
  
  return(combined)
}


# =============================================================================
# 4. NETWORK VISUALIZATION FUNCTIONS
# =============================================================================

#' Create and configure igraph network
#'
#' @param edge_list Combined edge list
#' @param tf_genes Optional vector of transcription factor gene IDs
#' @return Configured igraph object
create_network_graph <- function(edge_list, tf_genes = NULL) {
  
  # Build graph
  g <- graph_from_data_frame(d = edge_list, directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  # Set node types
  V(g)$type <- ifelse(V(g)$name %in% edge_list$Description, 
                      "Description", "Gene")
  
  # Mark transcription factors if provided
  if (!is.null(tf_genes)) {
    V(g)$type[V(g)$name %in% tf_genes] <- "Transcription Factor"
  }
  
  # Assign classes to description nodes
  desc_class <- edge_list %>% distinct(Description, Class)
  V(g)$Class <- desc_class$Class[match(V(g)$name, desc_class$Description)]
  
  # Set node categories for coloring
  V(g)$Category <- case_when(
    V(g)$type == "Transcription Factor" ~ "Transcription Factor",
    V(g)$type == "Gene" ~ "Gene",
    V(g)$Class == "GO" ~ "GO term",
    V(g)$Class == "KEGG" ~ "KEGG pathway",
    TRUE ~ NA_character_
  )
  
  # Set labels (hide gene labels, show term labels)
  V(g)$label <- ifelse(V(g)$type == "Gene", "", V(g)$name)
  
  # Scale label size by node degree
  degrees <- degree(g)
  V(g)$label.cex <- rescale(degrees, to = c(0.6, 0.8))
  
  return(g)
}


#' Plot gene-term network
#'
#' @param graph igraph object
#' @param title Plot title
#' @param layout Layout algorithm (default: 'tree')
#' @return ggraph plot object
plot_gene_network <- function(graph, title = "Gene-Functional Term Network",
                              layout = "tree") {
  
  # Define color palette
  color_palette <- c(
    "Gene" = "gray80",
    "GO term" = "#66c2a5",
    "KEGG pathway" = "#fc8d62",
    "Transcription Factor" = "#FF0033"
  )
  
  # Get categories present in this graph
  categories_present <- levels(factor(V(graph)$Category))
  colors_to_use <- color_palette[names(color_palette) %in% categories_present]
  
  # Create plot
  p <- ggraph(graph, layout = layout) +
    geom_edge_link(alpha = 0.4, colour = "grey70") +
    geom_node_point(aes(color = Category), size = 4, show.legend = TRUE) +
    geom_node_text(aes(label = label), repel = TRUE, 
                   angle = 60, size = 3, color = "black") +
    scale_color_manual(name = "Node Type", values = colors_to_use) +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}


#' Save network plots and data
#'
#' @param plot ggraph plot object
#' @param graph igraph object
#' @param prefix Output file prefix
#' @param width Plot width in inches (default: 16)
#' @param height Plot height in inches (default: 9)
#' @param dpi Resolution for PNG (default: 300)
save_network_outputs <- function(plot, graph, prefix, 
                                width = 16, height = 9, dpi = 300) {
  
  # Save plots
  ggsave(paste0(prefix, ".svg"), plot, 
         width = width, height = height, bg = "white")
  ggsave(paste0(prefix, ".png"), plot, 
         width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paste0(prefix, ".pdf"), plot, 
         width = width, height = height)
  
  # Save edge list with metadata
  edges <- as_data_frame(graph, what = "edges") %>%
    mutate(
      source_Class = "gene",
      target_Class = V(graph)$Class[match(to, V(graph)$name)]
    )
  
  write.table(edges, paste0(prefix, "_edges.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Network outputs saved with prefix: ", prefix)
}


# =============================================================================
# 5. MAIN PIPELINE FUNCTION
# =============================================================================

#' Run complete enrichment and network analysis pipeline
#'
#' @param config List with all file paths and parameters
#' @return List with enrichment results and networks
run_enrichment_pipeline <- function(config) {
  
  message("=== Starting Enrichment Analysis Pipeline ===\n")
  
  # 1. Load data
  message("Loading eggNOG annotations...")
  annot_data <- read_eggnog_annotations(config$annot_file)
  
  message("\nLoading differentially expressed genes...")
  deg_data <- read_deg_files(config$up_file, config$down_file)
  
  message("\nLoading GO data...")
  go_up <- read_go_data(config$annot_file, config$go_up_file)
  go_down <- read_go_data(config$annot_file, config$go_down_file)
  
  # 2. Prepare enrichment data
  message("\nPreparing enrichment data...")
  
  # Extract KOs for up and down regulated genes
  up_kos <- annot_data$kos_genes %>%
    filter(gene %in% deg_data$up) %>%
    pull(KO) %>%
    unique()
  
  down_kos <- annot_data$kos_genes %>%
    filter(gene %in% deg_data$down) %>%
    pull(KO) %>%
    unique()
  
  universe_kos <- unique(annot_data$kos)
  
  message("  Up-regulated genes: ", length(deg_data$up))
  message("  Up-regulated KOs: ", length(up_kos))
  message("  Down-regulated genes: ", length(deg_data$down))
  message("  Down-regulated KOs: ", length(down_kos))
  message("  Universe KOs: ", length(universe_kos))
  
  # 3. Perform KEGG enrichment
  message("\nPerforming KEGG enrichment...")
  kegg_up <- perform_kegg_enrichment(up_kos, universe_kos)
  kegg_down <- perform_kegg_enrichment(down_kos, universe_kos)
  
  message("  Up-regulated: ", nrow(kegg_up@result[kegg_up@result$p.adjust < 0.05,]), " significant terms")
  message("  Down-regulated: ", nrow(kegg_down@result[kegg_down@result$p.adjust < 0.05,]), " significant terms")
  
  # 4. Create visualizations
  message("\nCreating dotplots...")
  print(dotplot(kegg_up, title = "KEGG Enrichment - Up-regulated"))
  print(dotplot(kegg_down, title = "KEGG Enrichment - Down-regulated"))
  
  # 5. Prepare and plot networks
  message("\n=== Building Networks ===")
  
  # Read DEA data
  dea_up <- read.table(config$up_file, sep = ",", header = TRUE)
  dea_up$gene <- dea_up$gene_id

#    tibble::rownames_to_column("gene") %>%
#    mutate(gene = sub("\\.v2\\.1$", "", gene))
  
  dea_down <- read.table(config$down_file, sep = ",", header = TRUE)
  dea_down$gene <- dea_down$gene_id
#    tibble::rownames_to_column("gene") %>%
#    mutate(gene = sub("\\.v2\\.1$", "", gene))
  
  # Up-regulated network
  message("\nBuilding up-regulated gene network...")
  kegg_edges_up <- prepare_kegg_edges(kegg_up, annot_data$kos_genes, dea_up)
  go_edges_up <- prepare_go_edges(go_up$results, go_up$annotations, dea_up)
  combined_up <- combine_edge_lists(go_edges_up, kegg_edges_up, dea_up)
  
  graph_up <- create_network_graph(combined_up, config$tf_genes)
  plot_up <- plot_gene_network(graph_up, "Up-Regulated Gene–Functional Term Network")
  
  save_network_outputs(plot_up, graph_up, "gene_network_up")
  
  # Down-regulated network
  message("\nBuilding down-regulated gene network...")
  kegg_edges_down <- prepare_kegg_edges(kegg_down, annot_data$kos_genes, dea_down)
  go_edges_down <- prepare_go_edges(go_down$results, go_down$annotations, dea_down)
  combined_down <- combine_edge_lists(go_edges_down, kegg_edges_down, dea_down)
  
  graph_down <- create_network_graph(combined_down, config$tf_genes)
  plot_down <- plot_gene_network(graph_down, "Down-Regulated Gene–Functional Term Network")
  
  save_network_outputs(plot_down, graph_down, "gene_network_down")
  
  message("\n=== Pipeline Complete ===")
  
  # Return results
  return(list(
    kegg_up = kegg_up,
    kegg_down = kegg_down,
    networks = list(
      up = list(graph = graph_up, plot = plot_up),
      down = list(graph = graph_down, plot = plot_down)
    )
  ))
}


# =============================================================================
# 6. EXAMPLE USAGE
# =============================================================================

# Configure file paths and parameters
config <- list(
  # eggNOG annotation file
  annot_file = "/home/joao/RNAseq_RNAi/genome/files_bianca/annotation/eggnog.emapper.annotations",
  
  # DEG files
  up_file = "/home/joao/RNAseq_RNAi/rnaseq/test1/star_salmon/deseq2_qc/in_planta_IAC66_6_vs_in_vitro/DEGs_upregulated_in_planta_IAC66_6_vs_in_vitro.csv",
  down_file = "/home/joao/RNAseq_RNAi/rnaseq/test1/star_salmon/deseq2_qc/in_planta_IAC66_6_vs_in_vitro/DEGs_downregulated_in_planta_IAC66_6_vs_in_vitro.csv",
  
  # GO files
  go_up_file = "/home/joao/RNAseq_RNAi/rnaseq/test1/star_salmon/deseq2_qc/in_planta_IAC66_6_vs_in_vitro/GO_enrichment/GO_BP_upregulated.csv",
  go_down_file = "/home/joao/RNAseq_RNAi/rnaseq/test1/star_salmon/deseq2_qc/in_planta_IAC66_6_vs_in_vitro/GO_enrichment/GO_BP_downregulated.csv",
  
  # Optional: Transcription factors to highlight
  tf_genes = NULL  # Or: read.table("tf_in_degs.txt")$V1
)

# Run the pipeline
results <- run_enrichment_pipeline(config)
