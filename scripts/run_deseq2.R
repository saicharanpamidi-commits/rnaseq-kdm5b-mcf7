#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Full Analysis: shKDM5B vs shSCR (MCF-7 breast cancer)
# Counts: featureCounts reverse-stranded (-s 2), GENCODE v44
# =============================================================================

# Use Cairo-based bitmap rendering (no X11 needed on HPC)
options(bitmapType = "cairo")

# ─── 0. Package installation ─────────────────────────────────────────────────
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
    if (bioc) {
      BiocManager::install(pkg, lib = user_lib, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
    }
  }
}

ensure_pkg("BiocManager",      bioc = FALSE)
ensure_pkg("DESeq2",           bioc = TRUE)
ensure_pkg("apeglm",           bioc = TRUE)
ensure_pkg("org.Hs.eg.db",     bioc = TRUE)
ensure_pkg("AnnotationDbi",    bioc = TRUE)
ensure_pkg("clusterProfiler",  bioc = TRUE)
ensure_pkg("ReactomePA",       bioc = TRUE)
ensure_pkg("enrichplot",       bioc = TRUE)
ensure_pkg("ggplot2",          bioc = FALSE)
ensure_pkg("ggrepel",          bioc = FALSE)
ensure_pkg("pheatmap",         bioc = FALSE)
ensure_pkg("RColorBrewer",     bioc = FALSE)
ensure_pkg("dplyr",            bioc = FALSE)
ensure_pkg("tibble",           bioc = FALSE)
ensure_pkg("msigdbr",          bioc = FALSE)

# ─── 1. Libraries ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(msigdbr)
})

# ─── 2. Paths ─────────────────────────────────────────────────────────────────
COUNTS_FILE <- "/N/slate/bkota/slate_project/results/featurecounts/gene_counts.txt"
OUT_DIR     <- "/N/slate/bkota/slate_project/deseq2_results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== DESeq2: shKDM5B vs shSCR | MCF-7 ===\n\n")

# ─── 3. Load featureCounts ────────────────────────────────────────────────────
cat("Loading counts...\n")
raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                  skip = 1, stringsAsFactors = FALSE)

count_mat <- as.matrix(raw[, 7:ncol(raw)])
rownames(count_mat) <- raw$Geneid
# Strip path and .sorted.bam suffix from column names
colnames(count_mat) <- gsub(".*/", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.sorted\\.bam$", "", colnames(count_mat))

cat("Genes:", nrow(count_mat), "| Samples:", paste(colnames(count_mat), collapse = ", "), "\n\n")

# ─── 4. Sample metadata (shSCR = reference) ───────────────────────────────────
sample_info <- data.frame(
  sample    = c("SRR31992838","SRR31992839","SRR31992840",
                "SRR31992841","SRR31992842","SRR31992843"),
  condition = factor(c("shSCR","shSCR","shSCR",
                        "shKDM5B","shKDM5B","shKDM5B"),
                     levels = c("shSCR","shKDM5B")),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample

# Re-order columns to match sample_info
count_mat <- count_mat[, rownames(sample_info)]
stopifnot(all(colnames(count_mat) == rownames(sample_info)))

# ─── 5. DESeq2 ────────────────────────────────────────────────────────────────
cat("Building DESeqDataSet...\n")
dds  <- DESeqDataSetFromMatrix(countData = count_mat,
                                colData   = sample_info,
                                design    = ~ condition)
# Pre-filter: keep genes with >=10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes after low-count filter:", nrow(dds), "\n")

cat("Running DESeq2...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","shKDM5B","shSCR"), alpha = 0.05)
cat("\nResults summary:\n"); summary(res)

# LFC shrinkage with apeglm
cat("\nApplying apeglm LFC shrinkage...\n")
res_shrunk <- lfcShrink(dds, coef = "condition_shKDM5B_vs_shSCR", type = "apeglm")

# ─── 6. Gene symbol mapping ───────────────────────────────────────────────────
cat("\nMapping ENSEMBL IDs → gene symbols...\n")
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  mutate(ensembl_clean = sub("\\..*$", "", gene_id))  # strip version suffix

res_df$gene_symbol <- mapIds(org.Hs.eg.db,
  keys = res_df$ensembl_clean, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first")

res_df$entrez_id <- mapIds(org.Hs.eg.db,
  keys = res_df$ensembl_clean, column = "ENTREZID",
  keytype = "ENSEMBL", multiVals = "first")

# Add Wald statistic from un-shrunk results (used later for GSEA ranking)
res_df$stat <- results(dds, contrast = c("condition","shKDM5B","shSCR"))$stat[
  match(res_df$gene_id, rownames(results(dds, contrast = c("condition","shKDM5B","shSCR"))))]

res_df <- res_df %>%
  arrange(padj) %>%
  dplyr::select(gene_id, gene_symbol, entrez_id, baseMean,
                log2FoldChange, lfcSE, stat, pvalue, padj)

# ─── 7. Save annotated results ────────────────────────────────────────────────
write.csv(res_df,
          file.path(OUT_DIR, "deseq2_results_annotated.csv"),
          row.names = FALSE)

sig_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange >  0, na.rm = TRUE)
sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange <  0, na.rm = TRUE)
cat(sprintf("\nDE summary: %d up | %d down | %d total (padj<0.05)\n\n",
            sig_up, sig_down, sig_up + sig_down))

cat("Top 10 upregulated:\n")
print(res_df %>% filter(log2FoldChange > 0) %>%
      head(10) %>% dplyr::select(gene_symbol, log2FoldChange, padj))
cat("\nTop 10 downregulated:\n")
print(res_df %>% filter(log2FoldChange < 0) %>%
      head(10) %>% dplyr::select(gene_symbol, log2FoldChange, padj))

# ─── 8. PCA plot (VST, colored by condition) ──────────────────────────────────
cat("\nGenerating PCA plot...\n")
vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.2, show.legend = FALSE, box.padding = 0.4) +
  scale_color_manual(values = c("shSCR" = "#2166AC", "shKDM5B" = "#D6604D")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA — VST Normalized", subtitle = "shKDM5B vs shSCR | MCF-7") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "pca_plot.png"), pca_p, width = 7, height = 6, dpi = 150)
cat("  pca_plot.png saved.\n")

# ─── 9. MA plot (significant in red) ─────────────────────────────────────────
cat("Generating MA plot...\n")
ma_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  mutate(sig = padj < 0.05)

ma_p <- ggplot(ma_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D6604D"),
                     labels = c("NS", "padj<0.05")) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", linewidth = 0.4) +
  xlab("log10(Mean Normalized Counts + 1)") +
  ylab("log2 Fold Change (shKDM5B / shSCR)") +
  ggtitle("MA Plot — shKDM5B vs shSCR",
          subtitle = "Red = significant (padj<0.05) | dotted = ±1 LFC") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "ma_plot.png"), ma_p, width = 8, height = 6, dpi = 150)
cat("  ma_plot.png saved.\n")

# ─── 10. Volcano plot (padj<0.05 & |LFC|>1 labeled) ─────────────────────────
cat("Generating volcano plot...\n")
vol_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj + 1e-300),
    label    = ifelse(is.na(gene_symbol), sub("\\..*$", "", gene_id), gene_symbol),
    category = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                               ~ "NS"
    )
  )

top_lbl <- vol_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
  arrange(padj) %>%
  head(25)

vol_p <- ggplot(vol_df, aes(x = log2FoldChange, y = neg_log10_padj, color = category)) +
  geom_point(size = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#D6604D", "Down" = "#2166AC", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  geom_text_repel(data = top_lbl, aes(label = label),
                  size = 2.8, color = "black",
                  box.padding = 0.35, max.overlaps = 25,
                  show.legend = FALSE) +
  xlab("log2 Fold Change (shKDM5B / shSCR)") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Volcano Plot — shKDM5B vs shSCR",
          subtitle = sprintf("MCF-7 KDM5B knockdown  |  Up: %d  Down: %d  (padj<0.05, |LFC|>1)",
                             sum(vol_df$category == "Up"),
                             sum(vol_df$category == "Down"))) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "volcano_plot.png"), vol_p, width = 9, height = 7, dpi = 150)
cat("  volcano_plot.png saved.\n")

# ─── 11. Heatmap — top 50 DE genes, Z-score, gene symbol labels ──────────────
cat("Generating heatmap...\n")
top50_ids <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  head(50) %>%
  pull(gene_id)

if (length(top50_ids) >= 2) {
  vst_mat    <- assay(vsd)[top50_ids, ]
  vst_scaled <- t(scale(t(vst_mat)))

  # Use gene symbols as row labels (fall back to ENSEMBL if NA)
  sym_map <- res_df %>%
    filter(gene_id %in% top50_ids) %>%
    dplyr::select(gene_id, gene_symbol) %>%
    mutate(label = ifelse(is.na(gene_symbol), sub("\\..*$", "", gene_id), gene_symbol))
  rownames(vst_scaled) <- sym_map$label[match(rownames(vst_scaled), sym_map$gene_id)]

  col_ann    <- data.frame(Condition = sample_info$condition,
                           row.names = rownames(sample_info))
  ann_colors <- list(Condition = c("shSCR" = "#2166AC", "shKDM5B" = "#D6604D"))

  pheatmap(vst_scaled,
           annotation_col    = col_ann,
           annotation_colors = ann_colors,
           color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 8, fontsize_col = 10,
           main = "Top 50 DE Genes — shKDM5B vs shSCR\n(Z-score of VST counts | MCF-7)",
           filename = file.path(OUT_DIR, "heatmap_top50.png"),
           width = 8.5, height = 12)
  cat("  heatmap_top50.png saved.\n")
} else {
  cat("  WARNING: fewer than 2 significant genes — heatmap skipped.\n")
}

# =============================================================================
# PATHWAY ANALYSIS
# =============================================================================

sig_genes   <- res_df %>% filter(padj < 0.05, !is.na(entrez_id))
up_entrez   <- sig_genes %>% filter(log2FoldChange >  1) %>% pull(entrez_id)
down_entrez <- sig_genes %>% filter(log2FoldChange < -1) %>% pull(entrez_id)
all_entrez  <- sig_genes %>% pull(entrez_id)
bg_entrez   <- res_df %>% filter(!is.na(entrez_id)) %>% pull(entrez_id)

cat(sprintf("\nPathway inputs: %d up | %d down | %d total DE | %d background\n\n",
            length(up_entrez), length(down_entrez),
            length(all_entrez), length(bg_entrez)))

# ─── 12. GO Biological Process ORA ───────────────────────────────────────────
cat("Running GO BP ORA...\n")
go_bp <- enrichGO(
  gene          = all_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (!is.null(go_bp) && nrow(go_bp) > 0) {
  write.csv(as.data.frame(go_bp),
            file.path(OUT_DIR, "go_bp_results.csv"), row.names = FALSE)
  p <- dotplot(go_bp, showCategory = 20,
               title = "GO Biological Process — All DE Genes", font.size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12))
  ggsave(file.path(OUT_DIR, "go_bp_dotplot.png"), p, width = 10, height = 10, dpi = 150)
  cat(sprintf("  GO BP: %d significant terms.\n", nrow(go_bp)))
} else {
  cat("  WARNING: No significant GO BP terms found.\n")
}

# ─── 13. KEGG ORA ─────────────────────────────────────────────────────────────
cat("Running KEGG ORA...\n")
kegg_res <- enrichKEGG(
  gene          = all_entrez,
  universe      = bg_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg_res),
            file.path(OUT_DIR, "kegg_results.csv"), row.names = FALSE)
  p <- dotplot(kegg_res, showCategory = 20,
               title = "KEGG Pathways — All DE Genes", font.size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12))
  ggsave(file.path(OUT_DIR, "kegg_dotplot.png"), p, width = 10, height = 9, dpi = 150)
  cat(sprintf("  KEGG: %d significant pathways.\n", nrow(kegg_res)))
} else {
  cat("  WARNING: No significant KEGG pathways found.\n")
}

# ─── 14. Reactome ORA ─────────────────────────────────────────────────────────
cat("Running Reactome ORA...\n")
reactome_res <- enrichPathway(
  gene          = all_entrez,
  universe      = bg_entrez,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

if (!is.null(reactome_res) && nrow(reactome_res) > 0) {
  write.csv(as.data.frame(reactome_res),
            file.path(OUT_DIR, "reactome_results.csv"), row.names = FALSE)
  p <- dotplot(reactome_res, showCategory = 20,
               title = "Reactome Pathways — All DE Genes", font.size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12))
  ggsave(file.path(OUT_DIR, "reactome_dotplot.png"), p, width = 12, height = 10, dpi = 150)
  cat(sprintf("  Reactome: %d significant pathways.\n", nrow(reactome_res)))
} else {
  cat("  WARNING: No significant Reactome pathways found.\n")
}

# ─── 15. GSEA — MSigDB Hallmarks ──────────────────────────────────────────────
cat("Running GSEA (MSigDB Hallmarks)...\n")

# Use Wald test statistic as ranking metric — continuous, no ties
ranked_df <- res_df %>%
  filter(!is.na(stat), !is.na(entrez_id)) %>%
  arrange(desc(stat)) %>%
  filter(!duplicated(entrez_id))

gene_list        <- ranked_df$stat
names(gene_list) <- ranked_df$entrez_id

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

gsea_res <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = hallmarks,
  pvalueCutoff  = 0.2,    # BH-adjusted; use 0.2 for hallmarks discovery
  pAdjustMethod = "BH",
  minGSSize     = 15,
  maxGSSize     = 500,
  eps           = 0,
  seed          = 42,
  verbose       = FALSE
)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  write.csv(as.data.frame(gsea_res),
            file.path(OUT_DIR, "gsea_hallmarks_results.csv"), row.names = FALSE)

  p <- dotplot(gsea_res, showCategory = 20, split = ".sign",
               title = "GSEA — MSigDB Hallmarks") +
    facet_grid(. ~ .sign) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 9))
  ggsave(file.path(OUT_DIR, "gsea_dotplot.png"), p, width = 14, height = 10, dpi = 150)

  cat(sprintf("  GSEA: %d significant hallmarks.\n", nrow(gsea_res)))
  cat("\nTop activated hallmarks:\n")
  print(as.data.frame(gsea_res) %>% filter(NES > 0) %>%
        arrange(p.adjust) %>% head(5) %>%
        dplyr::select(ID, NES, p.adjust))
  cat("\nTop suppressed hallmarks:\n")
  print(as.data.frame(gsea_res) %>% filter(NES < 0) %>%
        arrange(p.adjust) %>% head(5) %>%
        dplyr::select(ID, NES, p.adjust))
} else {
  cat("  WARNING: No significant GSEA hallmarks at padj<0.05.\n")
}

# ─── 16. Summary ──────────────────────────────────────────────────────────────
cat("\n=== COMPLETE ===\n")
cat("Output directory:", OUT_DIR, "\n")
cat(sprintf("DE: %d up | %d down | %d total (padj<0.05)\n", sig_up, sig_down, sig_up + sig_down))
cat("\nFiles saved:\n")
for (f in list.files(OUT_DIR, full.names = FALSE)) cat(" ", f, "\n")
