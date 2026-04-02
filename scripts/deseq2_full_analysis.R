#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Full Analysis Pipeline
# Experiment: KDM5B knockdown vs scramble control in MCF-7
# Outputs:
#   - Results table with gene symbols
#   - PCA plot
#   - MA plot
#   - Volcano plot (with gene symbols)
#   - Heatmap top 50 DE genes
#   - ORA: GO Biological Process, KEGG, Reactome
#   - GSEA: MSigDB Hallmark gene sets
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Install missing packages (run once if needed)
# -----------------------------------------------------------------------------
# Uncomment and run this block once if packages are not installed:
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "apeglm", "org.Hs.eg.db",
#                        "clusterProfiler", "ReactomePA", "enrichplot"))
# install.packages(c("ggplot2", "pheatmap", "RColorBrewer",
#                    "dplyr", "tibble", "msigdbr", "ggrepel"))

# -----------------------------------------------------------------------------
# 1. Load Libraries
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(msigdbr)
})

# -----------------------------------------------------------------------------
# 2. Paths
# -----------------------------------------------------------------------------
COUNTS_FILE <- "/N/slate/bkota/slate_project/results/counts/gene_counts.txt"
OUTPUT_DIR  <- "/N/slate/bkota/slate_project/results/deseq2"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "pathway"), recursive = TRUE, showWarnings = FALSE)

cat("=== DESeq2 Full Analysis: shKDM5B vs shSCR (MCF-7) ===\n\n")

# -----------------------------------------------------------------------------
# 3. Load and Parse featureCounts Output
# -----------------------------------------------------------------------------
cat("Loading counts...\n")
raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                  skip = 1, stringsAsFactors = FALSE)

count_matrix <- raw[, 7:ncol(raw)]
rownames(count_matrix) <- raw$Geneid
colnames(count_matrix) <- gsub(".*/", "", colnames(count_matrix))
colnames(count_matrix) <- gsub("\\.sorted\\.bam$", "", colnames(count_matrix))

cat("Genes:", nrow(count_matrix), "| Samples:", paste(colnames(count_matrix), collapse=", "), "\n\n")

# -----------------------------------------------------------------------------
# 4. Sample Metadata
# -----------------------------------------------------------------------------
sample_info <- data.frame(
  sample    = c("SRR31992838","SRR31992839","SRR31992840",
                "SRR31992841","SRR31992842","SRR31992843"),
  condition = c("shSCR","shSCR","shSCR",
                "shKDM5B","shKDM5B","shKDM5B"),
  stringsAsFactors = FALSE
)
sample_info$condition <- factor(sample_info$condition, levels = c("shSCR","shKDM5B"))
rownames(sample_info) <- sample_info$sample
count_matrix <- count_matrix[, rownames(sample_info)]
stopifnot(all(colnames(count_matrix) == rownames(sample_info)))

# -----------------------------------------------------------------------------
# 5. DESeq2
# -----------------------------------------------------------------------------
cat("Running DESeq2...\n")
dds  <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData   = sample_info,
                                design    = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","shKDM5B","shSCR"), alpha = 0.05)
cat("\nResults summary:\n"); summary(res)

# LFC shrinkage
cat("\nApplying apeglm LFC shrinkage...\n")
res_shrunk <- lfcShrink(dds, coef = "condition_shKDM5B_vs_shSCR", type = "apeglm")

# -----------------------------------------------------------------------------
# 6. Gene Symbol Mapping
# -----------------------------------------------------------------------------
cat("\nMapping ENSEMBL IDs to gene symbols...\n")
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    ensembl_clean = gsub("\\..*$", "", gene_id)  # strip version e.g. .3
  )

# Map to gene symbols, Entrez IDs, and gene names
res_df$gene_symbol <- mapIds(org.Hs.eg.db,
                              keys    = res_df$ensembl_clean,
                              column  = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")

res_df$entrez_id   <- mapIds(org.Hs.eg.db,
                              keys    = res_df$ensembl_clean,
                              column  = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")

res_df$gene_name   <- mapIds(org.Hs.eg.db,
                              keys    = res_df$ensembl_clean,
                              column  = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Final results table — ordered by padj
res_df <- res_df %>%
  arrange(padj) %>%
  mutate(
    significance = case_when(
      padj < 0.001 & abs(log2FoldChange) >= 2 ~ "padj<0.001 & |LFC|>=2",
      padj < 0.05  & abs(log2FoldChange) >= 1 ~ "padj<0.05 & |LFC|>=1",
      padj < 0.05                             ~ "padj<0.05",
      TRUE                                    ~ "NS"
    )
  ) %>%
  dplyr::select(gene_id, ensembl_clean, gene_symbol, gene_name, entrez_id,
                baseMean, log2FoldChange, lfcSE, pvalue, padj, significance)

write.csv(res_df,
          file.path(OUTPUT_DIR, "deseq2_results_annotated.csv"),
          row.names = FALSE)

sig_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange >  0, na.rm=TRUE)
sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange <  0, na.rm=TRUE)
cat(sprintf("\nUpregulated   (padj<0.05): %d\n", sig_up))
cat(sprintf("Downregulated (padj<0.05): %d\n", sig_down))
cat(sprintf("Total DE:                  %d\n\n", sig_up + sig_down))

# Preview top hits
cat("Top 10 upregulated genes:\n")
print(res_df %>% filter(log2FoldChange > 0, !is.na(padj)) %>%
      head(10) %>% dplyr::select(gene_symbol, log2FoldChange, padj))
cat("\nTop 10 downregulated genes:\n")
print(res_df %>% filter(log2FoldChange < 0, !is.na(padj)) %>%
      head(10) %>% dplyr::select(gene_symbol, log2FoldChange, padj))

# -----------------------------------------------------------------------------
# 7. PCA Plot
# -----------------------------------------------------------------------------
cat("\nGenerating PCA plot...\n")
vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x=PC1, y=PC2, color=condition, label=name)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(size=3.2, show.legend=FALSE, box.padding=0.4) +
  scale_color_manual(values=c("shSCR"="#2166AC","shKDM5B"="#D6604D")) +
  xlab(paste0("PC1: ",pct_var[1],"% variance")) +
  ylab(paste0("PC2: ",pct_var[2],"% variance")) +
  ggtitle("PCA — Sample Clustering",
          subtitle="shKDM5B vs shSCR | MCF-7 | VST normalized") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.title=element_blank(),
        panel.grid.minor=element_blank())

ggsave(file.path(OUTPUT_DIR,"pca_plot.pdf"), pca_plot, width=7, height=6)
ggsave(file.path(OUTPUT_DIR,"pca_plot.png"), pca_plot, width=7, height=6, dpi=150)

# -----------------------------------------------------------------------------
# 8. MA Plot
# -----------------------------------------------------------------------------
cat("Generating MA plot...\n")
ma_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  mutate(significant = padj < 0.05)

ma_plot <- ggplot(ma_df, aes(x=log10(baseMean+1), y=log2FoldChange, color=significant)) +
  geom_point(size=0.6, alpha=0.5) +
  scale_color_manual(values=c("FALSE"="grey60","TRUE"="#D6604D"),
                     labels=c("NS","padj<0.05")) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_hline(yintercept=c(-1,1), linetype="dotted", color="steelblue", linewidth=0.4) +
  xlab("log10(Mean Normalized Counts + 1)") +
  ylab("log2 Fold Change (shKDM5B / shSCR)") +
  ggtitle("MA Plot — shKDM5B vs shSCR",
          subtitle="MCF-7 | Red = significant (padj<0.05) | dotted = ±1 LFC") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.title=element_blank(),
        panel.grid.minor=element_blank())

ggsave(file.path(OUTPUT_DIR,"ma_plot.pdf"), ma_plot, width=8, height=6)
ggsave(file.path(OUTPUT_DIR,"ma_plot.png"), ma_plot, width=8, height=6, dpi=150)

# -----------------------------------------------------------------------------
# 9. Volcano Plot — with gene symbols
# -----------------------------------------------------------------------------
cat("Generating volcano plot...\n")
vol_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj + 1e-300),
    label          = ifelse(is.na(gene_symbol), ensembl_clean, gene_symbol),
    category       = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                               ~ "NS"
    )
  )

top_label <- vol_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
  arrange(padj) %>%
  head(20)

vol_plot <- ggplot(vol_df, aes(x=log2FoldChange, y=neg_log10_padj, color=category)) +
  geom_point(size=0.7, alpha=0.6) +
  scale_color_manual(values=c("Up"="#D6604D","Down"="#2166AC","NS"="grey70")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="black", linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", linewidth=0.4) +
  geom_text_repel(data=top_label, aes(label=label),
                  size=2.8, color="black",
                  box.padding=0.3, max.overlaps=20,
                  show.legend=FALSE) +
  xlab("log2 Fold Change (shKDM5B / shSCR)") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Volcano Plot — shKDM5B vs shSCR",
          subtitle=paste0("MCF-7 KDM5B knockdown  |  Up: ",
                          sum(vol_df$category=="Up"), "  Down: ",
                          sum(vol_df$category=="Down"),
                          "  (padj<0.05, |LFC|>1)")) +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.title=element_blank(),
        panel.grid.minor=element_blank())

ggsave(file.path(OUTPUT_DIR,"volcano_plot.pdf"), vol_plot, width=9, height=7)
ggsave(file.path(OUTPUT_DIR,"volcano_plot.png"), vol_plot, width=9, height=7, dpi=150)

# -----------------------------------------------------------------------------
# 10. Heatmap — Top 50 DE genes with gene symbols
# -----------------------------------------------------------------------------
cat("Generating heatmap...\n")
top50_ids <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  head(50) %>%
  pull(gene_id)

if (length(top50_ids) >= 2) {
  vst_mat    <- assay(vsd)[top50_ids, ]
  vst_scaled <- t(scale(t(vst_mat)))

  # Replace ENSEMBL IDs with gene symbols on row labels
  symbol_map <- res_df %>%
    filter(gene_id %in% top50_ids) %>%
    dplyr::select(gene_id, gene_symbol) %>%
    mutate(label = ifelse(is.na(gene_symbol), gene_id, gene_symbol))
  rownames(vst_scaled) <- symbol_map$label[match(rownames(vst_scaled), symbol_map$gene_id)]

  col_ann <- data.frame(Condition = sample_info$condition,
                        row.names = rownames(sample_info))
  ann_colors <- list(Condition = c("shSCR"="#2166AC","shKDM5B"="#D6604D"))

  pdf(file.path(OUTPUT_DIR,"heatmap_top50.pdf"), width=10, height=14)
  pheatmap(vst_scaled,
           annotation_col    = col_ann,
           annotation_colors = ann_colors,
           color             = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           cluster_rows=TRUE, cluster_cols=TRUE,
           show_rownames=TRUE, show_colnames=TRUE,
           fontsize_row=8, fontsize_col=10,
           main="Top 50 DE Genes — shKDM5B vs shSCR\n(Z-score of VST counts | MCF-7)")
  dev.off()
  png(file.path(OUTPUT_DIR,"heatmap_top50.png"), width=1000, height=1400, res=120)
  pheatmap(vst_scaled,
           annotation_col    = col_ann,
           annotation_colors = ann_colors,
           color             = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           cluster_rows=TRUE, cluster_cols=TRUE,
           show_rownames=TRUE, show_colnames=TRUE,
           fontsize_row=8, fontsize_col=10,
           main="Top 50 DE Genes — shKDM5B vs shSCR\n(Z-score of VST counts | MCF-7)")
  dev.off()
}

# =============================================================================
# PATHWAY ANALYSIS
# =============================================================================

# Prepare gene lists
sig_genes <- res_df %>%
  filter(padj < 0.05, !is.na(entrez_id))

up_entrez   <- sig_genes %>% filter(log2FoldChange >  1) %>% pull(entrez_id)
down_entrez <- sig_genes %>% filter(log2FoldChange < -1) %>% pull(entrez_id)
all_entrez  <- sig_genes %>% pull(entrez_id)

# Background = all tested genes
bg_entrez <- res_df %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id)

cat(sprintf("\nPathway analysis input:\n"))
cat(sprintf("  Upregulated genes   (|LFC|>1): %d\n", length(up_entrez)))
cat(sprintf("  Downregulated genes (|LFC|>1): %d\n", length(down_entrez)))
cat(sprintf("  All DE genes:                  %d\n", length(all_entrez)))
cat(sprintf("  Background genes:              %d\n\n", length(bg_entrez)))

# -----------------------------------------------------------------------------
# 11. ORA — GO Biological Process
# -----------------------------------------------------------------------------
cat("Running GO Biological Process ORA...\n")
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
            file.path(OUTPUT_DIR,"pathway","go_bp_results.csv"),
            row.names=FALSE)

  p <- dotplot(go_bp, showCategory=20, title="GO Biological Process — All DE Genes",
               font.size=10) +
    theme(plot.title=element_text(face="bold", size=12))
  ggsave(file.path(OUTPUT_DIR,"pathway","go_bp_dotplot.pdf"), p, width=10, height=10)
  ggsave(file.path(OUTPUT_DIR,"pathway","go_bp_dotplot.png"), p, width=10, height=10, dpi=150)

  # Separate up vs down
  go_up <- enrichGO(gene=up_entrez, universe=bg_entrez, OrgDb=org.Hs.eg.db,
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
  go_dn <- enrichGO(gene=down_entrez, universe=bg_entrez, OrgDb=org.Hs.eg.db,
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)

  if (!is.null(go_up) && nrow(go_up) > 0) {
    write.csv(as.data.frame(go_up),
              file.path(OUTPUT_DIR,"pathway","go_bp_upregulated.csv"), row.names=FALSE)
    p2 <- dotplot(go_up, showCategory=15, title="GO BP — Upregulated in shKDM5B") +
      theme(plot.title=element_text(face="bold"))
    ggsave(file.path(OUTPUT_DIR,"pathway","go_bp_up_dotplot.pdf"), p2, width=10, height=9)
    ggsave(file.path(OUTPUT_DIR,"pathway","go_bp_up_dotplot.png"), p2, width=10, height=9, dpi=150)
  }

  if (!is.null(go_dn) && nrow(go_dn) > 0) {
    write.csv(as.data.frame(go_dn),
              file.path(OUTPUT_DIR,"pathway","go_bp_downregulated.csv"), row.names=FALSE)
    p3 <- dotplot(go_dn, showCategory=15, title="GO BP — Downregulated in shKDM5B") +
      theme(plot.title=element_text(face="bold"))
    ggsave(file.path(OUTPUT_DIR,"pathway","go_bp_dn_dotplot.pdf"), p3, width=10, height=9)
    ggsave(file.path(OUTPUT_DIR,"pathway","go_bp_dn_dotplot.png"), p3, width=10, height=9, dpi=150)
  }
  cat("GO BP analysis complete.\n")
} else {
  cat("WARNING: No significant GO BP terms found.\n")
}

# -----------------------------------------------------------------------------
# 12. ORA — KEGG Pathways
# -----------------------------------------------------------------------------
cat("Running KEGG ORA...\n")
kegg_res <- enrichKEGG(
  gene          = all_entrez,
  universe      = bg_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  kegg_res <- setReadable(kegg_res, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
  write.csv(as.data.frame(kegg_res),
            file.path(OUTPUT_DIR,"pathway","kegg_results.csv"), row.names=FALSE)

  p <- dotplot(kegg_res, showCategory=20, title="KEGG Pathways — All DE Genes",
               font.size=10) +
    theme(plot.title=element_text(face="bold", size=12))
  ggsave(file.path(OUTPUT_DIR,"pathway","kegg_dotplot.pdf"), p, width=10, height=9)
  ggsave(file.path(OUTPUT_DIR,"pathway","kegg_dotplot.png"), p, width=10, height=9, dpi=150)
  cat("KEGG analysis complete.\n")
} else {
  cat("WARNING: No significant KEGG pathways found.\n")
}

# -----------------------------------------------------------------------------
# 13. ORA — Reactome Pathways
# -----------------------------------------------------------------------------
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
            file.path(OUTPUT_DIR,"pathway","reactome_results.csv"), row.names=FALSE)

  p <- dotplot(reactome_res, showCategory=20, title="Reactome Pathways — All DE Genes",
               font.size=10) +
    theme(plot.title=element_text(face="bold", size=12))
  ggsave(file.path(OUTPUT_DIR,"pathway","reactome_dotplot.pdf"), p, width=12, height=10)
  ggsave(file.path(OUTPUT_DIR,"pathway","reactome_dotplot.png"), p, width=12, height=10, dpi=150)
  cat("Reactome analysis complete.\n")
} else {
  cat("WARNING: No significant Reactome pathways found.\n")
}

# -----------------------------------------------------------------------------
# 14. GSEA — MSigDB Hallmark Gene Sets
# -----------------------------------------------------------------------------
cat("Running GSEA with MSigDB Hallmarks...\n")

# Build ranked gene list: signed -log10(padj) × sign(LFC)
ranked_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(entrez_id)) %>%
  mutate(rank_score = sign(log2FoldChange) * -log10(padj + 1e-300)) %>%
  arrange(desc(rank_score)) %>%
  filter(!duplicated(entrez_id))

gene_list        <- ranked_df$rank_score
names(gene_list) <- ranked_df$entrez_id

# Hallmark gene sets for Homo sapiens
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

gsea_res <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = hallmarks,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize    = 15,
  maxGSSize    = 500,
  eps          = 1e-10,
  seed         = 42,
  verbose      = FALSE
)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  write.csv(as.data.frame(gsea_res),
            file.path(OUTPUT_DIR,"pathway","gsea_hallmarks_results.csv"),
            row.names=FALSE)

  # Dotplot
  p <- dotplot(gsea_res, showCategory=20, split=".sign",
               title="GSEA — MSigDB Hallmarks") +
    facet_grid(.~.sign) +
    theme(plot.title=element_text(face="bold"),
          axis.text.y=element_text(size=9))
  ggsave(file.path(OUTPUT_DIR,"pathway","gsea_hallmarks_dotplot.pdf"), p, width=14, height=10)
  ggsave(file.path(OUTPUT_DIR,"pathway","gsea_hallmarks_dotplot.png"), p, width=14, height=10, dpi=150)

  # Top enrichment plots for top 3 activated and top 3 suppressed
  gsea_df   <- as.data.frame(gsea_res) %>% arrange(NES)
  top_act   <- tail(gsea_df$ID, 3)
  top_supp  <- head(gsea_df$ID, 3)
  top_paths <- c(top_supp, top_act)

  pdf(file.path(OUTPUT_DIR,"pathway","gsea_top_enrichment_plots.pdf"), width=10, height=6)
  for (pathway in top_paths) {
    tryCatch({
      p <- gseaplot2(gsea_res, geneSetID=pathway,
                     title=gsub("HALLMARK_","",pathway))
      print(p)
    }, error=function(e) cat("Skipping plot for", pathway, "\n"))
  }
  dev.off()

  cat("GSEA analysis complete.\n")

  # Print top results
  cat("\nTop activated hallmarks (positive NES):\n")
  print(as.data.frame(gsea_res) %>%
          filter(NES > 0) %>%
          arrange(p.adjust) %>%
          head(5) %>%
          dplyr::select(ID, NES, p.adjust))

  cat("\nTop suppressed hallmarks (negative NES):\n")
  print(as.data.frame(gsea_res) %>%
          filter(NES < 0) %>%
          arrange(p.adjust) %>%
          head(5) %>%
          dplyr::select(ID, NES, p.adjust))
} else {
  cat("WARNING: No significant GSEA hallmarks found at padj<0.05.\n")
}

# -----------------------------------------------------------------------------
# 15. Final Summary
# -----------------------------------------------------------------------------
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Main results:\n")
cat("  deseq2_results_annotated.csv\n\n")
cat("Plots:\n")
cat("  pca_plot.pdf/png\n")
cat("  ma_plot.pdf/png\n")
cat("  volcano_plot.pdf/png\n")
cat("  heatmap_top50.pdf/png\n\n")
cat("Pathway analysis (results/deseq2/pathway/):\n")
cat("  go_bp_results.csv + dotplots (all, up, down)\n")
cat("  kegg_results.csv + dotplot\n")
cat("  reactome_results.csv + dotplot\n")
cat("  gsea_hallmarks_results.csv + dotplot + enrichment plots\n\n")
cat(sprintf("DE summary: %d up | %d down | %d total (padj<0.05)\n",
            sig_up, sig_down, sig_up+sig_down))
