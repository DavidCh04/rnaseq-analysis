# ==============================================================================
# RNA-seq Differential Expression Analysis with DESeq2
# ==============================================================================
# Project: Sodium Butyrate Effects on Glucolipotoxic Renal Tubular Cells
# Author: David Chaparr Lara
# Date: January 2025
# Description: Analysis of RNA-seq data from HK-2 cells treated with sodium
#              butyrate under glucolipotoxic conditions to identify 
#              differentially expressed genes and enriched pathways.
# ==============================================================================

# Load required libraries ------------------------------------------------------
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# ==============================================================================
# 1. IMPORT SALMON QUANTIFICATIONS
# ==============================================================================

# Load sample metadata
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)
metadata$condition <- factor(metadata$condition, levels = c("GLT", "GLT_NaB"))

# Create file paths for Salmon quantification files
files <- file.path("results/salmon_quants", metadata$sample_id, "quant.sf")
names(files) <- metadata$sample_id

# Verify all files exist
if (!all(file.exists(files))) {
  stop("One or more quantification files not found. Check file paths.")
}

# Load transcript-to-gene mapping
tx2gene <- read.table("tx2gene.txt", header = FALSE, 
                      col.names = c("TXNAME", "GENEID"))

# Import transcript-level estimates with tximport
# txOut = FALSE aggregates to gene-level counts
txi <- tximport(files, 
                type = "salmon", 
                txOut = FALSE, 
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE,
                ignoreAfterBar = TRUE)

# ==============================================================================
# 2. CREATE DESEQ2 OBJECT
# ==============================================================================

# Create DESeqDataSet from tximport object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata,
                                design = ~ condition)

# ==============================================================================
# 3. PRE-FILTERING
# ==============================================================================

# Remove genes with very low expression (< 10 total counts)
# This improves statistical power and reduces multiple testing burden
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes retained after filtering:", nrow(dds), "\n")

# ==============================================================================
# 4. DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

# Run DESeq2 pipeline (normalization, dispersion estimation, statistical testing)
dds <- DESeq(dds)

# Extract results for GLT_NaB vs GLT comparison
res <- results(dds, contrast = c("condition", "GLT_NaB", "GLT"))

# Sort by adjusted p-value
res <- res[order(res$padj), ]

# Print summary
cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
summary(res, alpha = 0.05)

# ==============================================================================
# 5. GENE ANNOTATION
# ==============================================================================

# Add gene symbols to results
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(res),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

res$symbol <- gene_symbols

# ==============================================================================
# 6. SAVE RESULTS
# ==============================================================================

# Create output directory if it doesn't exist
dir.create("results/deseq2", showWarnings = FALSE, recursive = TRUE)

# Save annotated results
res_annotated <- as.data.frame(res)
write.csv(res_annotated, 
          "results/deseq2/differential_expression_results_annotated.csv",
          row.names = TRUE)

cat("\nResults saved to: results/deseq2/differential_expression_results_annotated.csv\n")

# ==============================================================================
# 7. VISUALIZATIONS
# ==============================================================================

# Create figures directory
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- MA Plot ---
pdf("results/figures/MA_plot.pdf", width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: GLT+NaB vs GLT")
dev.off()

# --- Volcano Plot ---
pdf("results/figures/volcano_plot.pdf", width = 10, height = 8)
EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'GLT+NaB vs GLT',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 4.0,
                legendPosition = 'right')
dev.off()

# --- PCA Plot ---
# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)

pdf("results/figures/PCA_plot.pdf", width = 8, height = 6)
plotPCA(vsd, intgroup = "condition") +
  theme_bw() +
  ggtitle("PCA - GLT vs GLT+NaB") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# --- Heatmap of Top 50 Genes ---
# Select top 50 most significant genes
top_gene_names <- rownames(res)[head(order(res$padj), 50)]

# Extract expression matrix
mat <- assay(vsd)[top_gene_names, ]

pdf("results/figures/heatmap_top50.pdf", width = 10, height = 12)
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         scale = "row",
         annotation_col = as.data.frame(colData(dds)[, "condition", drop = FALSE]),
         main = "Top 50 Differentially Expressed Genes")
dev.off()

cat("\nVisualization plots saved to: results/figures/\n")

# ==============================================================================
# 8. PATHWAY ENRICHMENT ANALYSIS
# ==============================================================================

cat("\n=== RUNNING PATHWAY ENRICHMENT ANALYSIS ===\n")

# --- Prepare gene lists ---

# Significant genes (padj < 0.05)
sig_genes <- rownames(res)[!is.na(res$padj) & res$padj < 0.05]

# Convert ENSEMBL to ENTREZ IDs
gene_entrez <- mapIds(org.Hs.eg.db,
                      keys = sig_genes,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
gene_entrez <- gene_entrez[!is.na(gene_entrez)]

# Separate upregulated and downregulated genes
res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

up_genes <- rownames(res_sig)[res_sig$log2FoldChange > 0]
up_entrez <- mapIds(org.Hs.eg.db, 
                    keys = up_genes, 
                    column = "ENTREZID", 
                    keytype = "ENSEMBL", 
                    multiVals = "first")
up_entrez <- up_entrez[!is.na(up_entrez)]

down_genes <- rownames(res_sig)[res_sig$log2FoldChange < 0]
down_entrez <- mapIds(org.Hs.eg.db, 
                      keys = down_genes, 
                      column = "ENTREZID", 
                      keytype = "ENSEMBL", 
                      multiVals = "first")
down_entrez <- down_entrez[!is.na(down_entrez)]

# --- Gene Ontology (GO) Enrichment ---

# GO: Biological Process (all significant genes)
go_bp <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# GO: Cellular Component
go_cc <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# GO: Molecular Function
go_mf <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# GO enrichment for upregulated genes
go_up <- enrichGO(gene = up_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# GO enrichment for downregulated genes
go_down <- enrichGO(gene = down_entrez,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)

# --- KEGG Pathway Enrichment ---
kegg <- enrichKEGG(gene = gene_entrez,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

# Print top results
cat("\n=== TOP 10 GO BIOLOGICAL PROCESS TERMS ===\n")
print(head(go_bp, 10))

cat("\n=== TOP 10 KEGG PATHWAYS ===\n")
print(head(kegg, 10))

# --- Save enrichment results ---

dir.create("results/enrichment", showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(go_bp), "results/enrichment/GO_biological_process.csv")
write.csv(as.data.frame(go_cc), "results/enrichment/GO_cellular_component.csv")
write.csv(as.data.frame(go_mf), "results/enrichment/GO_molecular_function.csv")
write.csv(as.data.frame(kegg), "results/enrichment/KEGG_pathways.csv")
write.csv(as.data.frame(go_up), "results/enrichment/GO_upregulated.csv")
write.csv(as.data.frame(go_down), "results/enrichment/GO_downregulated.csv")

cat("\nEnrichment results saved to: results/enrichment/\n")

# --- Enrichment visualizations ---

# Dot plot - GO Biological Process
pdf("results/figures/enrichment_GO_BP_dotplot.pdf", width = 12, height = 10)
dotplot(go_bp, showCategory = 20) + 
  ggtitle("GO Biological Process Enrichment")
dev.off()

# Bar plot - GO Biological Process
pdf("results/figures/enrichment_GO_BP_barplot.pdf", width = 12, height = 10)
barplot(go_bp, showCategory = 20) + 
  ggtitle("Top 20 Enriched GO Biological Processes")
dev.off()

# KEGG pathways dotplot
pdf("results/figures/enrichment_KEGG_dotplot.pdf", width = 12, height = 8)
dotplot(kegg, showCategory = 20) + 
  ggtitle("KEGG Pathway Enrichment")
dev.off()

# Gene-Concept Network (cnetplot)
if (nrow(go_bp) > 0) {
  pdf("results/figures/enrichment_cnetplot.pdf", width = 14, height = 14)
  cnetplot(go_bp, 
           showCategory = 5,
           circular = FALSE,
           colorEdge = TRUE)
  dev.off()
}

# Heatmap plot
if (nrow(go_bp) > 0) {
  pdf("results/figures/enrichment_heatplot.pdf", width = 12, height = 10)
  heatplot(go_bp, showCategory = 10)
  dev.off()
}

cat("\nEnrichment visualizations saved to: results/figures/\n")

# ==============================================================================
# 9. SEARCH FOR HYPOTHESIS-SPECIFIC TERMS
# ==============================================================================

cat("\n=== SEARCHING FOR AUTOPHAGY/LIPOPHAGY TERMS ===\n")

# Search for autophagy-related terms
autophagy_terms <- grep("autoph|lipoph|lysoso", 
                        go_bp@result$Description, 
                        ignore.case = TRUE, 
                        value = TRUE)

if (length(autophagy_terms) > 0) {
  cat("\nAUTOPHAGY/LIPOPHAGY-RELATED TERMS FOUND:\n")
  autophagy_results <- go_bp@result[go_bp@result$Description %in% autophagy_terms, ]
  print(autophagy_results[, c("Description", "pvalue", "qvalue", "Count")])
} else {
  cat("\nNo direct autophagy terms found in top enriched GO terms\n")
}

# Search for lipid metabolism terms
lipid_terms <- grep("lipid|fatty acid|triglyceride", 
                    go_bp@result$Description, 
                    ignore.case = TRUE, 
                    value = TRUE)

if (length(lipid_terms) > 0) {
  cat("\n\nLIPID METABOLISM TERMS FOUND:\n")
  lipid_results <- go_bp@result[go_bp@result$Description %in% lipid_terms, ]
  print(lipid_results[, c("Description", "pvalue", "qvalue", "Count")])
}

# Search for epithelial barrier terms
barrier_terms <- grep("junction|barrier|adhesion|epithelial", 
                      go_bp@result$Description, 
                      ignore.case = TRUE, 
                      value = TRUE)

if (length(barrier_terms) > 0) {
  cat("\n\nEPITHELIAL BARRIER TERMS FOUND:\n")
  barrier_results <- go_bp@result[go_bp@result$Description %in% barrier_terms, ]
  print(barrier_results[, c("Description", "pvalue", "qvalue", "Count")])
}

# ==============================================================================
# 10. SESSION INFO
# ==============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nSession Information:\n")
sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================