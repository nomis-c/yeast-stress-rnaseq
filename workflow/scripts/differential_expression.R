library(DESeq2)
library(pheatmap)

dsdata_file <- snakemake@input[["dsdata"]]
rldm_file <- snakemake@input[["rldm"]]
outdir <- snakemake@params[["outdir"]]
lfc_cutoff <- snakemake@params[["lfc_cutoff"]]
padj_cutoff <- as.numeric(snakemake@params[["padj_cutoff"]])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dsdata <- readRDS(dsdata_file)
rldm <- read.table(rldm_file, header = TRUE, row.names = 1)

# Fit DESeq2 model once
dds <- DESeq(dsdata)

contrast_label <- snakemake@wildcards[["contrast"]]
padj_label <- snakemake@wildcards[["padj"]]

parts <- strsplit(contrast_label, "_vs_")[[1]]
cond_a <- parts[1]
cond_b <- parts[2]
label <- paste0(cond_a, "_vs_", cond_b, "_padj", padj_label)

message("Processing contrast: ", label)

res <- results(dds, contrast = c("Condition", cond_a, cond_b))
resSig <- subset(res, padj <= padj_cutoff & log2FoldChange >= lfc_cutoff)

write.table(as.data.frame(resSig),
            file = file.path(outdir, paste0(label, ".tsv")),
            sep = "\t", quote = FALSE)

message(sprintf("  Significant genes: %d", nrow(resSig)))

# MA plot
png(file.path(outdir, paste0(label, "_MA.png")), width = 800, height = 600)
plotMA(res, ylim = c(-4, 4), main = label)
dev.off()

# Heatmap of all significant genes
if (nrow(resSig) > 1) {
    sigm <- rldm[rownames(resSig), , drop = FALSE]
    png(file.path(outdir, paste0(label, "_heatmap_sig.png")),
        width = 1200, height = 1000)
    pheatmap(sigm, scale = "row", show_rownames = FALSE, key = TRUE, main = label)
    dev.off()
}

# Top-10 by p-value
res_df <- as.data.frame(res)
n_top <- min(10, nrow(res_df[!is.na(res_df$pvalue), ]))

pval_best <- rownames(res_df[order(res_df$pvalue, na.last = TRUE), ])[seq_len(n_top)]
png(file.path(outdir, paste0(label, "_top10_pval.png")), width = 800, height = 600)
pheatmap(rldm[pval_best, , drop = FALSE],
         scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         key = TRUE, show_rownames = TRUE,
         main = paste(label, "- top", n_top, "by p-value"))
dev.off()

# Top-10 by fold change
fc_best <- rownames(
    res_df[order(res_df$log2FoldChange, decreasing = TRUE, na.last = TRUE), ]
)[seq_len(n_top)]
png(file.path(outdir, paste0(label, "_top10_fc.png")), width = 800, height = 600)
pheatmap(rldm[fc_best, , drop = FALSE],
         scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         key = TRUE, show_rownames = TRUE,
         main = paste(label, "- top", n_top, "by fold change"))
dev.off()
