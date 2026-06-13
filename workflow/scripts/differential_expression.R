library(DESeq2)
library(pheatmap)

dsdata_file <- snakemake@input[["dsdata"]]
rldm_file   <- snakemake@input[["rldm"]]
contrasts   <- snakemake@params[["contrasts"]]
outdir      <- snakemake@params[["outdir"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dsdata <- readRDS(dsdata_file)
rldm   <- read.table(rldm_file, header = TRUE, row.names = 1)

# Fit DESeq2 model once
dds <- DESeq(dsdata)

sig_results <- list()

for (contrast in contrasts) {
    parts  <- strsplit(contrast, ",")[[1]]
    cond_a <- parts[1]
    cond_b <- parts[2]
    label  <- paste0(cond_a, "_vs_", cond_b)

    message("Processing contrast: ", label)

    res    <- results(dds, contrast = c("Condition", cond_a, cond_b))
    resSig <- subset(res, padj <= 0.05 & log2FoldChange > 0)

    sig_results[[label]] <- resSig

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
    n_top  <- min(10, nrow(res_df[!is.na(res_df$pvalue), ]))

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
}

# Overlap heatmap across all contrasts
if (length(sig_results) >= 2) {
    shared_genes <- Reduce(intersect, lapply(sig_results, rownames))

    if (length(shared_genes) > 0) {
        sigm_shared <- rldm[shared_genes, , drop = FALSE]
        png(file.path(outdir, "overlap_all_contrasts_heatmap.png"),
            width = 1200, height = 1000)
        pheatmap(sigm_shared, scale = "row", show_rownames = TRUE,
                 key = TRUE, main = "Genes upregulated across all contrasts")
        dev.off()

        write.table(as.data.frame(sigm_shared),
                    file = file.path(outdir, "overlap_all_contrasts.tsv"),
                    sep = "\t", quote = FALSE)

        message(sprintf("Shared genes across all contrasts: %d", length(shared_genes)))
    } else {
        message("No shared genes across all contrasts.")
    }
}

