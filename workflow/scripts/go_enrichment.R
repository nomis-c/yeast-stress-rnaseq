library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)
library(enrichplot)

de_tsv_path <- snakemake@input[["de_results"]]
rldm_path <- snakemake@input[["rldm"]]
go_up_path <- snakemake@output[["go_up"]]
go_down_path <- snakemake@output[["go_down"]]
outdir <- snakemake@params[["outdir"]]
padj_cutoff <- as.numeric(snakemake@params[["padj_cutoff"]])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

contrast_label <- snakemake@wildcards[["contrast"]]
padj_label <- snakemake@wildcards[["padj"]]
label <- paste0(contrast_label, "_padj", padj_label)

message("GO enrichment for: ", label)

de_df <- read.table(de_tsv_path, header = TRUE, sep = "\t", row.names = 1)

# Background universe for ORA: the genes that survived the low-count filter in
# deseq2.R (rldm holds exactly these retained genes), i.e. the genes actually
# tested for DE. Without this, enrichGO defaults to the whole annotated genome,
# which inflates enrichment significance.
rldm <- read.table(rldm_path, header = TRUE, sep = "\t", row.names = 1)
universe <- rownames(rldm)

# Run enrichGO on one directional gene set and write its table + plots.
# Up- and down-regulated genes are kept separate so the resulting BP terms stay
# direction-interpretable (mixing them makes the larger set dominate).
run_go <- function(genes, direction, out_path) {
    if (length(genes) == 0) {
        message("  No ", direction, "-regulated genes - skipping.")
        write.table(data.frame(), file = out_path, sep = "\t", quote = FALSE)
        return(invisible(NULL))
    }

    ego <- enrichGO(
        gene = genes,
        universe = universe,
        OrgDb = org.Sc.sgd.db,
        keyType = "ORF",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = padj_cutoff,
        qvalueCutoff = 0.2,
        readable = FALSE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.table(as.data.frame(ego),
                    file = out_path, sep = "\t", quote = FALSE)

        # barplot
        png(file.path(outdir, paste0(label, "_GO_", direction, "_barplot.png")),
            width = 1000, height = 800)
        print(barplot(ego, showCategory = 20,
                      title = paste(label, "-", direction,
                                    "GO Enrichment (BP)")))
        dev.off()

        # cnetplot
        png(file.path(outdir, paste0(label, "_GO_", direction, "_cnetplot.png")),
            width = 1400, height = 1200)
        print(cnetplot(ego, showCategory = 5))
        dev.off()

        message(sprintf("  %s-regulated enriched GO terms: %d",
                        direction, nrow(as.data.frame(ego))))
    } else {
        message("  No significant ", direction, "-regulated GO terms found.")
        write.table(data.frame(), file = out_path, sep = "\t", quote = FALSE)
    }
}

# Split the significant DE set by direction of regulation.
up_genes <- rownames(de_df[!is.na(de_df$log2FoldChange) & de_df$log2FoldChange > 0, , drop = FALSE])
down_genes <- rownames(de_df[!is.na(de_df$log2FoldChange) & de_df$log2FoldChange < 0, , drop = FALSE])

tryCatch({
    run_go(up_genes, "up", go_up_path)
    run_go(down_genes, "down", go_down_path)
}, error = function(e) {
    stop("GO enrichment failed: ", e$message)
})
