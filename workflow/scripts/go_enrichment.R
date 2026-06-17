library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)
library(enrichplot)

de_tsv_path <- snakemake@input[["de_results"]]
go_tsv_path <- snakemake@output[["go_results"]]
outdir <- snakemake@params[["outdir"]]
padj_cutoff <- as.numeric(snakemake@params[["padj_cutoff"]])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

contrast_label <- snakemake@wildcards[["contrast"]]
padj_label <- snakemake@wildcards[["padj"]]
label <- paste0(contrast_label, "_padj", padj_label)

message("GO enrichment for: ", label)

de_df <- read.table(de_tsv_path, header = TRUE, sep = "\t", row.names = 1)
gene_list <- rownames(de_df)

if (length(gene_list) == 0) {
    message("  No significant genes - skipping.")
    write.table(data.frame(), file = go_tsv_path, sep = "\t", quote = FALSE)
} else {
    tryCatch({
        ego <- enrichGO(
            gene = gene_list,
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
                        file = go_tsv_path, sep = "\t", quote = FALSE)

            # barplot
            png(file.path(outdir, paste0(label, "_GO_barplot.png")),
                width = 1000, height = 800)
            print(barplot(ego, showCategory = 20,
                          title = paste(label, "- GO Enrichment (BP)")))
            dev.off()

            # cnetplot
            png(file.path(outdir, paste0(label, "_GO_cnetplot.png")),
                width = 1400, height = 1200)
            print(cnetplot(ego, showCategory = 5))
            dev.off()

            message(sprintf("  Enriched GO terms: %d", nrow(as.data.frame(ego))))
        } else {
            message("  No significant GO terms found.")
            write.table(data.frame(), file = go_tsv_path, sep = "\t", quote = FALSE)
        }

    }, error = function(e) {
        stop("GO enrichment failed: ", e$message)
    })
}
