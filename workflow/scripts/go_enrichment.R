library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)

de_files  <- snakemake@input[["de_results"]]
contrasts <- snakemake@params[["contrasts"]]
outdir    <- snakemake@params[["outdir"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Named list: label -> path to DE TSV
de_paths <- setNames(
    de_files,
    sapply(contrasts, function(c) gsub(",", "_vs_", c))
)

sig_genes_per_contrast <- list()


# Per-contrast GO enrichment

for (contrast in contrasts) {
    label       <- gsub(",", "_vs_", contrast)
    de_tsv_path <- de_paths[[label]]
    go_tsv_path <- file.path(outdir, paste0(label, "_GO.tsv"))

    message("GO enrichment for: ", label)

    de_df     <- read.table(de_tsv_path, header = TRUE, sep = "\t", row.names = 1)
    gene_list <- rownames(de_df)

    sig_genes_per_contrast[[label]] <- gene_list

    if (length(gene_list) == 0) {
        message("  No significant genes - skipping.")
        write.table(data.frame(), file = go_tsv_path, sep = "\t", quote = FALSE)
        next
    }

    tryCatch({
        ego <- enrichGO(
            gene          = gene_list,
            OrgDb         = org.Sc.sgd.db,
            keyType       = "ORF",       # e.g. YAL001C
            ont           = "BP",        # Biological Process
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.05,
            qvalueCutoff  = 0.2,
            readable      = TRUE
        )

        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
            write.table(as.data.frame(ego),
                        file = go_tsv_path, sep = "\t", quote = FALSE)

            png(file.path(outdir, paste0(label, "_GO_dotplot.png")),
                width = 1000, height = 800)
            print(dotplot(ego, showCategory = 20,
                          title = paste(label, "- GO Enrichment (BP)")))
            dev.off()

            message(sprintf("  Enriched GO terms: %d", nrow(as.data.frame(ego))))
        } else {
            message("  No significant GO terms found.")
            write.table(data.frame(), file = go_tsv_path, sep = "\t", quote = FALSE)
        }

    }, error = function(e) {
        message("  GO enrichment failed: ", e$message)
        write.table(data.frame(), file = go_tsv_path, sep = "\t", quote = FALSE)
    })
}


# Overlap: GO enrichment on genes significant across ALL contrasts

if (length(sig_genes_per_contrast) >= 2) {
    shared_genes <- Reduce(intersect, sig_genes_per_contrast)

    if (length(shared_genes) > 0) {
        message(sprintf("GO enrichment on %d shared genes across all contrasts.",
                        length(shared_genes)))

        tryCatch({
            ego_shared <- enrichGO(
                gene          = shared_genes,
                OrgDb         = org.Sc.sgd.db,
                keyType       = "ORF",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE
            )

            if (!is.null(ego_shared) && nrow(as.data.frame(ego_shared)) > 0) {
                write.table(as.data.frame(ego_shared),
                            file = file.path(outdir, "overlap_all_contrasts_GO.tsv"),
                            sep = "\t", quote = FALSE)

                png(file.path(outdir, "overlap_all_contrasts_GO_dotplot.png"),
                    width = 1000, height = 800)
                print(dotplot(ego_shared, showCategory = 20,
                              title = "GO Enrichment - genes significant across all contrasts"))
                dev.off()
            }
        }, error = function(e) {
            message("GO enrichment for shared genes failed: ", e$message)
        })
    } else {
        message("No shared genes across all contrasts.")
    }
}
