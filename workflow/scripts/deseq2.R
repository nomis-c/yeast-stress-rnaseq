library(DESeq2)
library(tximport)
library(pheatmap)

quant_files <- snakemake@input[["quant_files"]]
samples_file <- snakemake@input[["samples"]]
outdir <- snakemake@params[["outdir"]]
sample_ids <- snakemake@params[["sample_ids"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load sample metadata
samples <- read.table(samples_file, header = TRUE, comment.char = "#")
samples <- samples[, c("sample", "condition")]
colnames(samples) <- c("Sample", "Condition")
rownames(samples) <- samples$Sample

names(quant_files) <- sample_ids

# Build tx2gene from first quant.sf
# Ensembl fungi cDNA transcript IDs: YAL001C_mRNA -> gene YAL001C
first_quant <- read.table(quant_files[1], header = TRUE, sep = "\t")
tx_names <- first_quant$Name
tx2gene <- data.frame(
    tx_id = tx_names,
    gene_id = sub("_[^_]+$", "", tx_names),
    stringsAsFactors = FALSE
)
unmapped <- tx2gene$gene_id == tx2gene$tx_id
if (any(unmapped)) {
    warning(sprintf(
        "%d transcript(s) have no underscore suffix and will use the full transcript name as gene ID: %s",
        sum(unmapped),
        paste(head(tx2gene$tx_id[unmapped], 5), collapse = ", ")
    ))
}

# Import Salmon quantification
txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE
)

# Guard: all sample IDs from quant files must appear in the sample sheet
stopifnot(
    "Sample IDs in quant files do not match sample sheet" =
        all(colnames(txi$counts) %in% rownames(samples))
)

# Build DESeqDataSet
samples_ordered <- samples[colnames(txi$counts), ]

dsdata <- DESeqDataSetFromTximport(
    txi = txi,
    colData = samples_ordered,
    design = ~ Condition
)

# Filter low-count genes: >=10 reads in >=2 samples
keep <- rowSums(counts(dsdata) >= 10) >= 2
dsdata <- dsdata[keep, ]
message(sprintf("Genes retained after low-count filtering: %d", sum(keep)))

# Fit model once here so downstream contrast/padj jobs only call results()
dds <- DESeq(dsdata)

# Transformations
rld <- rlog(dds, blind = FALSE)
vst <- vst(dds, blind = FALSE)

# Sample-distance heatmap (VST)
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
png(file.path(outdir, "sample_distance_heatmap.png"), width = 800, height = 600)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    annotation_col = samples_ordered["Condition"],
    main = "Sample-to-sample distances (VST)"
)
dev.off()

# PCA plots
png(file.path(outdir, "pca_normalized.png"), width = 800, height = 600)
print(plotPCA(rld, intgroup = c("Sample", "Condition")))
dev.off()

png(file.path(outdir, "pca_unnormalized.png"), width = 800, height = 600)
print(plotPCA(DESeqTransform(dds), intgroup = c("Sample", "Condition")))
dev.off()

# Export
rldm <- assay(rld)
write.table(as.data.frame(rldm),
            file = snakemake@output[["rldm"]], sep = "\t", quote = FALSE)

saveRDS(dds, file = snakemake@output[["dds"]])
message("deseq2.R complete.")
