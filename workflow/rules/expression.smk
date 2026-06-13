rule deseq2_pca:
    """
    Loads Salmon quant.sf files via tximport, builds a DESeqDataSet,
    applies rlog and VST normalization, and generates QC plots.

    Input:  quant.sf files from rule salmon_quant + samples.tsv
    Output: rlog matrix (TSV) + serialized DESeqDataSet (RDS)
    """
    input:
        quant_files = expand("results/salmon/{sample}/quant.sf",
                             sample=samples.index),
        samples     = config["samples"]
    output:
        rldm   = "results/deseq2/rldm.tsv",
        dsdata = "results/deseq2/dsdata.rds"
    params:
        outdir     = "results/deseq2/plots",
        sample_ids = list(samples.index)
    conda:
        "workflow/envs/r.yaml"
    script:
        "../scripts/deseq2.R"


rule differential_expression:
    """
    Fits the DESeq2 model and tests for differential expression per contrast.
    For each contrast in config["contrasts"]:
      - Extracts DE results (padj <= 0.05, log2FC > 0)
      - Exports significant results as TSV
      - Generates MA plot, heatmap of significant genes,
        top-10 by p-value, top-10 by fold change
      - Produces overlap heatmap for genes significant across all contrasts

    Input:  DESeqDataSet RDS + rlog matrix from rule deseq2_pca
    Output: Per-contrast TSV files in results/deseq2/
    """
    input:
        dsdata = "results/deseq2/dsdata.rds",
        rldm   = "results/deseq2/rldm.tsv"
    output:
        expand("results/deseq2/{contrast}.tsv",
               contrast=[c.replace(",", "_vs_") for c in config["contrasts"]])
    params:
        contrasts = config["contrasts"],
        outdir    = "results/deseq2"
    conda:
        "workflow/envs/r.yaml"
    script:
        "../scripts/differential_expression.R"