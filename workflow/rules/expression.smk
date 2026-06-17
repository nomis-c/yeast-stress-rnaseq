rule deseq2_pca:
    input:
        quant_files=expand("results/salmon/{sample}/quant.sf", sample=samples.index),
        samples=config["samples"],
    output:
        rldm="results/deseq2/rldm.tsv",
        dsdata="results/deseq2/dsdata.rds",
    params:
        outdir="results/deseq2/plots",
        sample_ids=list(samples.index),
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/deseq2.R"


rule differential_expression:
    input:
        dsdata=rules.deseq2_pca.output.dsdata,
        rldm=rules.deseq2_pca.output.rldm,
    output:
        "results/deseq2/{contrast}_padj{padj}.tsv",
    params:
        outdir="results/deseq2",
        lfc_cutoff=config["lfc_cutoff"],
        padj_cutoff="{padj}",
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/differential_expression.R"
