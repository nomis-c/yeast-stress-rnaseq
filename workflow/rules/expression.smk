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
        expand(
            "results/deseq2/{contrast}.tsv",
            contrast=[c.replace(",", "_vs_") for c in config["contrasts"]],
        ),
    params:
        contrasts=config["contrasts"],
        outdir="results/deseq2",
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/differential_expression.R"
