rule go_enrichment:
    input:
        de_results=rules.differential_expression.output[0],
        rldm=rules.deseq2_pca.output.rldm,
    output:
        go_up="results/enrichment/{contrast}_padj{padj}_GO_up.tsv",
        go_down="results/enrichment/{contrast}_padj{padj}_GO_down.tsv",
    params:
        outdir="results/enrichment",
        padj_cutoff="{padj}",
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/go_enrichment.R"
