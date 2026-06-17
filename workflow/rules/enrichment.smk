rule go_enrichment:
    input:
        de_results=rules.differential_expression.output[0],
    output:
        go_results="results/enrichment/{contrast}_padj{padj}_GO.tsv",
    params:
        outdir="results/enrichment",
        padj_cutoff="{padj}",
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/go_enrichment.R"
