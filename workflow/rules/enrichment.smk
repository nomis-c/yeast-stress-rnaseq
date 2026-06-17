rule go_enrichment:
    input:
        de_results="results/deseq2/{contrast}_padj{padj}.tsv",
    output:
        go_results="results/enrichment/{contrast}_padj{padj}_GO.tsv",
    params:
        outdir="results/enrichment",
        padj_cutoff="{padj}",
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/go_enrichment.R"
