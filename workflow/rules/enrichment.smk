rule go_enrichment:
    input:
        de_results = rules.differential_expression.output
    output:
        go_results = expand("results/enrichment/{contrast}_GO.tsv",
                            contrast=[c.replace(",", "_vs_") for c in config["contrasts"]])
    params:
        contrasts = config["contrasts"],
        outdir    = "results/enrichment"
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/go_enrichment.R"
