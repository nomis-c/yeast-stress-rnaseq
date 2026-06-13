rule go_enrichment:
    """
    Runs GO enrichment analysis on significant DE genes per contrast.

    For each contrast:
      - Reads the significant DE TSV from rule differential_expression
      - Tests for over-represented GO Biological Process terms (BH correction)
      - Exports GO results as TSV
      - Saves a dot-plot (gene ratio vs GO term, colored by p.adjust)

    Also runs enrichment on the intersection of significant genes across
    all contrasts if multiple contrasts are defined.

    Input:  Per-contrast DE TSV files from rule differential_expression
    Output: Per-contrast GO TSV files in results/enrichment/
    """
    input:
        de_results = expand("results/deseq2/{contrast}.tsv",
                            contrast=[c.replace(",", "_vs_") for c in config["contrasts"]])
    output:
        go_results = expand("results/enrichment/{contrast}_GO.tsv",
                            contrast=[c.replace(",", "_vs_") for c in config["contrasts"]])
    params:
        contrasts = config["contrasts"],
        outdir    = "results/enrichment"
    conda:
        "workflow/envs/r.yaml"
    script:
        "../scripts/go_enrichment.R"