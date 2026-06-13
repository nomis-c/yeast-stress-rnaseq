rule salmon_index:
    """
    Builds a Salmon quasi-mapping index from the cDNA FASTA.
    Input:  cDNA FASTA from rule download_cdna
    Output: Salmon index directory at results/salmon_index/
    """
    input:
        cdna = "data/reference/cdna.fa"
    output:
        directory("results/salmon_index")
    conda:
        "workflow/envs/salmon.yaml"
    threads: 8
    shell:
        "salmon index -t {input.cdna} -i {output} -p {threads}"


rule salmon_quant:
    """
    Quantifies transcript expression from trimmed paired-end reads.

    Key flags:
      -l A               -- auto-detect library type
      --validateMappings -- stricter quasi-mapping
      --gcBias           -- GC content bias correction
      --seqBias          -- sequence-specific bias correction
    """
    input:
        left  = "data/rnaseq_trimmed/{sample}_1.fastq",
        right = "data/rnaseq_trimmed/{sample}_2.fastq",
        index = "results/salmon_index"
    output:
        quant = "results/salmon/{sample}/quant.sf"
    params:
        outdir = "results/salmon/{sample}"
    conda:
        "workflow/envs/salmon.yaml"
    threads: 8
    shell:
        """
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.left} \
            -2 {input.right} \
            -p {threads} \
            -o {params.outdir} \
            --validateMappings \
            --gcBias \
            --seqBias
        """