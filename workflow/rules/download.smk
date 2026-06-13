rule download_rnaseq:
    """
    Downloads paired-end RNA-seq reads from NCBI SRA.
    Accession is looked up via the sample sheet; files are renamed
    to match the sample wildcard after download.

    Input:  Sample name wildcard -> accession from samples.tsv
    Output: Paired FASTQs at data/rnaseq/{sample}_1/2.fastq
    """
    output:
        left  = "data/rnaseq/{sample}_1.fastq",
        right = "data/rnaseq/{sample}_2.fastq"
    params:
        accession = lambda wc: samples.loc[wc.sample, "accession"]
    conda:
        "workflow/envs/main.yaml"
    threads: 4
    shell:
        """
        fasterq-dump {params.accession} -O data/rnaseq/ \
            --threads {threads} --split-files
        mv data/rnaseq/{params.accession}_1.fastq {output.left}
        mv data/rnaseq/{params.accession}_2.fastq {output.right}
        """


rule download_cdna:
    """
    Downloads the S. cerevisiae cDNA FASTA from Ensembl Fungi.
    Used as input to salmon index (transcript-level quantification).
    Input:  URL from config["cdna_url"]
    Output: Decompressed FASTA at data/reference/cdna.fa
    """
    output:
        "data/reference/cdna.fa"
    params:
        url = config["cdna_url"]
    conda:
        "workflow/envs/main.yaml"
    shell:
        "wget -O {output}.gz {params.url} && gunzip {output}.gz"