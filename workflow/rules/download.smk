rule download_rnaseq:
    output:
        left="resources/data/rnaseq/{sample}_1.fastq",
        right="resources/data/rnaseq/{sample}_2.fastq",
    params:
        accession=lambda wc: samples.loc[wc.sample, "accession"],
    conda:
        "../envs/env.yaml"
    threads: 4
    resources:
        tmpdir="tmp",
    shell:
        """
        fasterq-dump {params.accession} -O resources/data/rnaseq/ \
            --threads {threads} --split-files \
            --temp {resources.tmpdir}
        mv resources/data/rnaseq/{params.accession}_1.fastq {output.left}
        mv resources/data/rnaseq/{params.accession}_2.fastq {output.right}
        """


rule download_cdna:
    output:
        "resources/data/reference/cdna.fa",
    params:
        url=config["cdna_url"],
    conda:
        "../envs/env.yaml"
    shell:
        """
        mkdir -p resources/data/reference
        curl -A "Mozilla/5.0" -L -o resources/data/reference/cdna.fa.gz {params.url}
        gunzip resources/data/reference/cdna.fa.gz
        """
