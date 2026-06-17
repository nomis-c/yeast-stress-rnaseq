rule salmon_index:
    input:
        cdna=rules.download_cdna.output[0],
    output:
        directory("results/salmon_index"),
    conda:
        "../envs/env.yaml"
    threads: 8
    shell:
        "salmon index -t {input.cdna} -i {output} -p {threads}"


rule salmon_quant:
    input:
        left=rules.fastp.output.left,
        right=rules.fastp.output.right,
        index=rules.salmon_index.output[0],
    output:
        quant="results/salmon/{sample}/quant.sf",
    params:
        outdir="results/salmon/{sample}",
    conda:
        "../envs/env.yaml"
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
