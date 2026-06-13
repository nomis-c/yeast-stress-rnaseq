rule fastp:
    input:
        left=rules.download_rnaseq.output.left,
        right=rules.download_rnaseq.output.right,
    output:
        left="results/rnaseq_trimmed/{sample}_1.fastq",
        right="results/rnaseq_trimmed/{sample}_2.fastq",
        json="results/qc/fastp/{sample}.json",
        html="results/qc/fastp/{sample}.html",
    conda:
        "../envs/env.yaml"
    threads: 4
    shell:
        """
        fastp \
            -i {input.left} -I {input.right} \
            -o {output.left} -O {output.right} \
            -j {output.json} -h {output.html} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --correction \
            --thread {threads}
        """
