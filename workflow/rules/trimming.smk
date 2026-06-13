rule fastp:
    """
    Trims adapters and low-quality bases from paired-end RNA-seq reads.

    Key flags:
      --detect_adapter_for_pe       -- auto-detect adapters
      --qualified_quality_phred 20  -- discard bases with Q<20
      --length_required 36          -- discard reads shorter than 36 bp
      --correction                  -- correct mismatched bases in overlaps
    """
    input:
        left  = "data/rnaseq/{sample}_1.fastq",
        right = "data/rnaseq/{sample}_2.fastq"
    output:
        left  = "data/rnaseq_trimmed/{sample}_1.fastq",
        right = "data/rnaseq_trimmed/{sample}_2.fastq",
        json  = "results/qc/fastp/{sample}.json",
        html  = "results/qc/fastp/{sample}.html"
    conda:
        "workflow/envs/main.yaml"
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