rule multiqc:
    input:
        fastp_json  = expand("results/qc/fastp/{sample}.json",
                             sample=samples.index),
        salmon_dirs = expand("results/salmon/{sample}/quant.sf",
                             sample=samples.index),
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc results/ -o results/qc/ -n multiqc_report --force"
