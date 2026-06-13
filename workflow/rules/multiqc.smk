rule multiqc:
    input:
        fastp_json  = expand("results/qc/fastp/{sample}.json",
                             sample=samples.index),
        salmon_dirs = expand("results/salmon/{sample}/quant.sf",
                             sample=samples.index),
        fastqc_html = expand("results/qc/fastqc/{accession}_fastqc.html",
                             accession=config["accession"])
    output:
        "results/qc/multiqc_report.html"
    conda:
        "workflow/envs/main.yaml"
    shell:
        "multiqc results/ -o results/qc/ -n multiqc_report --force"