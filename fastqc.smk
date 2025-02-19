import os
from os.path import join

# Dynamically fetch sample IDs from the given directory.
IDS, = glob_wildcards(join(config["fastq_files_path"], "{sample}_R1.fastq.gz"))

# Rule to perform FastQC on raw FASTQ files.
rule fastqc_raw:
    input:
        r1 = lambda wc: join(config["fastq_files_path"], f"{wc.sample}_R1.fastq.gz"),
        r2 = lambda wc: join(config["fastq_files_path"], f"{wc.sample}_R2.fastq.gz")
    output:
        html_r1 = "{sample}_R1_fastqc.html",
        zip_r1  = "{sample}_R1_fastqc.zip",
        html_r2 = "{sample}_R2_fastqc.html",
        zip_r2  = "{sample}_R2_fastqc.zip"
    threads: 8
    shell:
        """
        {config[fastqc_bin]} --quiet -t {threads} --outdir . {input.r1}
        {config[fastqc_bin]} --quiet -t {threads} --outdir . {input.r2}
        """

# Rule to perform FastQC on trimmed FASTQ files.
rule fastqc_trimmed:
    input:
        r1 = "trimmed/{sample}_trimmed_R1.fastq.gz",
        r2 = "trimmed/{sample}_trimmed_R2.fastq.gz"
    output:
        html_r1 = "trimmed/{sample}_trimmed_R1_fastqc.html",
        zip_r1  = "trimmed/{sample}_trimmed_R1_fastqc.zip",
        html_r2 = "trimmed/{sample}_trimmed_R2_fastqc.html",
        zip_r2  = "trimmed/{sample}_trimmed_R2_fastqc.zip"
    threads: 8
    shell:
        """
        {config[fastqc_bin]} --quiet -t {threads} --outdir . {input.r1}
        {config[fastqc_bin]} --quiet -t {threads} --outdir . {input.r2}
        """

# The 'all' rule specifies the final targets for both raw and trimmed analyses.
rule all:
    input:
        # Raw files outputs
        expand("{sample}_R1_fastqc.html", sample=IDS),
        expand("{sample}_R2_fastqc.html", sample=IDS),
        expand("{sample}_R1_fastqc.zip",  sample=IDS),
        expand("{sample}_R2_fastqc.zip",  sample=IDS),
        # Trimmed files outputs
        expand("trimmed/{sample}_trimmed_R1_fastqc.html", sample=IDS),
        expand("trimmed/{sample}_trimmed_R2_fastqc.html", sample=IDS),
        expand("trimmed/{sample}_trimmed_R1_fastqc.zip",  sample=IDS),
        expand("trimmed/{sample}_trimmed_R2_fastqc.zip",  sample=IDS)
