import os

configfile: "config/config.yaml"

# ====================== CONFIG ======================
SAMPLES = config["samples"]
OUTDIR = config["outdir"]
THREADS = config["threads"]
REF = config["ref"]["fasta"]
ADAPTERS = config["adapter_file"]

# ====================== TARGETS ======================
rule all:
    input:
        expand(f"{OUTDIR}/annotated/{{sample}}.{config['annovar']['buildver']}_multianno.txt", sample=SAMPLES),
        f"{OUTDIR}/multiqc_report.html"

# ====================== RULES ======================
rule fastqc:
    input:
        r1 = lambda w: f"{config['fastq']['dir']}/{w.sample}{config['fastq']['suffix_R1']}",
        r2 = lambda w: f"{config['fastq']['dir']}/{w.sample}{config['fastq']['suffix_R2']}",
    output:
        html1 = f"{OUTDIR}/fastqc/{{sample}}_R1_fastqc.html",
        zip1  = f"{OUTDIR}/fastqc/{{sample}}_R1_fastqc.zip",
        html2 = f"{OUTDIR}/fastqc/{{sample}}_R2_fastqc.html",
        zip2  = f"{OUTDIR}/fastqc/{{sample}}_R2_fastqc.zip",
    log:
        f"{OUTDIR}/logs/fastqc/{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/fastqc
        fastqc {input.r1} {input.r2} -o {OUTDIR}/fastqc -t {threads} > {log} 2>&1
        """

rule trim:
    input:
        r1 = lambda w: f"{config['fastq']['dir']}/{w.sample}{config['fastq']['suffix_R1']}",
        r2 = lambda w: f"{config['fastq']['dir']}/{w.sample}{config['fastq']['suffix_R2']}",
        adapters = ADAPTERS
    output:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{OUTDIR}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        r1u = f"{OUTDIR}/trimmed/{{sample}}_R1.unpaired.fastq.gz",
        r2u = f"{OUTDIR}/trimmed/{{sample}}_R2.unpaired.fastq.gz",
    log:
        f"{OUTDIR}/logs/trimmomatic/{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/trimmed {OUTDIR}/logs/trimmomatic
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.r1u} \
            {output.r2} {output.r2u} \
            ILLUMINACLIP:{input.adapters}:2:30:10:2:True \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            2> {log}
        """

rule bwa_align:
    input:
        r1 = rules.trim.output.r1,
        r2 = rules.trim.output.r2,
        ref = REF
    output:
        bam = f"{OUTDIR}/aligned/{{sample}}.sorted.bam"
    log:
        f"{OUTDIR}/logs/bwa/{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/aligned {OUTDIR}/logs/bwa
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} | \
        samtools sort -@ {threads} -o {output.bam} -
        """

rule samtools_index:
    input:
        bam = rules.bwa_align.output.bam
    output:
        bai = f"{OUTDIR}/aligned/{{sample}}.sorted.bam.bai"
    shell:
        "samtools index {input.bam}"

rule haplotype_caller:
    input:
        bam = rules.bwa_align.output.bam,
        bai = rules.samtools_index.output.bai,
        ref = REF
    output:
        vcf = f"{OUTDIR}/variants/{{sample}}.vcf"
    log:
        f"{OUTDIR}/logs/gatk/{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/variants {OUTDIR}/logs/gatk
        gatk HaplotypeCaller \
            --reference {input.ref} \
            --input {input.bam} \
            --output {output.vcf} \
            --native-pair-hmm-threads {threads} \
            2> {log}
        """

rule annovar:
    input:
        vcf = rules.haplotype_caller.output.vcf
    output:
        anno = f"{OUTDIR}/annotated/{{sample}}.{config['annovar']['buildver']}_multianno.txt"
    params:
        table_annovar = config["annovar"]["table_annovar"],
        humandb = config["annovar"]["humandb"],
        buildver = config["annovar"]["buildver"],
        protocol = config["annovar"]["protocol"],
        operation = config["annovar"]["operation"]
    log:
        f"{OUTDIR}/logs/annovar/{{sample}}.log"
    shell:
        """
        mkdir -p {OUTDIR}/annotated {OUTDIR}/logs/annovar
        perl {params.table_annovar} {input.vcf} {params.humandb} \
            -buildver {params.buildver} \
            -out {OUTDIR}/annotated/{{wildcards.sample}} \
            -remove -protocol {params.protocol} \
            -operation {params.operation} \
            -nastring . -vcfinput \
            > {log} 2>&1
        """

rule multiqc:
    input:
        expand(f"{OUTDIR}/fastqc/{{sample}}_R1_fastqc.zip", sample=SAMPLES)
    output:
        f"{OUTDIR}/multiqc_report.html"
    log:
        f"{OUTDIR}/logs/multiqc.log"
    shell:
        """
        mkdir -p {OUTDIR}
        multiqc {OUTDIR}/fastqc -o {OUTDIR} --filename multiqc_report.html > {log} 2>&1
        """
