# NGS Variant Calling Workflow

⭐ **Essential** reproducible Snakemake pipeline for WES/WGS  
**FastQC → Trimming (Trimmomatic) → BWA alignment → GATK HaplotypeCaller → ANNOVAR annotation**

**Includes:** config files + MultiQC run report

**Tags:** Linux / Bash • NGS • Snakemake

## Workflow
1. FastQC (raw reads)
2. Trimming with Trimmomatic (PE, adapter clipping + quality trimming)
3. BWA-MEM alignment + sorting
4. GATK HaplotypeCaller (germline variants → VCF per sample)
5. ANNOVAR annotation (refGene + common databases)
6. MultiQC report (final run report)

## Quick Start
```bash
# 1. Clone
git clone https://github.com/YOURUSER/ngs-variant-calling-snakemake.git
cd ngs-variant-calling-snakemake

# 2. Prepare reference (once)
#    bwa index ref.fa
#    samtools faidx ref.fa
#    gatk CreateSequenceDictionary -R ref.fa -O ref.dict

# 3. Prepare ANNOVAR (once)
#    Download ANNOVAR → unzip → set path in config/config.yaml
#    Download databases (humandb/): perl annotate_variation.pl -buildver hg38 -downdb ...

# 4. Create conda environment
conda env create -f environment.yaml -n ngs-pipeline
conda activate ngs-pipeline

# 5. Edit config/config.yaml (paths + samples)

# 6. Run pipeline
snakemake --cores 16 --keep-going --printshellcmds


Generate HTML run report (optional):
Bashsnakemake --report results/run_report.html
Output directory: results/

fastqc/
trimmed/
aligned/ (sorted BAM + BAI)
variants/ (*.vcf)
annotated/ (*.hg38_multianno.txt)
multiqc_report.html (run report)

Requirements

Conda / Mamba
Reference genome (hg38/hg19) + indexes
ANNOVAR + humandb (configured in config.yaml)

Configuration
All parameters are in config/config.yaml (samples, paths, threads, ANNOVAR protocol, etc.).


