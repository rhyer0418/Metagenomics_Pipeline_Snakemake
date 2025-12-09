# Automated Metagenomics Pipeline (Snakemake)
**Version:** v1  
**Description:** Auto pipeline for metagenomics binning, refinement, taxonomy, and annotation.

## 1. Project Directory Structure
Before running, organize your working directory as follows:
```text
Metagenomics_Pipeline_Snakemake/
├── config.yaml          # Configuration file (sample names, paths, parameters)
├── Snakefile            # Core workflow script
├── metagenome_v1.sif     # The Singularity container
└── data/                 # Folder containing raw sequencing reads
    ├── SampleA_1.fastq.gz
    ├── SampleA_2.fastq.gz
    ├── SampleB_1.fastq.gz
    └── ...
## 2. How to Run
Step 1: Dry Run (Sanity Check)
Always perform a dry run first to verify that Snakemake detects the correct number of samples.

```bash
snakemake -np
```
Output check: Look at Job counts.

Step 2: Production Run
CRITICAL: Use the -j flag to limit the total number of CPU cores used by the pipeline. This prevents server overload (Out Of Memory crashes).

Example: If your server has 64 CPUs, set -j 60 to leave a small buffer. Snakemake will automatically manage parallel tasks (e.g., running 2 Assemblies [2x24=48 threads] simultaneously).

## Run in background (nohup)
```bash
nohup snakemake -j 60 --latency-wait 60 --keep-going &
```
-j 60: Max parallel cores.

--latency-wait 60: Waits 60s for files to appear (useful for network storage).

--keep-going: If one sample fails, continue processing the others.

## 3.Output Files
| **Directory**        | **Content**                                                  |
| -------------------- | ------------------------------------------------------------ |
| `assembly/{sample}/` | Final contigs (`final_assembly.fasta`)                       |
| `binning/{sample}/`  | Refined bins (`metawrap_bins/`)                              |
| `drep/`              | Dereplicated genome set & cluster info (`GenomeInformation.csv`) |
| `checkm/`            | Quality assessment results (`bin_stats_ext.tsv`)             |
| `gtdbtk/`            | Taxonomy classification (`gtdbtk.bac120.summary.tsv`)        |
