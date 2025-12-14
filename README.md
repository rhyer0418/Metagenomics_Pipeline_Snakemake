# Automated Metagenomics Pipeline (Snakemake)
**Version:** v1  
**Description:** Auto pipeline for metagenomics binning, refinement, taxonomy, and annotation.

## Software
| **Software**        | **Website**                 |
| --------------------         | -------------------------------------|
| `metaWRAP(v1.3.2)`            | https://github.com/bxlab/metaWRAP   |
| `metaSPAdes (v4.2.0)`         | https://github.com/ablab/spades  |
| `MEGAHIT (v1.1.3) `           |https://github.com/voutcn/megahit|
| `MaxBin2 (v2.2.7`             | https://github.com/assemblerflow/flowcraft/blob/master/docs/user/components/maxbin2.rst |
| `metaBAT2 (v2.12.1)`          | https://bitbucket.org/berkeleylab/metabat/src/master/|
| `CONCOCT (v1.1.0)`            |https://github.com/BinPro/CONCOCT|
| `dRep (v3.6.2)`               | https://github.com/MrOlm/drep |
| `GTDB-tk (v2.1.1)`            | http://gtdb.ecogenomic.org/ |
| `Salmon (v1.10.3)`            | https://github.com/COMBINE-lab/salmon |
| `PhyloPhlAn (v3.1.68)`        | https://github.com/biobakery/phylophlan|
| `MicrobeAnnotator (v2.0.5)`   | https://github.com/cruizperez/MicrobeAnnotator |



## 1. Project Directory Structure
Before running, organize your working directory as follows:
```text
Metagenomics_Pipeline_Snakemake/
├── config.yaml          # Configuration file (sample names, paths, parameters)
├── Snakefile            # Core workflow script
##├── metagenome_v1.sif     # The Singularity container
└── data/                 # Folder containing raw sequencing reads
    ├── SampleA_1.fq.gz
    ├── SampleA_2.fq.gz
    ├── SampleB_1.fq.gz
    └── ...
```

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
nohup snakemake --cores 40 --resources mem_mb=110000 --latency-wait 60 --keep-going --rerun-incomplete &
```

--latency-wait 60: Waits 60s for files to appear (useful for network storage).

--keep-going: If one sample fails, continue processing the others.

## 3.Output Files
```
results/
├── QC/
│   ├── sample1/
│   │   ├── clean_1.fastq
│   │   └── clean_2.fastq
│   ├── sample2/
│   │   ├── clean_1.fastq
│   │   └── clean_2.fastq
│   └── ...
├── Assembly/
│   ├── sample1/
│   │   └── final_assembly.fasta
│   ├── sample2/
│   │   └── final_assembly.fasta
│   └── ...
├── Bins/
│   ├── sample1/
│   │   ├── metabat2_bins/
│   │   ├── maxbin2_bins/
│   │   └── concoct_bins/
│   ├── sample2/
│   │   ├── metabat2_bins/
│   │   ├── maxbin2_bins/
│   │   └── concoct_bins/
│   └── ...
├── bin_refinement/
│   ├── sample1/
│   │   └── metawrap_50_5_bins/
│   ├── sample2/
│   │   └── metawrap_50_5_bins/
│   └── ...
├── dRep/
│   ├── data_tables/
│   │   └── Wdb.csv
│   ├── dereplicated_genomes/
│   │   └── *.fa
│   └── clean_bins/
│       └── *.fa (rename bins)
├── gtdb/
│   └── gtdbtk.bac120.summary.tsv
├── function/
│   ├── bin_translated_genes/
│   │   └── *.faa
│   └── bin_annotations/
│       └── *.txt 
├── microbeannotator_results/
│   └── annotation_results.txt
├── humann/
│   ├── sample1_genefamilies.tsv
│   ├── sample2_genefamilies.tsv
│   └── ...
└── Phylogeny/
    └── phylophlan_output/
        └── phylogenetic_tree_files
```
        
| **Directory**        | **Content**                                                  |
| -------------------- | ------------------------------------------------------------ |
| `Assembly/{sample}/` | Final contigs (`final_assembly.fasta`)                       |
| `Bins/{sample}/`  | Refined bins (`metawrap_bins/`)                              |
| `checkm/`            | Quality assessment results (`bin_stats_ext.tsv`)             |
| `gtdb/`            | Taxonomy classification (`gtdbtk.bac120.summary.tsv`)        |
