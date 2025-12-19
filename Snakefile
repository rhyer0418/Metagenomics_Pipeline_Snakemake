import os

# Load configuration
configfile: "config.yaml"

# Define output root directory
OUT = config["outdir"]

# === 1. Parse Samples ===
# Detect samples from the data directory based on R1 suffix
SAMPLES, = glob_wildcards(os.path.join(config["samples_dir"], "{sample}" + config["r1_suffix"]))
print(f"Detected {len(SAMPLES)} samples: {SAMPLES}")

# === 2. Target Files ===
rule all:
    input:
        # QC results
        expand(f"{OUT}/QC/{{sample}}/clean_1.fastq", sample=SAMPLES),
        # dRep results
        f"{OUT}/dRep/data_tables/Wdb.csv",
        # GTDB-Tk results
        f"{OUT}/gtdb/gtdbtk.bac120.summary.tsv",
        # Annotation results
        f"{OUT}/microbeannotator_results/metabolic_summary__module_completeness.tab",
        # HUMAnN3 results
        expand(f"{OUT}/humann/{{sample}}_merged_genefamilies.tsv", sample=SAMPLES),
        # Phylophlan results
        f"{OUT}/Phylogeny/phylophlan_output/RAxML_bestTree.dereplicated_genomes_refined.tre"

# === 3. Quality Control (MetaWRAP) ===
rule metawrap_qc:
    input:
        r1 = os.path.join(config["samples_dir"], "{sample}" + config["r1_suffix"]),
        r2 = os.path.join(config["samples_dir"], "{sample}" + config["r2_suffix"])
    output:
        c1 = f"{OUT}/QC/{{sample}}/clean_1.fastq",
        c2 = f"{OUT}/QC/{{sample}}/clean_2.fastq"
    params:
        outdir = f"{OUT}/QC/{{sample}}",
        env = config["envs"]["metawrap"],
        skip = "--skip-bmtagger" if config.get("skip_bmtagger", True) else ""
    threads: config["threads"]["qc"]
    benchmark: f"{OUT}/benchmarks/qc/{{sample}}.txt"
    shell:
        """
        # Disable strict variable checking to prevent Conda script errors
        set +u
        source {params.env}
        set -u
        set -e
        
        if [ -d "{params.outdir}" ]; then rm -rf {params.outdir}; fi
        mkdir -p {params.outdir}

        start_time=$(date +%s)
        echo "[Time Log] QC for {wildcards.sample} STARTED at $(date)"

        metawrap read_qc -1 {input.r1} -2 {input.r2} -t {threads} -o {params.outdir} {params.skip}
  
        mv results/QC/{wildcards.sample}/final_pure_reads_1.fastq \
           results/QC/{wildcards.sample}/clean_1.fastq

        mv results/QC/{wildcards.sample}/final_pure_reads_2.fastq \
           results/QC/{wildcards.sample}/clean_2.fastq

        end_time=$(date +%s)
        echo "[Time Log] QC for {wildcards.sample} FINISHED at $(date)"
        echo "[Time Log] Duration: $((end_time - start_time)) seconds"
        """

# === 4. Assembly (MetaSPAdes) ===
rule assembly:
    input:
        c1 = f"{OUT}/QC/{{sample}}/clean_1.fastq",
        c2 = f"{OUT}/QC/{{sample}}/clean_2.fastq"
    output:
        contigs = f"{OUT}/Assembly/{{sample}}/final_assembly.fasta"
    params:
        outdir = f"{OUT}/Assembly/{{sample}}",
        env = config["envs"]["metawrap"]
    threads: config["threads"]["assembly"]
    resources:
        mem_mb = 100000
    benchmark: f"{OUT}/benchmarks/assembly/{{sample}}.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        if [ -d "{params.outdir}" ]; then rm -rf {params.outdir}; fi
        mkdir -p {params.outdir}

        start_time=$(date +%s)
        echo "[Time Log] Assembly for {wildcards.sample} STARTED at $(date)"

        metawrap assembly -1 {input.c1} -2 {input.c2} -m 100 -t {threads} --metaspades -o {params.outdir}

        end_time=$(date +%s)
        echo "[Time Log] Assembly for {wildcards.sample} FINISHED at $(date)"
        echo "[Time Log] Duration: $((end_time - start_time)) seconds"
        """

# === 5. Binning ===
rule binning:
    input:
        contigs = f"{OUT}/Assembly/{{sample}}/final_assembly.fasta",
        c1 = f"{OUT}/QC/{{sample}}/clean_1.fastq",
        c2 = f"{OUT}/QC/{{sample}}/clean_2.fastq"
    output:
        directory(f"{OUT}/Bins/{{sample}}/metabat2_bins"),
        directory(f"{OUT}/Bins/{{sample}}/maxbin2_bins"),
        directory(f"{OUT}/Bins/{{sample}}/concoct_bins")
    params:
        outdir = f"{OUT}/Bins/{{sample}}",
        env = config["envs"]["metawrap"]
    threads: config["threads"]["binning"]
    benchmark: f"{OUT}/benchmarks/binning/{{sample}}.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        if [ -d "{params.outdir}" ]; then rm -rf {params.outdir}; fi
        mkdir -p {params.outdir}

        start_time=$(date +%s)
        echo "[Time Log] Binning for {wildcards.sample} STARTED at $(date)"

        metawrap binning -o {params.outdir} -t {threads} -a {input.contigs} \
            --metabat2 --maxbin2 --concoct {input.c1} {input.c2}

        end_time=$(date +%s)
        echo "[Time Log] Binning for {wildcards.sample} FINISHED at $(date)"
        echo "[Time Log] Duration: $((end_time - start_time)) seconds"
        """

# === 6. Bin Refinement ===
rule bin_refinement:
    input:
        mb2 = f"{OUT}/Bins/{{sample}}/metabat2_bins",
        max = f"{OUT}/Bins/{{sample}}/maxbin2_bins",
        con = f"{OUT}/Bins/{{sample}}/concoct_bins"
    output:
        directory(f"{OUT}/bin_refinement/{{sample}}/metawrap_50_5_bins")
    params:
        outdir = f"{OUT}/bin_refinement/{{sample}}",
        env = config["envs"]["metawrap"]
    threads: config["threads"]["refinement"]
    benchmark: f"{OUT}/benchmarks/refinement/{{sample}}.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        if [ -d "{params.outdir}" ]; then rm -rf {params.outdir}; fi
        mkdir -p {params.outdir}

        start_time=$(date +%s)
        echo "[Time Log] Refinement for {wildcards.sample} STARTED at $(date)"

        metawrap bin_refinement -o {params.outdir} -t {threads} \
            -A {input.mb2} -B {input.max} -C {input.con} -c 50 -x 5

        end_time=$(date +%s)
        echo "[Time Log] Refinement for {wildcards.sample} FINISHED at $(date)"
        echo "[Time Log] Duration: $((end_time - start_time)) seconds"
        """

# === 7. dRep Dereplication (Aggregation) ===
rule run_drep:
    input:
        # Wait for all samples to finish refinement
        bins = expand(f"{OUT}/bin_refinement/{{sample}}/metawrap_50_5_bins", sample=SAMPLES)
    output:
        f"{OUT}/dRep/data_tables/Wdb.csv",
        directory(f"{OUT}/dRep/dereplicated_genomes")
    params:
        outdir = f"{OUT}/dRep",
        env = config["envs"]["drep"]
    threads: config["threads"]["drep"]
    benchmark: f"{OUT}/benchmarks/drep/drep_all.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        if [ -d "{params.outdir}" ]; then rm -rf {params.outdir}; fi
        mkdir -p {params.outdir}

        start_time=$(date +%s)
        echo "[Time Log] dRep STARTED at $(date)"

        dRep dereplicate {params.outdir} \
            -g {OUT}/bin_refinement/*/metawrap_50_5_bins/*.fa \
            -sa 0.99 -pa 0.9 -nc 0.30 -p {threads} -d -comp 50 -con 5

        end_time=$(date +%s)
        echo "[Time Log] dRep FINISHED at $(date)"
        echo "[Time Log] Duration: $((end_time - start_time)) seconds"
        """

# === 8. GTDB-Tk Taxonomy ===
rule gtdbtk:
    input: f"{OUT}/dRep/dereplicated_genomes"
    output: f"{OUT}/gtdb/gtdbtk.bac120.summary.tsv"
    params:
        outdir = f"{OUT}/gtdb",
        db = config["databases"]["gtdbtk"],
        env = config["envs"]["gtdbtk"]
    threads: config["threads"]["gtdbtk"]
    benchmark: f"{OUT}/benchmarks/gtdb/gtdbtk.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.db}

        start_time=$(date +%s)
        echo "[Time Log] GTDB-Tk STARTED at $(date)"

        gtdbtk classify_wf --cpus {threads} --genome_dir {input} \
            --out_dir {params.outdir} --extension fa --force

        end_time=$(date +%s)
        echo "[Time Log] GTDB-Tk FINISHED at $(date)"
        echo "[Time Log] Duration: $((end_time - start_time)) seconds"
        """

# === 9. Rename Contigs (SeqKit) ===
rule rename_contigs:
    input: f"{OUT}/dRep/dereplicated_genomes"
    output: directory(f"{OUT}/dRep/clean_bins")
    benchmark: f"{OUT}/benchmarks/seqkit/rename.txt"
    shell:
        """
        set -euo pipefail
        
        mkdir -p {output}

        for f in {input}/*.fa; do
            base=$(basename "$f" .fa)
            seqkit replace -p ".+" -r "${{base}}_contig_{{nr}}" "$f" > {output}/${{base}}.fa
        done
        """

# === 10. Functional Annotation (MetaWRAP) ===
rule annotate_bins:
    input: bins = f"{OUT}/dRep/clean_bins"
    output:
        directory(f"{OUT}/function/bin_translated_genes"),
        directory(f"{OUT}/function/bin_untranslated_genes")
    params: 
        outdir = f"{OUT}/function", 
        env = config["envs"]["metawrap"]
    threads: 15
    benchmark: f"{OUT}/benchmarks/function/annotate.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        if [ -d "{params.outdir}" ]; then rm -rf {params.outdir}; fi
        mkdir -p {params.outdir}

        metaWRAP annotate_bins -o {params.outdir} -t {threads} -b {input.bins}
        """

# === 11. MicrobeAnnotator ===
rule microbeannotator:
    input: proteins = f"{OUT}/function/bin_translated_genes"
    output: f"{OUT}/microbeannotator_results/metabolic_summary__module_completeness.tab"
    params:
        outdir = f"{OUT}/microbeannotator_results",
        db = config["databases"]["microbeannotator"],
        env = config["envs"]["microbeannotator"]
    threads: config["threads"]["microbeannotator"]
    benchmark: f"{OUT}/benchmarks/function/microbeannotator.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        mkdir -p {params.outdir}

        microbeannotator -i {input.proteins}/*.faa -d {params.db} -o {params.outdir} \
            -m diamond -p {threads} -t 4 --refine --light
        """

# === 12. HUMAnN3 Profiling ===
rule humann3:
    input:
        c1 = f"{OUT}/QC/{{sample}}/clean_1.fastq",
        c2 = f"{OUT}/QC/{{sample}}/clean_2.fastq"
    output: 
        gene = f"{OUT}/humann/{{sample}}_merged_genefamilies.tsv",
        path = f"{OUT}/humann/{{sample}}_merged_pathabundance.tsv",
        cov  = f"{OUT}/humann/{{sample}}_merged_pathcoverage.tsv"
    params:
        outdir = f"{OUT}/humann",
        env = config["envs"]["humann3"],
        mpadb = config["databases"]["metaphlan"],
        mpaidx = config["databases"]["metaphlan_index"]
    threads: config["threads"]["humann"]
    benchmark: f"{OUT}/benchmarks/humann/{{sample}}.txt"
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e
        
        mkdir -p {params.outdir}

        # Merge paired-end reads for HUMAnN3
        cat {input.c1} {input.c2} > {params.outdir}/{wildcards.sample}_merged.fq

        humann --input {params.outdir}/{wildcards.sample}_merged.fq \
            --threads {threads} \
            --output {params.outdir} \
            --metaphlan-options '--bowtie2db {params.mpadb} --index {params.mpaidx} --offline'

        rm {params.outdir}/{wildcards.sample}_merged.fq
        """

# === 13. Phylophlan Phylogeny ===
rule phylophlan:
    input: 
        genomes = f"{OUT}/dRep/dereplicated_genomes"
    output:
        tree = f"{OUT}/Phylogeny/phylophlan_output/RAxML_bestTree.dereplicated_genomes_refined.tre"
    params:
        config_file = config["phylophlan"]["config"],
        db = config["phylophlan"]["db"],
        env = config["envs"]["phylophlan"],
        out_dir = f"{OUT}/Phylogeny/phylophlan_output"
    threads: config["threads"]["phylophlan"]
    shell:
        """
        set +u
        source {params.env}
        set -u
        set -e

        mkdir -p {params.out_dir}

        phylophlan \
            -i $(readlink -f {input.genomes}) \
            -d phylophlan\
            -f $(readlink -f {params.config_file}) \
            --databases_folder {params.db} \
            -t a \
            --diversity high \
            --fast \
            -o {params.out_dir} \
            --nproc {threads} \
            --genome_extension .fa \
            --verbose
        """
