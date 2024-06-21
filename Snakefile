import hashlib
import os
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

import snakemake.utils

from scripts.nanosnk import *


configfile: "config/config.yaml"


snakemake.utils.validate(config, "config/config.schema.json")

# Destination for pipeline output
RESULTS_DIR = config["results_dir"]
# Location of custom scripts
SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")


if sys.version_info < (3, 7):
    # Parts of this script requires on ordered dicts
    abort("Python v3.7 or later required")
elif not Path(config["input_dir"]).is_dir():
    abort("Per sample FAST5 directories not found at", config["input_dir"])

for genome in config["genomes"].values():
    if not Path(genome["minimap2_fasta"]).is_file():
        abort("Reference FASTA not found at", genome["minimap2_fasta"])

    tandem_repeats = genome["sniffles_tandem_repeats"]
    if tandem_repeats and not Path(tandem_repeats).is_file():
        abort("Tandem repeats not found found at", tandem_repeats)

if not Path(config["dorado_models"]).is_file():
    abort("Dorado model table not found at", config["dorado_models"])


GENOMES = config["genomes"]

SAMPLES = collect_samples(
    source=config["input_dir"],
    batch_size=config["batch_size"],
)

# Names sorted to ensure stable input order
GENOME_NAMES = tuple(sorted(GENOMES))
SAMPLE_NAMES = tuple(sorted(SAMPLES))

# Samples for which genotyping should be performed using sniffles
SAMPLES_FOR_GENOTYPING = collect_samples_for_genotyping(
    samples=SAMPLES, blacklist=config["excluded_from_genotyping"]
)

#######################################################################################


def get_fast5_chunk_for_chunk(wildcards):
    return (
        SAMPLES[wildcards.sample].batches[wildcards.batch].chunks[wildcards.hash].files
    )


def get_fastq_chunks_for_batch(wildcards):
    return SAMPLES[wildcards.sample].fastq_chunks(RESULTS_DIR, wildcards.batch)


def get_passed_batches_for_sample(wildcards):
    return SAMPLES[wildcards.sample].bam_batches(RESULTS_DIR, wildcards.genome, "pass")


def get_failed_batches_for_sample(wildcards):
    return SAMPLES[wildcards.sample].bam_batches(RESULTS_DIR, wildcards.genome, "fail")


def get_fasta_from_wildcards(wildcards):
    return GENOMES[wildcards.genome]["minimap2_fasta"]


def get_fasta_mmi_from_wildcards(wildcards):
    return get_fasta_from_wildcards(wildcards) + ".mmi"


def get_tandem_repeats_from_wildcards(wildcards):
    # None is allowed in the config specification, but gets stringified as None in shell
    filepath = GENOMES[wildcards.genome]["sniffles_tandem_repeats"]
    if filepath:
        return os.path.abspath(filepath)

    return ""


#######################################################################################


wildcard_constraints:
    # Genome, sample, and batch names derived from paths should not span directories
    genome="[a-zA-Z0-9_-]+",
    sample="[a-zA-Z0-9_-]+",
    batch="[a-zA-Z0-9_-]+",
    # hashes are all hexencoded SHA256 hashes
    hash="[a-zA-Z0-9]+",
    # Filtered BAMs
    kind="pass|fail",


rule pipeline:
    input:
        # Cohort genotypes from filtered BAMs
        expand(f"{RESULTS_DIR}/genotypes.{{genome}}.vcf.gz", genome=GENOME_NAMES),
        # Alignments against each genome
        expand(
            f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.pass.bam",
            genome=GENOME_NAMES,
            sample=SAMPLE_NAMES,
        ),
        expand(
            f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.fail.bam",
            genome=GENOME_NAMES,
            sample=SAMPLE_NAMES,
        ),
        # BAI files are explicitly required, for samples excluded from genotyping
        expand(
            f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.pass.bam.bai",
            genome=GENOME_NAMES,
            sample=SAMPLE_NAMES,
        ),
        # Alignment metrics
        expand(f"{RESULTS_DIR}/statistics/metrics.{{genome}}.json", genome=GENOME_NAMES),
        expand(
            f"{RESULTS_DIR}/statistics/metrics.{{genome}}.experiments.tsv",
            genome=GENOME_NAMES,
        ),
        expand(
            f"{RESULTS_DIR}/statistics/metrics.{{genome}}.jpg.html",
            genome=GENOME_NAMES,
        ),
        expand(
            f"{RESULTS_DIR}/statistics/metrics.{{genome}}.mapping.tsv",
            genome=GENOME_NAMES,
        ),
        expand(
            f"{RESULTS_DIR}/statistics/metrics.{{genome}}.sequencing.tsv",
            genome=GENOME_NAMES,
        ),


rule dorado_model:
    priority: 50
    group:
        "setup"
    input:
        get_fast5_chunk_for_chunk,
    output:
        batch=f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}/{{hash}}.fast5s.txt",
        model=f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}/{{hash}}.model.txt",
    params:
        models=os.path.abspath(config["dorado_models"]),
        script=os.path.abspath("scripts/select_model.py"),
    resources:
        gpumisc=1,  # Used to limit how many tasks are queud on the GPU node
    shadow:
        "minimal"
    shell:
        """
        realpath {input:q} > {output.batch:q}
        python3 {params.script:q} --file-lists {params.models:q} {output.batch:q} \
            > {output.model:q}
        """


rule dorado:
    priority: 100
    input:
        fast5s=get_fast5_chunk_for_chunk,
        batch=f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}/{{hash}}.fast5s.txt",
        model=f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}/{{hash}}.model.txt",
    output:
        temporary(f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}/{{hash}}.fq.gz"),
    threads: 8
    resources:
        gputask=1,  # Used to limit how many tasks are queud on the GPU node
        slurm_partition="gpuqueue",
        slurm_extra="--gres=gpu:a100:1",
    envmodules:
        "cuda/11.8",
        "dorado/0.3.4",
        "libdeflate/1.18",
        "htslib/1.18",
        "samtools-libdeflate/1.18",
    shadow:
        "minimal"
    shell:
        """
        readonly MODEL=$(cat {input.model:q})
        readonly BATCH="{output}.batch"

        mkdir "${{BATCH}}"
        cat {input.batch:q} | xargs -I "{{}}" ln -s "{{}}" "${{BATCH}}/"

        dorado basecaller --verbose --recursive --device "cuda:all" "${{MODEL}}" ${{BATCH}} \
            | samtools bam2fq -T '*' \
            | bgzip -@ 6 \
            > {output:q}

        rm -r "${{BATCH}}"
        """


rule merge_fastq_chunks:
    input:
        get_fastq_chunks_for_batch,
    output:
        f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}.fq.gz",
    envmodules:
        "python/3.9.16",
    params:
        script=os.path.abspath("scripts/mergebgzip.py"),
    shadow:
        "minimal"
    shell:
        r"""
        python3 {params.script:q} {input:q} > {output:q}
        """


rule minimap2:
    input:
        fa=get_fasta_from_wildcards,
        mmi=get_fasta_mmi_from_wildcards,
        fastq=f"{RESULTS_DIR}/reads/{{sample}}/{{batch}}.fq.gz",
    output:
        passed=temporary(
            f"{RESULTS_DIR}/alignments/{{sample}}/{{batch}}.{{genome}}.pass.bam"
        ),
        failed=temporary(
            f"{RESULTS_DIR}/alignments/{{sample}}/{{batch}}.{{genome}}.fail.bam"
        ),
    params:
        qscore_filter=10,
        script=os.path.abspath("scripts/filter_bam.py"),
    threads: 10
    envmodules:
        "minimap2/2.26",
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shadow:
        "minimal"
    shell:
        """
        minimap2 -y -t 10 -ax map-ont -R '@RG\\tID:{wildcards.batch}\\tSM:{wildcards.sample}' {input.mmi:q} {input.fastq:q} \
        | python3 {params.script:q} \
            --input - \
            --output-fail {output.failed:q} \
            --output-pass - \
            --min-query-length 500 \
            --min-base-quality {params.qscore_filter} \
            --exclude-flags 0x704 \
        | samtools sort -@ 4 -u -T {output.passed:q} \
        | samtools calmd -Q -@ 4 -b - {input.fa:q} \
        > {output.passed:q}
        """


rule merge_bam_pass:
    input:
        get_passed_batches_for_sample,
    output:
        f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.pass.bam",
    threads: 4
    envmodules:
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shadow:
        "minimal"
    shell:
        """
        samtools merge -@ {threads} {output:q} {input:q}
        """


rule merge_bam_fail:
    input:
        get_failed_batches_for_sample,
    output:
        f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.fail.bam",
    envmodules:
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shadow:
        "minimal"
    shell:
        """
        samtools cat -o {output:q} {input:q}
        """


rule index_bam_files:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    threads: 4
    envmodules:
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shadow:
        "minimal"
    shell:
        """
        samtools index {input:q}
        """


rule sniffles2_snf:
    input:
        fa=get_fasta_from_wildcards,
        bam=f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.pass.bam",
        bai=f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.pass.bam.bai",
    output:
        f"{RESULTS_DIR}/genotypes/{{sample}}.{{genome}}.pass.snf",
    threads: 8
    params:
        # None is allowed in the config specification, but gets stringified as None in shell
        tandem_repeats=get_tandem_repeats_from_wildcards,
    envmodules:
        "sniffles/2.3.3",
    shadow:
        "minimal"
    shell:
        """
        args=()
        if [ -n "{params.tandem_repeats}" ]; then
            args+=("--tandem-repeats" {params.tandem_repeats:q})
        fi

        sniffles \
            --threads {threads} \
            --reference "{input.fa:q}" \
            --input "{input.bam:q}" \
            --snf "{output:q}" \
            "${{args[@]}}"
        """


rule sniffles2_tsv:
    input:
        snf=lambda wildcards: expand(
            f"{RESULTS_DIR}/genotypes/{{sample}}.{wildcards.genome}.pass.snf",
            sample=SAMPLES_FOR_GENOTYPING,
        ),
    output:
        tsv=f"{RESULTS_DIR}/genotypes/{{genome}}.tsv",
    run:
        with open(output.tsv, "wt") as handle:
            for filename in input:
                name, _ = os.path.basename(filename).split(".", 1)
                print(os.path.abspath(filename), name, sep="\t", file=handle)


rule sniffles2_vcf:
    input:
        fa=get_fasta_from_wildcards,
        tsv=f"{RESULTS_DIR}/genotypes/{{genome}}.tsv",
        snf=lambda wildcards: expand(
            f"{RESULTS_DIR}/genotypes/{{sample}}.{wildcards.genome}.pass.snf",
            sample=SAMPLES_FOR_GENOTYPING,
        ),
    output:
        f"{RESULTS_DIR}/genotypes.{{genome}}.vcf.gz",
    params:
        tandem_repeats=get_tandem_repeats_from_wildcards,
    threads: 32
    envmodules:
        "sniffles/2.3.3",
    shadow:
        "minimal"
    shell:
        """
        args=()
        if [ -n "{params.tandem_repeats}" ]; then
            args+=("--tandem-repeats" {params.tandem_repeats:q})
        fi

        sniffles \
            --threads {threads} \
            --reference "{input.fa:q}" \
            --input "{input.tsv:q}" \
            --vcf "{output:q}" \
            "${{args[@]}}"
        """


#######################################################################################
## QC


rule qc_metrics:
    input:
        passed=f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.pass.bam",
        failed=f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.fail.bam",
    output:
        f"{RESULTS_DIR}/alignments/{{sample}}.{{genome}}.json",
    params:
        sample=config["qc_sample"],
        script=os.path.abspath("scripts/qc_metrics.py"),
    threads: 4
    shell:
        r"""
        python3 {params.script:q} \
            --sample-size {params.sample} \
            --output {output:q} \
            --threads {threads} \
            {input:q}
        """


rule qc_metrics_join:
    input:
        lambda wildcards: expand(
            f"{RESULTS_DIR}/alignments/{{sample}}.{wildcards.genome}.json",
            sample=SAMPLE_NAMES,
        ),
    output:
        f"{RESULTS_DIR}/statistics/metrics.{{genome}}.json",
    shell:
        r"""
        cat {input:q} > {output:q}
        """


rule qc_metrics_report:
    input:
        f"{RESULTS_DIR}/statistics/metrics.{{genome}}.json",
    output:
        f"{RESULTS_DIR}/statistics/metrics.{{genome}}.experiments.tsv",
        f"{RESULTS_DIR}/statistics/metrics.{{genome}}.jpg.html",
        f"{RESULTS_DIR}/statistics/metrics.{{genome}}.mapping.tsv",
        f"{RESULTS_DIR}/statistics/metrics.{{genome}}.sequencing.tsv",
    params:
        script=os.path.abspath("scripts/qc_report.py"),
        prefix=f"{RESULTS_DIR}/statistics/metrics.{{genome}}",
    shell:
        r"""
        python3 {params.script:q} --image-format jpg \
            {input:q} {params.prefix:q}
        """
