import hashlib
import os
import sys
from collections import defaultdict
from pathlib import Path

import snakemake.utils

configfile: "config/config.yaml"

snakemake.utils.validate(config, "config/config.schema.json")

# Destination for pipeline output
RESULTS_DIR = config["results_dir"]
# Location of custom scripts
SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")

def abort(*args, **kwargs):
    kwargs["file"] = sys.stderr
    print("ERROR:", *args, **kwargs)
    sys.exit(1)


if sys.version_info < (3, 7):
    # Parts of this script requires on ordered dicts
    abort("Python v3.7 or later required")
elif not Path(config["input_dir"]).is_dir():
    abort("Per sample FAST5 directories not found at", config["input_dir"])
elif not Path(config["minimap2_fasta"]).is_file():
    abort("Reference FASTA not found at", config["minimap2_fasta"])
elif not (Path(config["dorado_model"]) / "config.toml").is_file():
    abort("Dorado model not found at", config["dorado_model"])


#######################################################################################


def fragment(size, items):
    """Faster alternative to grouper for lists/strings."""
    return (items[i : i + size] for i in range(0, len(items), size))


def sha256(items):
    """Calculates the SHA256 hash for a set of values"""
    hasher = hashlib.sha256()
    for it in items:
        hasher.update(it.encode("utf-8"))
    return hasher.hexdigest()


def _collect_fast5s(source, groups=None):
    # fast5s are grouped by dir, to ensure that adding/removing whole folders
    # does not result in the re-processing of existing folders
    if groups is None:
        groups = defaultdict(list)

    for it in source.iterdir():
        if it.is_dir():
            _collect_fast5s(it, groups)
        elif it.suffix.lower() in (".fast5", ".pod5"):
            groups[it.parent].append(str(it))

    return groups.values()


def _collect_samples(destination, source, batch_size=25):
    samples = {}
    for it in Path(source).iterdir():
        samples[it.name] = sample = {}

        for group in _collect_fast5s(it):
            group.sort()  # Ensure stable batches even filesystem order changes

            if batch_size is None:
                sample[sha256(group)] = group
            else:
                for batch in fragment(batch_size, group):
                    sample[sha256(batch)] = batch

    return samples


def _generate_chunks(destination, samples):
    return {
        sample: {
            extension: [
                os.path.join(destination, sample + ".cache", f"{key}.{extension}")
                for key in batches
            ]
            for extension in ("bam", "fq.gz")
        }
        for sample, batches in sorted(samples.items())
    }


SAMPLES = _collect_samples(
    destination=config["results_dir"],
    source=config["input_dir"],
    batch_size=config["batch_size"],
)
CHUNKS = _generate_chunks(destination=config["results_dir"], samples=SAMPLES)

#######################################################################################


wildcard_constraints:
    # Sample names derived from paths should not span directories
    sample="[^/]+",
    # hashes are all hexencoded SHA256 hashes
    hash="[a-zA-Z0-9]+",


rule pipeline:
    input:
        f"{RESULTS_DIR}/genotypes.vcf.gz",
        expand(f"{RESULTS_DIR}/{{sample}}.fq.gz", sample=SAMPLES),
        # MultiQC reports from FastQC reports
        f"{RESULTS_DIR}/statistics/premap/multiqc.html",
        f"{RESULTS_DIR}/statistics/postmap/multiqc.html",
        # Alignment metrics
        f"{RESULTS_DIR}/statistics/metrics.json",


rule dorado_symlinks:
    priority: 50
    group:
        "dorado"
    input:
        lambda wildcards: SAMPLES[wildcards.sample][wildcards.hash],
    output:
        batch=temporary(directory(f"{RESULTS_DIR}/{{sample}}.cache/{{hash}}.fast5s")),
    resources:
        gpuqueue=1,  # Used to limit how many tasks are queud on the GPU node
        slurm_partition="gpuqueue",
    run:
        os.makedirs(output.batch, exist_ok=True)
        for filepath in input:
            src = os.path.realpath(filepath)
            dst = os.path.join(output.batch, os.path.basename(filepath))
            os.symlink(src, dst)


rule dorado:
    priority: 100
    group:
        "dorado"
    input:
        batch=f"{RESULTS_DIR}/{{sample}}.cache/{{hash}}.fast5s",
    output:
        fastq=temporary(f"{RESULTS_DIR}/{{sample}}.cache/{{hash}}.fq.gz"),
    params:
        model=config["dorado_model"],
    threads: 6
    resources:
        gpuqueue=1,  # Used to limit how many tasks are queud on the GPU node
        slurm_partition="gpuqueue",
        slurm_extra="--gres=gpu:a100:1",
    envmodules:
        "cuda/11.8",
        "dorado/0.2.4",
        "libdeflate/1.18",
        "htslib/1.18",
        "samtools-libdeflate/1.18",
    shell:
        """
        dorado basecaller --recursive --device "cuda:all" {params.model} {input.batch} \
            | samtools bam2fq -T '*' \
            | bgzip -@ 4 \
            > {output.fastq}
        """


rule minimap2:
    input:
        fa=config["minimap2_fasta"],
        mmi=config["minimap2_fasta"] + ".mmi",
        fastq=f"{RESULTS_DIR}/{{sample}}.cache/{{hash}}.fq.gz",
    output:
        bam=temporary(f"{RESULTS_DIR}/{{sample}}.cache/{{hash}}.bam"),
    threads: 14
    envmodules:
        "minimap2/2.26",
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shell:
        """
            minimap2 -y -t 10 -ax map-ont -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' {input.mmi} {input.fastq} \
            | samtools sort -@ 4 -u \
            | samtools calmd -Q -@ 4 -b - {input.fa} \
            > {output.bam}
        """


rule merge_fq:
    input:
        fastq=lambda wildcards: CHUNKS[wildcards.sample]["fq.gz"],
    output:
        fastq=f"{RESULTS_DIR}/{{sample}}.fq.gz",
    envmodules:
        "python/3.9.16",
    shell:
        r"""
        python3 scripts/mergebgzip.py {input.fastq} > {output.fastq}
        """


rule merge_bam:
    params:
        qscore_filter=10,
    input:
        bams=lambda wildcards: CHUNKS[wildcards.sample]["bam"],
    output:
        passed=f"{RESULTS_DIR}/{{sample}}.pass.bam",
        passed_csi=f"{RESULTS_DIR}/{{sample}}.pass.bam.csi",
        failed=f"{RESULTS_DIR}/{{sample}}.fail.bam",
        failed_csi=f"{RESULTS_DIR}/{{sample}}.fail.bam.csi",
    threads: 8
    envmodules:
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shell:
        """
        THREADS=1
        if [ {threads} -gt 2 ]; then
            THREADS=$(({threads} / 2))
        fi

        samtools merge -@ ${{THREADS}} -u - {input} \
            | samtools view -@ ${{THREADS}} --write-index -e '[qs] >= {params.qscore_filter}' --output {output.passed} --unoutput {output.failed} -b -
        """


rule sniffles2_snf:
    input:
        fa=config["minimap2_fasta"],
        bam=f"{RESULTS_DIR}/{{sample}}.pass.bam",
    output:
        snf=f"{RESULTS_DIR}/{{sample}}.pass.snf",
    threads: 8
    envmodules:
        "sniffles/2.0.7",
    shell:
        """
        sniffles \
            --threads {threads} \
            --reference "{input.fa}" \
            --input "{input.bam}" \
            --snf "{output.snf}"
        """


# TODO: Belongs in temporary or genotyping folder
rule sniffles2_tsv:
    input:
        snf=expand(f"{RESULTS_DIR}/{{sample}}.pass.snf", sample=SAMPLES),
    output:
        tsv=f"{RESULTS_DIR}/genotypes.tsv",
    run:
        with open(output.tsv, "wt") as handle:
            for filename in input:
                name, _ = os.path.basename(filename).split(".", 1)
                print(filename, name, sep="\t", file=handle)


# TODO: Better name / location? Place in `genotyping/`?
rule sniffles2_vcf:
    input:
        fa=config["minimap2_fasta"],
        tsv=f"{RESULTS_DIR}/genotypes.tsv",
        snf=expand(f"{RESULTS_DIR}/{{sample}}.pass.snf", sample=SAMPLES),
    output:
        vcf=f"{RESULTS_DIR}/genotypes.vcf.gz",
    threads: 8
    shell:
        """
        sniffles \
            --threads {threads} \
            --reference "{input.fa}" \
            --input "{input.tsv}" \
            --vcf "{output.vcf}"
        """


#######################################################################################
## QC


rule fastqc_fastq:
    input:
        fq=f"{RESULTS_DIR}/{{sample}}.fq.gz",
    output:
        html=f"{RESULTS_DIR}/statistics/premap/{{sample}}_fastqc.html",
        data=f"{RESULTS_DIR}/statistics/premap/{{sample}}_fastqc.zip",
    params:
        outdir=os.path.join(RESULTS_DIR, "statistics", "premap"),
        sample=config["qc_sample"],
    envmodules:
        "perl/5.26.3",
        "openjdk/20.0.0",
        "fastqc/0.11.9", # requies perl and openjdk
        "seqtk/1.4",
    shell:
        """
        # -t is used to force the allocation of additional memory
        seqtk sample {input} {params.sample} | fastqc -t 8 -f fastq -o {params.outdir} stdin:{wildcards.sample}
        """


rule fastqc_bam:
    input:
        bam=f"{RESULTS_DIR}/{{sample}}.{{kind}}.bam",
    output:
        html=f"{RESULTS_DIR}/statistics/postmap/{{sample}}_{{kind}}_fastqc.html",
        data=f"{RESULTS_DIR}/statistics/postmap/{{sample}}_{{kind}}_fastqc.zip",
    threads: 4
    params:
        outdir=os.path.join(RESULTS_DIR, "statistics", "postmap"),
        sample=config["qc_sample"],
    envmodules:
        "perl/5.26.3",
        "openjdk/20.0.0",
        "fastqc/0.11.9", # requires perl and openjdk
        "libdeflate/1.18",
        "samtools-libdeflate/1.18",
    shell:
        r"""
        # fastqc -t is used to force the allocation of additional memory
        samtools bam2fq -@ {threads} {input} \
            | seqtk sample - {params.sample} \
            | fastqc -t 8 --java {SCRIPTS_DIR}/java_wrapper.sh -f fastq -o {params.outdir} stdin:{wildcards.sample}_{wildcards.kind}
        """


#######################################################################################


rule multiqc_fastq:
    input:
        zip=expand(
            f"{RESULTS_DIR}/statistics/premap/{{sample}}_fastqc.zip", sample=SAMPLES
        ),
    output:
        html=f"{RESULTS_DIR}/statistics/premap/multiqc.html",
        data=f"{RESULTS_DIR}/statistics/premap/multiqc_data.zip",
    # FIXME: Needs envmodule
    shell:
        r"""
        multiqc --zip-data-dir \
            --filename {output.html} \
            {input.zip}
        """


rule multiqc_bam:
    input:
        zip=expand(
            f"{RESULTS_DIR}/statistics/postmap/{{sample}}_{{kind}}_fastqc.zip",
            sample=SAMPLES,
            kind=("pass", "fail"),
        ),
    output:
        html=f"{RESULTS_DIR}/statistics/postmap/multiqc.html",
        data=f"{RESULTS_DIR}/statistics/postmap/multiqc_data.zip",
    # FIXME: Needs envmodule
    shell:
        r"""
        multiqc --zip-data-dir \
            --filename {output.html} \
            {input.zip}
        """


#######################################################################################


rule qc_stats:
    input:
        passed=f"{RESULTS_DIR}/{{sample}}.pass.bam",
    output:
        json=temporary(f"{RESULTS_DIR}/{{sample}}.cache/metrics.json"),
    envmodules:
        "python/3.9.16",
    params:
        sample=config["qc_sample"],
    shell:
        r"""
        python3 scripts/qc_metrics.py stats --nsample {params.sample} \
            {input:q} \
            > {output:q}
        """


rule qc_stats_join:
    input:
        expand(f"{RESULTS_DIR}/{{sample}}.cache/metrics.json", sample=SAMPLES),
    output:
        f"{RESULTS_DIR}/statistics/metrics.json",
    envmodules:
        "python/3.9.16",
    shell:
        r"""
        python3 scripts/qc_metrics.py join \
            {input:q} \
            > {output:q}
        """
