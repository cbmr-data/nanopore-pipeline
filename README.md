# Nanopore pipeline

This repository contains a work-in-progress pipeline for genotyping Oxford Nanopore data using [Dorado](https://github.com/nanoporetech/dorado) and [sniffles](https://github.com/fritzsedlazeck/Sniffles).

## Prerequisites

This pipeline is designed to run on the [Esrum cluster](https://cbmr-data.github.io) using Slurm for task execution and using envmodules for software requirements. Modifying `profile/config.yaml` before running the pipeline will be required. See `Snakefile` for the expected program versions.

The various scripts rely on a number of python modules. It is recommended to use [uv](https://github.com/astral-sh/uv) to install these:

    $ pip install uv
    $ uv venv
    $ uv pip install -r requirements.txt

## Input data

The `projects/example` folder shows the expected layout for sequencing projects:

- A `data` folder containing any number of subfolders with `fast5` or `pod5` files. Data is automatically batched by folder. Can be overridden in `config.yaml`.
- A `results` folder containing the output from the pipeline. Can be overridden in `config.yaml`.
- A `config.yaml` file containing per-project settings and the reference genome(s) to use. See `config/config.yaml` for a template.
- A `models.tsv` file containing a table of dorado models to be used depending on the sample rate and bases per second values of the input data. These parameters are identified using the `select_model.py` script.

## Usage

To run the `example` project:

    $ . .venv/bin/activate.sh
    $ python3 run_pipeline.py example

Note that this project does not include actual data and therefore cannot be run to completion.
