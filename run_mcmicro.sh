#!/usr/bin/env bash

# activate the conda environment set up for mcmicro
conda activate mcmicro

# run mcmicro in the directory containing files
# directory structure should contain markers.csv, params.yml, and subdirectory raw/ with separate .tif for each channel
nextflow run labsyspharm/mcmicro --in /home/austin/mcmicro/L2-Section2/
