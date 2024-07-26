#! /bin/bash

# Wrapper script to run LZs in "bulk" by modifying the template SLURM script
#
# Requires a parameter file where it has the relevant parameters for
# /home/rtakei/handy_scripts/locuszooms/scripts/lz_scripts/locus_zoom.cheaha.R
#
# Note that the parameter file needs to be passed to this script in absolute
# path format

mkdir -p /data/scratch/$(whoami)/lz_outputs/{lz_plots,slurm_logs}

PAR_FILE=$1

sed -e "s/USER/$(whoami)/g" -e "s;NUM;$(wc -l ${PAR_FILE} | awk '{print $1}');g" -e "s;FILE;${PAR_FILE};g" /home/$(whoami)/handy_scripts/locuszooms/scripts/slurm/run_lz.template.build38.sl > /data/scratch/$(whoami)/lz_outputs/slurm_logs/run_lz.$(date +%d%b%y.%H%M%p).sl

sbatch /data/scratch/$(whoami)/lz_outputs/slurm_logs/run_lz.$(date +%d%b%y.%H%M%p).sl

