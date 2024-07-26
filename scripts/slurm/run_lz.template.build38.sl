#!/bin/bash
#
#SBATCH --job-name=locuszoom_b38
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --time=02:00:00
#SBATCH --output=/data/scratch/USER/lz_outputs/slurm_logs/%x_%A_%a.out
#SBATCH --error=/data/scratch/USER/lz_outputs/slurm_logs/%x_%A_%a.err
#SBATCH --array=1-NUM

module load R/4.3.1-foss-2022b

mkdir -p /data/scratch/USER/lz_outputs/lz_plots

INPAR=($(cat FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1))

Rscript /home/rtakei/handy_scripts/locuszooms/scripts/lz_scripts/locus_zoom.cheaha.build38.R ${INPAR[@]}

