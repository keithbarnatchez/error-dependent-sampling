#!/bin/bash
#SBATCH -J error-dep-causal
#SBATCH -o error-dep-causal-output.txt
#SBATCH -e error-dep-causal-errors.txt
#SBATCH -t 70:55:00
#SBATCH --priority=4294967293
#SBATCH --mem=16G
#SBATCH --no-requeue
#SBATCH -n 4
#SBATCH --mail-user=keithbarnatchez@g.harvard.edu
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2

module load R

Rscript sim-main.R
