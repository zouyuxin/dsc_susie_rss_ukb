#!/bin/bash

#SBATCH --time=5:00:00
#SBATCH --partition=mstephens
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=pip_compare_initNA_lamb1e-4
#SBATCH --mail-type=BEGIN,END,FAIL

# load the openmpi module
module load R

Rscript pip_compare_initNA_lamb1e-4.R
