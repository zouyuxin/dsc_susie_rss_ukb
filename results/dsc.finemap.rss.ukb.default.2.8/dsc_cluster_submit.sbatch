#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=dsc_finemap_rss_ukb
#SBATCH --mail-type=BEGIN,END,FAIL

# load the openmpi module
module load R

dsc susie_rss_ukb.dsc --host yuxin_midway.yml -c 4 &> dsc.finemap.rss.ukb.default.job.log
