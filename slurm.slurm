#!/bin/bash
#SBATCH --account=ebf11_c
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=10gb
#SBATCH --job-name=NEID_RV_ccf
#SBATCH --output /storage/home/efg5335/work/GRASS/NEID_RV_ccf.out
#SBATCH --error /storage/home/efg5335/work/GRASS/NEID_RV_ccf.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia NEID_RV_ccf.jl

##now=$(date +"%T")
##echo "end time ... $now"