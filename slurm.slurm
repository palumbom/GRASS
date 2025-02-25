#!/bin/bash
#SBATCH --account=open
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=10gb
#SBATCH --job-name=projected_RV
#SBATCH --output /storage/home/efg5335/work/GRASS/projected_RV.out
#SBATCH --error /storage/home/efg5335/work/GRASS/projected_RV.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia projected_RV.jl

##now=$(date +"%T")
##echo "end time ... $now"