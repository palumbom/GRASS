#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --time=48:00:00
#SBATCH --nodes=20
#SBATCH --ntasks=20
#SBATCH --mem=25gb
#SBATCH --job-name=ext_reduction
#SBATCH --output /storage/home/efg5335/work/GRASS/ext_reduction.out
#SBATCH --error /storage/home/efg5335/work/GRASS/ext_reduction.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia ext_reduction.jl

##now=$(date +"%T")
##echo "end time ... $now"