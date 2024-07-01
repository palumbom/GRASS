#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=15gb
#SBATCH --job-name=NEID_GRASS_all_lines_off_full
#SBATCH --output /storage/home/efg5335/work/GRASS/NEID_GRASS_all_lines_off_full.out
#SBATCH --error /storage/home/efg5335/work/GRASS/NEID_GRASS_all_lines_off_full.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia NEID_GRASS_all_lines.jl

##now=$(date +"%T")
##echo "end time ... $now"