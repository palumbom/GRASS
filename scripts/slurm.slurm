#!/bin/bash
#SBATCH --account=open
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=50gb
#SBATCH --job-name=set_geometry
#SBATCH --output /storage/home/efg5335/work/GRASS/set_geometry.out
#SBATCH --error /storage/home/efg5335/work/GRASS/set_geometry.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia set_geometry.jl

##now=$(date +"%T")
##echo "end time ... $now"