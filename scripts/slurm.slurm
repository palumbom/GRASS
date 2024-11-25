#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=10gb
#SBATCH --job-name=boulder_october_eclipse
#SBATCH --output /storage/home/efg5335/work/GRASS/boulder_october_eclipse.out
#SBATCH --error /storage/home/efg5335/work/GRASS/boulder_october_eclipse.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia set_geometry.jl

##now=$(date +"%T")
##echo "end time ... $now"