#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=10gb
#SBATCH --job-name=neid_all_lines_rv_off_KSSD_ext
#SBATCH --output /storage/home/efg5335/work/GRASS/neid_all_lines_rv_off_KSSD_ext.out
#SBATCH --error /storage/home/efg5335/work/GRASS/neid_all_lines_rv_off_KSSD_ext.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia NEID_GRASS_all_lines.jl

##now=$(date +"%T")
##echo "end time ... $now"