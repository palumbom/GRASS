#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=15gb
#SBATCH --job-name=neid_RVlinebyline_Kostogryz_LD_300
#SBATCH --output /storage/home/efg5335/work/GRASS/neid_RVlinebyline_Kostogryz_LD_300.out
#SBATCH --error /storage/home/efg5335/work/GRASS/neid_RVlinebyline_Kostogryz_LD_300.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia NEID_RV_ccf.jl

##now=$(date +"%T")
##echo "end time ... $now"