#!/bin/bash
#SBATCH --account=open
#SBATCH --time=48:00:00
#SBATCH --nodes=15
#SBATCH --ntasks=15
#SBATCH --mem=40gb
#SBATCH --job-name=exp_meter_proj_RV
#SBATCH --output /storage/home/efg5335/work/GRASS/exp_meter_proj_RV.out
#SBATCH --error /storage/home/efg5335/work/GRASS/exp_meter_proj_RV.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/GRASS
julia projected_RV.jl

##now=$(date +"%T")
##echo "end time ... $now"