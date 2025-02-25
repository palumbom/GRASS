#!/bin/bash
#SBATCH -A dfc13_mri
#SBATCH -p mgc-mri
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --gpus-per-task=2
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=48:00:00
#SBATCH --job-name=SSD_2
#SBATCH --chdir=/storage/home/efg5335/work/GRASS
#SBATCH --output=/storage/home/efg5335/work/GRASS/icymoons.out

echo "About to start: $SLURM_JOB_NAME"
date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR
echo "About to start Julia"
julia NEID_GRASS_all_lines.jl
echo "Julia exited"
date
