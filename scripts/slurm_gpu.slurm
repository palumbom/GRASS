#!/bin/bash
#SBATCH -A dfc13_mri
#SBATCH -p mgc-mri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=24:00:00
#SBATCH --job-name=N_convergence
#SBATCH --chdir=/storage/home/efg5335/work/GRASS
#SBATCH --output=/storage/home/efg5335/work/GRASS/N_convergence.out

echo "About to start: $SLURM_JOB_NAME"
date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR
echo "About to start Julia"
julia NEID_GRASS_all_lines.jl
echo "Julia exited"
date
