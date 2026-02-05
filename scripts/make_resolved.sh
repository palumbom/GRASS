#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=01:00:00
#SBATCH --job-name=for_eduardo
#SBATCH --chdir=/mnt/home/mpalumbo/work/GRASS
#SBATCH --output=/mnt/home/mpalumbo/work/logs/%x.%j.%a.out

echo "About to start: $SLURM_JOB_NAME"
date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR

echo "About to start Julia"
julia /mnt/home/mpalumbo/work/GRASS/scripts/make_resolved.jl $SLURM_ARRAY_TASK_ID
echo "Julia exited"
date