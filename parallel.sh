#!/bin/bash
#SBATCH --job-name=MSI_power
#SBATCH --account=shenx
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --array=1            # single task is fine
#SBATCH --output=Results/Run_%a.slurm.out
#SBATCH --error=Results/Run_%a.slurm.err
#SBATCH --mail-user=pyrrhanikoss@outlook.com
#SBATCH --mail-type=END

module load R/4.4.0-openblas-rocky8

# guarantee paths exist
cd /panfs/jay/groups/27/shenx/zhao1118/distributed_algorithm || exit 1
mkdir -p Results

# run, writing R’s own output to a separate file
R CMD BATCH --vanilla power_final.R Results/Run_${SLURM_ARRAY_TASK_ID}.Rout
