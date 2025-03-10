#!/bin/bash
#SBATCH --job-name=t3_2 # Job name
#SBATCH --output=data/out/t3_2_%a.out   # Standard output and error log
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --partition=icelake             # Partition to submit to
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --cpus-per-task=2         # Number of CPU cores per task
#SBATCH --mem=24000                     # Memory per node in MB
#SBATCH --mail-type=BEGIN,FAIL,END      # Notifications for job done & fail
#SBATCH --mail-user=sb2690@cam.ac.uk    # Email to which notifications are sent
#SBATCH --account=AMWOOD-SL3-CPU        # Account name
#SBATCH --array=1-12              # Number of jobs to run

# Specify the path to the config file
config=/rds/user/sb2690/hpc-work/R-multi-state-modelling/data/config/bench_test3_2.tsv

# Extract the parameters for the current $SLURM_ARRAY_TASK_ID
ptest=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
pmodel=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
pnelements=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
pncovariates=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
pntrain=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
pntest=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
pncpu=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)
prep=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $config)

source ~/.bashrc
module load R/4.3.1-icelake
micromamba activate msm
Rscript benchmark.R ${ptest} ${pmodel} ${pnelements} ${pncovariates} ${pntrain} ${pntest} ${pncpu} ${prep}
