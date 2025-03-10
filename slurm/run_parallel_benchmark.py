import os
import pandas as pd
from pathlib import Path


def print_file(n_cpu, n_jobs, filename="msm_benchmark_parallel.slurm"):
    template_params = f"""#!/bin/bash
#SBATCH --job-name=t3_{n_cpu} # Job name
#SBATCH --output=data/out/t3_{n_cpu}_%a.out   # Standard output and error log
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --partition=icelake-himem             # Partition to submit to
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --cpus-per-task={n_cpu}         # Number of CPU cores per task
#SBATCH --mem={6760*n_cpu}              # Memory per node in MB
#SBATCH --mail-type=BEGIN,FAIL,END      # Notifications for job done & fail
#SBATCH --mail-user=sb2690@cam.ac.uk    # Email to which notifications are sent
#SBATCH --account=AMWOOD-SL3-CPU        # Account name
#SBATCH --array=1-{n_jobs}              # Number of jobs to run

# Specify the path to the config file
config=/rds/user/sb2690/hpc-work/R-multi-state-modelling/data/config/bench_test3_{n_cpu}.tsv
"""
    template_body = """
# Extract the parameters for the current $SLURM_ARRAY_TASK_ID
ptest=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
pmodel=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
pnelements=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
pncovariates=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
pntrain=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
pntest=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
pncpu=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)
prep=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $config)

cd /rds/user/sb2690/hpc-work/R-multi-state-modelling
source ~/.bashrc
module load R/4.3.1-icelake
micromamba activate msm
Rscript benchmark.R ${ptest} ${pmodel} ${pnelements} ${pncovariates} ${pntrain} ${pntest} ${pncpu} ${prep}
"""
    with open(filename, "w") as file:
        file.write(template_params)
        file.write(template_body)


config_dir = Path("/rds/user/sb2690/hpc-work/R-multi-state-modelling/data/config")
configs = [_ for _ in os.listdir(config_dir) if _.endswith(".tsv") and "test3_" in _]
for config in configs:
    contents = pd.read_csv(config_dir / config, sep="\t")
    n_jobs = contents.shape[0]
    n_cpu = contents["n_cpu"].unique()[0]

    print_file(n_cpu, n_jobs)
    os.system("sbatch msm_benchmark_parallel.slurm")
