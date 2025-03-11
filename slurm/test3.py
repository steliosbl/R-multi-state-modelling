import os
import pandas as pd
from pathlib import Path


def print_file(cpus, row_start, row_stop, filename="msm_benchmark_parallel.slurm"):
    template_params = f"""#!/bin/bash
#SBATCH --job-name=t3_{cpus} # Job name
#SBATCH --output=data/out/t3_{cpus}_%a.out   # Standard output and error log
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --partition=icelake-himem             # Partition to submit to
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --cpus-per-task={cpus}         # Number of CPU cores per task
#SBATCH --mem={6760*max(cpus, 4)}             # Memory per node in MB
#SBATCH --mail-type=BEGIN,FAIL,END      # Notifications for job done & fail
#SBATCH --mail-user=sb2690@cam.ac.uk    # Email to which notifications are sent
#SBATCH --account=AMWOOD-SL3-CPU        # Account name
#SBATCH --array={row_start}-{row_stop}              # Number of jobs to run

# Specify the path to the config file
config=/rds/user/sb2690/hpc-work/R-multi-state-modelling/slurm/temp/test3_{cpus}.tsv
"""
    template_body = """
# Extract the parameters for the current $SLURM_ARRAY_TASK_ID
ptest=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
pmodel=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
pnelements=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
pncovariates=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
pntrain=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
pntest=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
ptimes=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)
pncpu=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $config)
prep=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $10}' $config)

cd /rds/user/sb2690/hpc-work/R-multi-state-modelling
source ~/.bashrc
module load R/4.3.1-icelake
micromamba activate msm
Rscript benchmark.R ${ptest} ${pmodel} ${pnelements} ${pncovariates} ${pntrain} ${pntest} ${ptimes} ${pncpu} ${prep}
"""
    with open(filename, "w") as file:
        file.write(template_params)
        file.write(template_body)


base_dir = Path("/rds/user/sb2690/hpc-work/R-multi-state-modelling")
temp_dir = base_dir / "slurm/temp"
config = pd.read_csv(base_dir / "data/config/test3.tsv", sep="\t")

for cpus in config["cpus"].unique():
    subset = config[config["cpus"] == cpus].sort_values("row")
    row_start, row_stop = subset.row.values[0], subset.row.values[-1]

    subset_filename = temp_dir / f"test3_{cpus}.tsv"
    subset.to_csv(subset_filename, sep="\t", index=False)

    slurm_filename = temp_dir / f"test3_{cpus}.slurm"
    print_file(cpus, row_start, row_stop, slurm_filename)
    os.system(f"sbatch {slurm_filename}")
