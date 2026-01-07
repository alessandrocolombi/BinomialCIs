#!/bin/bash
#SBATCH --job-name=Exp_n1small
#SBATCH --time=96:00:00                   # Maximum wall time (hh:mm:ss)
#SBATCH --nodes=1                         # Request one node
#SBATCH --ntasks=1                        # One task (process) total
#SBATCH --cpus-per-task=34                # One CPU core per task
#SBATCH --partition=g100_usr_pmem         # Partition (queue) to submit to
#SBATCH --qos=g100_qos_lprod              # Quality of Service
#SBATCH --mem=32G                         # Memory per node (adjust as needed)
#SBATCH --account=uBG25_Argiento          # Project account number

echo "Job started on $(hostname) at $(date)"

source /g100/home/userexternal/acolombi/miniconda3/etc/profile.d/conda.sh
conda activate r_env

Rscript Rscripts/CancerData/Analysis_parLapply.R

echo "Job finished at $(date)"