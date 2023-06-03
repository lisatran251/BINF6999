#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-08:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8000M
#SBATCH --job-name=trial2_job
#SBATCH --output=trial2_job.out

# Run the script 
python3 -m venv ~/my_venv
source ~/my_venv/bin/activate
pip install pandas numpy biopython 
python3 sample.py 
