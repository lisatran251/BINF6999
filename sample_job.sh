#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-04:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=800M
#SBATCH --job-name=trial3_job
#SBATCH --output=trial3_job.out

# Run the script 
python3 -m venv ~/my_venv
source ~/my_venv/bin/activate
pip install pandas numpy biopython 
python3 sample.py contigs_short.fasta
