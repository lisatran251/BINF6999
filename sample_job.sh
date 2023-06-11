#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-168:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=5
#SBATCH --mem=8G
#SBATCH --job-name=fulljob3
#SBATCH --output=fulljob3.out

# Run the script 
python3 -m venv ~/my_venv
source ~/my_venv/bin/activate
pip install pandas numpy biopython argparse
python3 updatedscript.py contigs_ex.fasta
