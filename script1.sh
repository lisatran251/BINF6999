#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-168:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=5
#SBATCH --mem=6G
#SBATCH --job-name=trial1
#SBATCH --output=trial1.out

# Run the script 
python3 -m venv ~/my_venv
source ~/my_venv/bin/activate
pip install pandas numpy biopython argparse
python3 script1.py $1
