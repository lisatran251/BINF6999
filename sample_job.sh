#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-24:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=5
#SBATCH --mem=3G
#SBATCH --job-name=trial4
#SBATCH --output=trial4.out

# Run the script 
python3 -m venv ~/my_venv
source ~/my_venv/bin/activate
pip install pandas numpy biopython argparse
python3 script3.py $1
