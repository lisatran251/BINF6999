#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-168:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=5
#SBATCH --mem=6G
#SBATCH --job-name=trial1
#SBATCH --output=trial1.out

# Input fasta file
fasta_file="contigs_ex.fasta"

# Prefix for the output files
prefix="chunk_"

# Number of sequences per file
seqs_per_file=10000

# Split the fasta file into chunks and submit a job for each chunk
awk -v size="$seqs_per_file" -v pre="$prefix" 'BEGIN {n_seq=0; n_file=1; file=sprintf("%s%02d.fasta", pre, n_file);} /^>/ {if(n_seq%size==0 && n_seq>0){close(file); print n_file-1, file >> "info.txt"; system("sbatch script1.sh " file); file=sprintf("%s%02d.fasta", pre, ++n_file);} n_seq++;} {print >> file; } END {print n_file, file >> "info.txt"; system("sbatch script1.sh " file);}' $fasta_file
