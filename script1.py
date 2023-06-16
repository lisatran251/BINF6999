# python -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython argparse

# Check the status of the jobs 
# squeue -u thuyduye or squeue -j job_id

# Test 1: find sequences start with forward primer and end with reverse primer (PASS)
# Test 2: find sequences start with any forward primer and end with any reverse primer (PASS)
# Test 3: find sequences start with any reverse complementary and end with any of reverse primer or its reverse complementary (PASS)
# Test 4: find the start and end position, and the length of the product (PASS)
# Test 5: find target gene that matched length and product (PASS)
# Test 6 (del): count occurence of ea. target gene (PASS)
# Test 7: Overlaps 

# brew install brewsci/bio/abricate
# abricate --check
# abricate --list

# grep -c '^>' ex.fasta

#rm *chunk chunk* *chunk* info.txt *out non*

#How to run: python3 sample.py ex.fastq

#How to check runtime: sacct -j 7157312 --format=JobID,Elapsed

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import time

# Create the parser
parser = argparse.ArgumentParser(description="Parse input and output file names")
parser.add_argument('input_file', type=str, help='Input file name')

# Parse the arguments
args = parser.parse_args()
input_file = args.input_file

# Read the primer csv file into a df
df = pd.read_csv('DARTE-QM_primer_design.csv')

# Create lists for F_seq and R_seq
F_primers = set(df['F_truseq'].tolist())
R_primers = set(df['R_truseq'].tolist())

# Create reverse complement for each primer 
F_reverse_complement = {str(Seq(primer).reverse_complement()) for primer in F_primers}
R_reverse_complement = {str(Seq(primer).reverse_complement()) for primer in R_primers}

# Create an empty DataFrame to store the results
results_df = pd.DataFrame(columns=['Product', 'Start position', 'End position', 'Length', 'Start Primer', 'End Primer', 'From reverse complement'])

# Start the timer
start_time = time.time()

# Read the sequence
with open(input_file,'r') as handle:
    records = SeqIO.parse(handle, 'fasta')

    # Iterate over each sequence record in the input file
    for record in records:
        # Convert the sequence record to a string
        sequence = str(record.seq)

        # Iterate over each combination of forward and reverse primers
        for f_primer in F_primers:
            for r_primer in R_primers:
                # Find the positions of the primers in the sequence
                f_pos = sequence.find(f_primer)
                r_pos = sequence.find(r_primer)

                # If the forward primer comes before the reverse primer and both primers are found
                if f_pos < r_pos and f_pos != -1 and r_pos != -1:
                    # Extract the product sequence
                    product = sequence[f_pos:r_pos+len(r_primer)]
                    # Add the information to the DataFrame
                    results_df = results_df.append({'Product': product, 'Start position': f_pos, 'End position': r_pos + len(r_primer), 'Length': len(product), 'Start Primer': f_primer, 'End Primer': r_primer, 'From reverse complement': 'No'}, ignore_index=True)

        # Repeat the process with reverse complement primers
        for f_rev_comp in F_reverse_complement:
            for r_rev_comp in R_reverse_complement:
                r_pos = sequence.find(r_rev_comp)
                f_pos = sequence.find(f_rev_comp)

                if r_pos < f_pos and r_pos != -1 and f_pos != -1:
                    product = sequence[r_pos:f_pos+len(f_rev_comp)]
                    results_df = results_df.append({'Product': product, 'Start position': r_pos, 'End position': f_pos + len(f_rev_comp), 'Length': len(product), 'Start Primer': r_rev_comp, 'End Primer': f_rev_comp, 'From reverse complement': 'Yes'}, ignore_index=True)

# End the timer
end_time = time.time()

# Write the DataFrame to a csv file
results_df.to_csv('test_wo_df.csv', index=False)

# Print the time it took to execute the program
print("Execution Time: ", end_time - start_time, "seconds")
