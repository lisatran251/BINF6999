# python -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython argparse

# Check the status of the job 
# squeue -u thuyduye or squeue -j job_id

# Test 1: find sequences start with forward primer and end with reverse primer (PASS)
# Test 2: find sequences start with any forward primer and end with any reverse primer (PASS)
# Test 3: find sequences start with any forward primer or its reverse complementary and end with any of reverse primer or its reverse complementary (PASS)

#How to run: python3 sample.py ex.fastq

import pandas as pd
import re
from Bio.Seq import Seq
from Bio import SeqI
import argparse

# Create the parser
parser = argparse.ArgumentParser(description="Parse input and output file names")
parser.add_argument('input_file', type=str, help='Input file name')

# Parse the arguments
args = parser.parse_args()
input_file = args.input_file

# Read the primer CSV file into a DataFrame
df = pd.read_csv('DARTE-QM_primer_design.csv')
primer_dict = {}
for _, row in df.iterrows():
    key = (row['product'], int(row['Product_length']))
    value = {'target_gene': row['target_gene'], 'target_locus': row['target_locus']}
    primer_dict[key] = value

# Create lists for F_seq and R_seq
F_primers = df['F_truseq'].unique().tolist()
R_primers = df['R_truseq'].unique().tolist()

# Create reverse complement for each primer 
F_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in F_primers]
R_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in R_primers]

# Precompile regular expressions
patterns = [(re.compile(f'({F_primer})(.*?)({R_primer})'), F_primer, R_primer) 
            for F_primer in F_primers for R_primer in R_primers]

rev_patterns = [(re.compile(f'({F_rev_comp})(.*?)({R_rev_comp})'), F_rev_comp, R_rev_comp) 
                for F_rev_comp in F_reverse_complement for R_rev_comp in R_reverse_complement]

# Store results in lists first
results = []
rev_results = []

# Read the sequence
with open(input_file,'r') as handle:
    records = SeqIO.parse(handle, 'fasta')

    for record in records:
        sequence = str(record.seq)
        
        # Loop through each pattern
        for pattern, F_primer, R_primer in patterns:
            matches = pattern.finditer(sequence)
            for match in matches:
                # Append the results to the list
                results.append({
                    'F_primer': F_primer,
                    'R_primer': R_primer,
                    'Product': match.group(),
                    'Start': match.start(),
                    'End': match.end(),
                    'Length': match.end() - match.start()
                })

        for rev_pattern, F_rev_comp, R_rev_comp in rev_patterns:
            rev_matches = rev_pattern.finditer(sequence)
            for match in rev_matches:
                rev_results.append({
                    'F_reverse_complement': F_rev_comp,
                    'R_reverse_complement': R_rev_comp,
                    'Product': match.group(),
                    'Start': match.start(),
                    'End': match.end(),
                    'Length': match.end() - match.start()
                })

# Convert results to DataFrame after the loop
result_df = pd.DataFrame(results)
result_reverse_df = pd.DataFrame(rev_results)

print(result_df)
print(result_reverse_df)
# result_df.to_csv('result.csv', sep='\t', index=False)
# result_reverse_df.to_csv('result_reverse.csv', sep='\t', index=False)

# # Iterate through result_df and result_reverse_df to check the matching keys in primer_dict
# for _, row in result_df.iterrows():
#     key = (row['Product'], row['Length'])
#     if key in primer_dict:
#         print("Product: ", row['Product'])
#         print("Target Gene: ", primer_dict[key]['target_gene'])
#         print("Target Locus: ", primer_dict[key]['target_locus'])
#         print("\n")

# for _, row in result_reverse_df.iterrows():
#     key = (row['Product'], row['Length'])
#     if key in primer_dict:
#         print("Product: ", row['Product'])
#         print("Target Gene: ", primer_dict[key]['target_gene'])
#         print("Target Locus: ", primer_dict[key]['target_locus'])
#         print("\n")
