# python -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython

# Check the status of the job 
# squeue -u thuyduye or squeue -j job_id

# Test 1: find sequences start with forward primer and end with reverse primer (PASS)
# Test 2: find sequences start with any forward primer and end with any reverse primer (PASS)
# Test 3: find sequences start with any forward primer or its reverse complementary and end with any of reverse primer or its reverse complementary (PASS)

import pandas as pd
import re
from Bio.Seq import Seq
from Bio import SeqIO

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

# Create empty df for storing the results
result_df = pd.DataFrame(columns=['F_primer','R_primer','Product', 'Start', 'End', 'Length'])
result_reverse_df = pd.DataFrame(columns=['F_reverse_complement','R_reverse_complement','Product', 'Start', 'End', 'Length'])

# Define the sequence
with open('contigs_short2.fasta','r') as handle:
    records = SeqIO.parse(handle, 'fasta')
    for record in records:
        sequence = str(record.seq)
        
        # Loop through each forward and reverse primer
        for F_primer in F_primers:
            for R_primer in R_primers:
                pattern = f'({F_primer})(.*?)({R_primer})'
                # Search the sequence with the pattern
                matches = re.finditer(pattern, sequence)
                
                for match in matches:
                    # If found, append the result to the result df
                    result_df = result_df.append({
                        'F_primer': F_primer,
                        'R_primer': R_primer,
                        'Product': match.group(),
                        'Start': match.start(),
                        'End': match.end(),
                        'Length': match.end() - match.start()
                    }, ignore_index=True)

        for F_rev_comp in F_reverse_complement:
            for R_rev_comp in R_reverse_complement:
                rev_pattern = f'({F_rev_comp})(.*?)({R_rev_comp})' 
                # Search the sequence with the reverse pattern
                rev_matches = re.finditer(rev_pattern, sequence)
                
                for match in rev_matches:
                    # If found, append the result to the result reverse df
                    result_reverse_df = result_reverse_df.append({
                        'F_reverse_complement': F_rev_comp,
                        'R_reverse_complement': R_rev_comp,
                        'Product': match.group(),
                        'Start': match.start(),
                        'End': match.end(),
                        'Length': match.end() - match.start()
                    }, ignore_index=True)

result_df.to_csv('result.csv', sep='\t', index=False)
result_reverse_df.to_csv('result_reverse.csv', sep='\t', index=False)

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

