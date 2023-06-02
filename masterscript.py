import pandas as pd
import csv
import re
from itertools import product

# Test 1: find sequences start with forward primer and end with reverse primer (PASS)
# Test 2: find sequences start with any forward primer and end with any reverse primer (PASS)

# Read the primer CSV file into a DataFrame
df = pd.read_csv('DARTE-QM_primer_design.csv')

# Create lists for F_seq and R_seq
F_primers = df['F_truseq'].unique().tolist()
R_primers = df['R_truseq'].unique().tolist()

# Create empty df for storing the results
result_df = pd.DataFrame(columns=['F_truseq','R_truseq','Product'])

# Define the sequence
with open('ex.fastq','r') as file:
    sequence = file.read() # replace with your actual sequence or read from file

# Loop through each forward and reverse primer
for F_primer in F_primers:
    for R_primer in R_primers:
        pattern = f'({F_primer})(.*?)({R_primer})'
        # Search the sequence with the pattern
        matches = re.findall(pattern, sequence)
        
        for match in matches:
            # If found, append the result to the result df
            result_df = result_df.append({
                'F_truseq': match[0],
                'R_truseq': match[2],
                'Product': ''.join(match)
            }, ignore_index=True)

print(result_df)
