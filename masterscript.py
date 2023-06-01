import pandas as pd
import csv
import re

# Read the primer CSV file into a DataFrame
df = pd.read_csv('DARTE-QM_primer_design.csv')

# Create new df for F_seq and R_seq
prim_list = df[['F_truseq','R_truseq']]

# Create empty df for storing the results
result_df = pd.DataFrame(columns=['F_truseq','R_truseq','Product'])

# Define the sequence
with open('filtered_contigs_ex.fastq','r') as file:
	sequence = file.read()# replace with your actual sequence or read from file

# Loop through each row in the primer df
for _, row in prim_list.iterrows():
    F_primer = row['F_truseq']
    R_primer = row['R_truseq']
    
    pattern = f'{F_primer}(.*?){R_primer}'
    # Search the sequence with the pattern
    match = re.search(pattern, sequence)

    if match:
        # If found, append the result to the result df
        result_df = result_df.append({
            'F_truseq': F_primer,
            'R_truseq': R_primer,
            'Product': match.group(0)
        }, ignore_index=True)    

print(result_df)
