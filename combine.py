import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import time
from Bio.SeqRecord import SeqRecord

# How to run: python3 file.py file.fasta file.csv

def extract_products(input_file, primer_file_path):
    # Read the primer csv file into a df
    df = pd.read_csv(primer_file_path)

    # Convert columns into a set, only unique primers are included
    F_primers = set(df['F_truseq'])
    R_primers = set(df['R_truseq'])

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
                        # Add the information to the DataFrame
                        results_df = results_df.append({'ID': record.id, 'Product': product, 'Start position': f_pos, 'End position': r_pos + len(r_primer), 'Length': len(product), 'Start Primer': f_primer, 'End Primer': r_primer, 'From reverse complement': 'No'}, ignore_index=True)


            # Repeat the process with reverse complement primers
            for f_rev_comp in F_reverse_complement:
                for r_rev_comp in R_reverse_complement:
                    r_pos = sequence.find(r_rev_comp)
                    f_pos = sequence.find(f_rev_comp)

                    if f_pos > r_pos and f_pos != -1 and r_pos != -1:
                        product = str(Seq(sequence[r_pos:f_pos+len(f_rev_comp)]).reverse_complement())
                        # Add the information to the DataFrame
                        results_df = results_df.append({'ID': record.id, 'Product': product, 'Start position': r_pos, 'End position': f_pos + len(f_rev_comp), 'Length': len(product), 'Start Primer': r_rev_comp, 'End Primer': f_rev_comp, 'From reverse complement': 'Yes'}, ignore_index=True)


    # End the timer
    end_time = time.time()


    # Print the time it took to execute the program
    print("Execution Time: ", end_time - start_time, "seconds")

    return df, results_df


def match_products(primer_file, result_file):
    # Create empty lists to hold SeqRecords for the new fasta file
    # and matching and non-matching rows for new CSV files
    records = []
    matching_rows = []
    nonmatching_rows = []

    # Iterate over the rows in result_file
    for _, row2 in result_file.iterrows():
        # Check if the product in result_file matches a product in primer_file
        match = primer_file[primer_file['product'] == row2['Product']]
    
        # If a match is found
        if not match.empty:
            # Add target_gene and target_locus to result_file
            row2['target_gene'] = match['target_gene'].values[0]
            row2['target_locus'] = match['target_locus'].values[0]
        
            # Append the row to matching_rows list
            matching_rows.append(pd.DataFrame([row2]))
        else:
            # If no match is found, append this product to the list for the new fasta file
            records.append(SeqRecord(Seq(row2['Product']), id=row2['ID']))
        
            # Append the row to nonmatching_rows list
            nonmatching_rows.append(pd.DataFrame([row2]))
        
    # Concatenate all dataframes in the matching and non-matching rows lists if the lists are not empty
    matching_df = pd.concat(matching_rows) if matching_rows else None
    nonmatching_df = pd.concat(nonmatching_rows) if nonmatching_rows else None

    # Append to the existing CSV files if the dataframes are not None
    if matching_df is not None:
        with open('foundGene.csv', 'a') as f:
            matching_df.to_csv(f, header=f.tell()==0, index=False)
    
    if nonmatching_df is not None:
        with open('nonmatchingGene.csv', 'a') as f:
            nonmatching_df.to_csv(f, header=f.tell()==0, index=False)

    # Append to the existing fasta file
    if records:
        with open('nonmatchingGene.fasta', 'a') as f:
            SeqIO.write(records, f, "fasta")

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Parse input and output file names")
    parser.add_argument('input_file', type=str, help='Input fasta file name')

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    primer_file_path = 'DARTE-QM_primer_design.csv'

    primer_df, product_df = extract_products(input_file, primer_file_path)
    match_products(primer_df, product_df)
