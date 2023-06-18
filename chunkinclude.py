from Bio import SeqIO
import math
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os 
import time

def chunkify_fasta(input_file, output_prefix, seqs_per_chunk):
    records = list(SeqIO.parse(input_file, 'fasta'))
    total_records = len(records)
    chunk_count = math.ceil(total_records / seqs_per_chunk)
    
    for i in range(chunk_count):
        chunk_records = records[i * seqs_per_chunk : (i+1) * seqs_per_chunk]
        SeqIO.write(chunk_records, f"{output_prefix}_{i}.fasta", "fasta")

def extract_products(input_file, primer_file):
    df = pd.read_csv(primer_file)

    F_primers = list(df['F_truseq'])
    R_primers = list(df['R_truseq'])

    F_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in F_primers]
    R_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in R_primers]

    # Start the timer
    start_time = time.time()

    results = []

    with open(input_file,'r') as handle:
        records = SeqIO.parse(handle, 'fasta')

        for record in records:
            sequence = str(record.seq)
            for f_primer, r_primer, f_rev_comp, r_rev_comp in zip(F_primers, R_primers, F_reverse_complement, R_reverse_complement):
                f_pos, r_pos, f_pos_rc, r_pos_rc = sequence.find(f_primer), sequence.find(r_primer), sequence.find(f_rev_comp), sequence.find(r_rev_comp)
                
                # Add condition for product length
                if f_pos < r_pos != -1:
                    product = sequence[f_pos:r_pos+len(r_primer)]
                    if 150 <= len(product) <= 400:  # Ensure the length of product is within 150 to 400
                        result = {'ID': record.id, 'Product': product, 'Start position': f_pos, 'End position': r_pos + len(r_primer), 'Length': len(product), 'Start Primer': f_primer, 'End Primer': r_primer, 'From reverse complement': 'No'}
                        results.append(result)
                if r_pos_rc < f_pos_rc != -1:
                    product = str(Seq(sequence[r_pos_rc:f_pos_rc+len(f_rev_comp)]).reverse_complement())
                    if 150 <= len(product) <= 400:  # Ensure the length of product is within 150 to 400
                        result = {'ID': record.id, 'Product': product, 'Start position': r_pos_rc, 'End position': f_pos_rc + len(f_rev_comp), 'Length': len(product), 'Start Primer': r_rev_comp, 'End Primer': f_rev_comp, 'From reverse complement': 'Yes'}
                        results.append(result)

    # End the timer
    end_time = time.time()

    # Print the time it took to execute the program
    print("Execution Time: ", end_time - start_time, "seconds")

    results_df = pd.DataFrame(results)
    return df, results_df

def match_products(primer_file, result_file):
    matching_rows = []
    nonmatching_rows = []
    records = []

    primer_products = set(primer_file['product'])

    for _, row in result_file.iterrows():
        product = row['Product']
        if product in primer_products:
            match = primer_file[primer_file['product'] == product]
            row['target_gene'] = match['target_gene'].values[0]
            row['target_locus'] = match['target_locus'].values[0]
            matching_rows.append(row)
        else:
            #records.append(SeqRecord(Seq(row['Product']), id=row['ID']))
            nonmatching_rows.append(row)

    matching_df = pd.DataFrame(matching_rows)
    nonmatching_df = pd.DataFrame(nonmatching_rows)

    if not matching_df.empty:
        if not os.path.isfile('foundGene.csv'):
            matching_df.to_csv('foundGene.csv', index=False)
        else: 
            matching_df.to_csv('foundGene.csv', mode='a', header=False, index=False)

    if not nonmatching_df.empty:
        if not os.path.isfile('nonmatchingGene.csv'):
            nonmatching_df.to_csv('nonmatchingGene.csv', index=False)
        else: 
            nonmatching_df.to_csv('nonmatchingGene.csv', mode='a', header=False, index=False)


    # with open('nonmatchingGene.fasta', 'a') as f:
    #     SeqIO.write(records, f, "fasta")
        
def main(input_file, primer_file):
    # Chunkify the fasta file
    chunk_prefix = "chunk"
    chunkify_fasta(input_file, chunk_prefix, 5)

    # Apply extraction and matching for each chunk
    i = 0
    while os.path.exists(f"{chunk_prefix}_{i}.fasta"):
        primer_df, product_df = extract_products(f"{chunk_prefix}_{i}.fasta", primer_file)
        match_products(primer_df, product_df)
        i += 1
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str)
    args = parser.parse_args()
    input_file = args.input_file
    primer_file = 'DARTE-QM_primer_design.csv'

    main(input_file, primer_file)
