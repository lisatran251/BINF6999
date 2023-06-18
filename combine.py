import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os 

def extract_products(input_file, primer_file):
    df = pd.read_csv(primer_file)

    F_primers = list(df['F_truseq'])
    R_primers = list(df['R_truseq'])

    F_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in F_primers]
    R_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in R_primers]

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

    matching_df.to_csv('foundGene.csv', index=False, mode='a')
    nonmatching_df.to_csv('nonmatchingGene.csv', index=False, mode='a')

    # with open('nonmatchingGene.fasta', 'a') as f:
    #     SeqIO.write(records, f, "fasta")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str)
    args = parser.parse_args()
    input_file = args.input_file
    primer_file = 'DARTE-QM_primer_design.csv'

    primer_df, product_df = extract_products(input_file, primer_file)
    match_products(primer_df, product_df)
