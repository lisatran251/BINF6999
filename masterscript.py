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

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import csv

# Create the parser
parser = argparse.ArgumentParser(description="Parse input and output file names")
parser.add_argument('input_file', type=str, help='Input file name')

# Parse the arguments
args = parser.parse_args()
input_file = args.input_file

# Read the primer csv file into a df
df = pd.read_csv('DARTE-QM_primer_design.csv')

# Add a new column to the df with the reverse complement of the product
df['product_reverse_complement'] = df['product'].apply(lambda x: str(Seq(x).reverse_complement()))

# Create dictionaries for primers in regular and reverse complement primers
primer_dict = {}
reverse_primer_dict = {}

# Loop through df to find product, length, target gene and gene locus 
for _, row in df.iterrows():
    key = (row['product'], int(row['Product_length']))
    value = {'target_gene': row['target_gene'], 'target_locus': row['target_locus']}
    primer_dict[key] = value

    #Creating keys and values for the reverse dictionary 
    rev_key = (row['product_reverse_complement'], int(row['Product_length']))
    reverse_primer_dict[rev_key] = value

# Create lists for F_seq and R_seq
F_primers = set(df['F_truseq'].unique().tolist())
R_primers = set(df['R_truseq'].unique().tolist())

# Create reverse complement for each primer 
F_reverse_complement = {str(Seq(primer).reverse_complement()) for primer in F_primers}
R_reverse_complement = {str(Seq(primer).reverse_complement()) for primer in R_primers}

# Merge 4 lists of primers 
all_primers = F_primers | R_primers | F_reverse_complement | R_reverse_complement

# Read the sequence
with open(input_file,'r') as handle:
    records = SeqIO.parse(handle, 'fasta')
    with open('output3.csv', 'a', newline='') as result_file:
        csv_writer = csv.writer(result_file)
        #csv_writer.writerow(['Product', 'Start position', 'End position', 'Length', 'Start Primer', 'End Primer', 'From reverse complement', 'Target Gene'])

        # List to store non-matching products
        non_matching_products = []

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
                        # Get the target gene name from the primer dictionary, or 'N/A' if not found
                        target_gene = primer_dict.get((product, len(product)), {}).get('target_gene', 'N/A')
                        # Write the information to the csv file
                        csv_writer.writerow([product, f_pos, r_pos + len(r_primer), len(product), f_primer, r_primer, 'No', target_gene])

                        # If the target gene is 'N/A' and product length is between 150 and 400 base pairs, add the product to the non_matching_products list
                        if target_gene == 'N/A' and 150 <= len(product) <= 400:
                            non_matching_products.append(SeqRecord(Seq(product), id=f'{record.id}_{f_pos}_{r_pos}', description='Non-matching product'))

            # Repeat the process with reverse complement primers
            for f_rev_comp in F_reverse_complement:
                for r_rev_comp in R_reverse_complement:
                    r_pos = sequence.find(r_rev_comp)
                    f_pos = sequence.find(f_rev_comp)

                    if r_pos < f_pos and r_pos != -1 and f_pos != -1:
                        product = sequence[r_pos:f_pos+len(f_rev_comp)]
                        target_gene = reverse_primer_dict.get((product, len(product)), {}).get('target_gene', 'N/A')
                        csv_writer.writerow([product, r_pos, f_pos + len(f_rev_comp), len(product), r_rev_comp, f_rev_comp, 'Yes', target_gene])

                        # If the target gene is 'N/A' and product length is between 150 and 400 base pairs, add the product to the non_matching_products list
                        if target_gene == 'N/A' and 150 <= len(product) <= 400:
                            non_matching_products.append(SeqRecord(Seq(product), id=f'{record.id}_{r_pos}_{f_pos}', description='Non-matching product'))

# Write non-matching products to a FASTA file
with open('non_matching_products.fasta', 'a') as output_handle:
    SeqIO.write(non_matching_products, output_handle, 'fasta')

# bwa index target_seq.fa

# bwa mem target_seq.fa non_matching_products.fasta > aligned_reads.sam

# samtools view -bS aligned_reads.sam > aligned_reads.bam

# samtools sort aligned_reads.bam -o sorted_reads.bam

# samtools index sorted_reads.bam



# gatk RealignerTargetCreator -R target_seq.fa -I sorted_reads.bam -o target_intervals.list

# gatk IndelRealigner -R reference.fasta -I sorted_reads.bam -targetIntervals target_intervals.list -o realigned_reads.bam

# gatk HaplotypeCaller -R reference.fasta -I realigned_reads.bam -O raw_variants.vcf

# gatk VariantRecalibrator -R reference.fasta -V raw_variants.vcf -O recal_data.table

# gatk ApplyVQSR -R reference.fasta -V raw_variants.vcf --recal-file recal_data.table -O filtered_variants.vcf


