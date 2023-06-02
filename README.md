# Evaluating Theoretical Abundance Of Antimicrobial Resistance Gene Targets In Publicly Available Metagenomic Dataset

## 05/27 
Some terminologies for review 
Primer is a short DNA or RNA sequence that serves as starting point for DNA synthesis. It is typically designed to be complementary to a specific target region of the DNA template that you want to amplify or analyze. 

In term of PCR, which is a technique to amplify specific DNA sequences, two primers are used: a forward primer and a reverse primer. These primers are typically designed to bind to opposite ends of the target DNA sequence. 

Forward primer: complementary to sense (coding) strand of DNA template. Same sequence as RNA transcript (except T instead of U) 

Reverse primer: complementary to antisense (non-coding) strand of DNA template. Same sequence as sense strand and has opposite sequence.

During PCR, the DNA template is denatured (heated) to separate the two DNA strands. Then, the forward and reverse primers anneal or bind to their complementary sequences on the template DNA. The primers provide a starting point for DNA polymerase to begin synthesizing new DNA strands, using the template as a guide.

Primer length: 18-25 nucleotides. Short: lack specificity. Long: non-specific binding. 

Target gene: GOI
- specific GOI within an organism’s genome
- represent particular DNA seq that encode a functional product, such as protein or non-coding RNA molecule
- chosen based on relevance to biological process, disease, or research qs
- to study gene expression, genetic variation or mutations, gene function or disease mechanism 
Target locus: genomic location or region where gene located
- specific genomic location or region within the genome
- physical position on chromosome where particular gene or DNA seq located. 
- include not only coding region, but also promoters, enhancers etc.

GC content:
- provide insights into various aspects of DNA or RNA seqs 
- high GC content affect DNA stability, melting temp, and formation of secondary structures
- influence gene expression, DNA-protein interactions, primer designs for techniques like PCR

Assembly:
- taking large number of short DNA sequences, called “reads”, and putting them back together to recreate the original chromosomes from which the DNA originated. 

## 05/28
### DARTE-QM_primer_design
target_gene: This is the identifier of the gene that the primer is designed to amplify. It's often a combination of the GenBank accession number and the location of the gene in the sequence.
target_locus: Similar to target_gene, it provides a reference to the location of the target gene within a specific genome sequence.
F_truseq and R_truseq: These are the sequences of the forward (F) and reverse (R) primers respectively. These sequences will bind to complementary sequences in the DNA that flank the target sequence.
Product_length: This is the length of the DNA fragment that would be produced if the primers successfully bind to their target and a DNA polymerase extends the primer. It includes the length of the sequence between the primers as well as the primers themselves.
F_length and R_length: These are the lengths of the forward and reverse primers.
F_GC and R_GC: These values represent the GC content of the forward and reverse primers respectively. The GC content of a DNA sequence is the percentage of nucleotides in the sequence that are either guanine (G) or cytosine (C). A higher GC content in a primer can increase its melting temperature, potentially affecting the conditions needed for successful PCR.
product: This represents the expected sequence of the entire amplicon – that is, the sequence that would be produced if a PCR were conducted with these primers. It includes the sequences of the primers themselves (at the beginning and end of the sequence) and the sequence that lies between them in the DNA.
  
### How to download the data from SRA NCBI 
```
module spider sra-toolkit

prefetch SRRXXXXXX

fastq-dump --split-files SRRXXXXXX
```
Create a txt file which contain the names of the name of sequence identifier. This will help avoid mistake when we input the name of the identifier everytime and easier to keep track because there are 56 sequence identifier. 
To load the metagenomic data into Graham: using SRA toolkits which will automatically downloaded both strands of each identifier. 
Since each sequence will take apx 15-20 mins to load. I plan to write a Shell script and submit it using sbatch which will simultaneously download the data. Sometimes the program may miss a few sequences, we can do wc -l to see if all the sequences are downloaded, I can also write a small script to print out the missing data sequences, so we can try to reload them again.

## 05/29
Fastq analysis
@SRR8931189.1 1 length=150
NCCTTGGAGGGATGTTTACCCTGCAAAATCATTTTGAAAACTGCATGAGGCTTTCTATTGGCTCTTGGTCTGAAGAAGTT
GAACAAAAATTAAAACTATCAGGGAAATTGGCTGCAAAGATGTAAAACCCTAAGCGGATCCAAAAACCTC
+SRR8931189.1 1 length=150
#AAFFJJJJJJJJAJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJAJFJJJFJJJJJAJFJ
FJJJJJJJJJJJJJJJFJJJJJJJFFFJA<-FFJAFFJFFJJJJJJJJFFJJFJJFFFJFAFFFFJJAFJ

1. SRR8931189.93: identifier with 93 is the internal identifier?, length 148 bps
2. NCC…CTC: DNA sequence, N is unknown or ‘any’ nucleotide
3. +: start of same seq identifier, separate seq from quality score
4. AAF…AFJ: quality score 

There are 56 FASTQ files in total 

### Create new environment to run Python on Graham (with packages) 
*consult profs if there are other methods*
```
module load python/3.x
python -m venv ~/my_venv
source ~/my_venv/bin/activate
pip install pandas numpy
deactivate 
```

## 05/31
Contigs: contiguous sequences of DNA that are assembled from the shorter sequences generated during sequencing. 
Genome assembly: generate the fewest and longest contigs possible, ideally matching the number and size of the chromosomes in organism’s genome.
Sucessfully find sequences start with forward primer and end with reverse primer
Draft project proposal

## 06/01 
Adjust the code to only find the unique values in the both primer sets (see masterscript.py version 3) 
Completed project proposal, waiting for feedback now 
