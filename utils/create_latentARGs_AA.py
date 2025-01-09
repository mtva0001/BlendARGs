# Table S5 from the paper https://doi.org/10.1186/s40168‑023‑01479‑0 was used to compile the Latent ARGs database.

import csv
from Bio.Seq import Seq

# Function to translate DNA sequences to protein sequences and write to FASTA
def translate_dna_to_protein(input_csv, output_fasta):
    with open(input_csv, 'r') as csvfile, open(output_fasta, 'w') as fastafile:
        reader = csv.reader(csvfile)
        next(reader)
        
        for row in reader:
            header = row[0].replace("20k_concatenated-long-orfs_", "").replace("@@@", "_")  # Header from the 1st column, way too long
            dna_sequence = row[3]  # DNA sequence from the 4th column
            
            # Translate DNA to protein and remove stop codons (*)
            protein_sequence = str(Seq(dna_sequence).translate()).replace('*', '')
            
            fastafile.write(f">{header}\n")
            fastafile.write(f"{protein_sequence}\n")


translate_dna_to_protein('Latent_ARGs.csv', 'Latent_ARGs_AA.fasta')
