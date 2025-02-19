#This script maps all geNomad-identified MGEs to CARD db, using RGI.
import os
import subprocess
import pandas as pd
import csv
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define paths
query_csv = "genomad_summary.csv"
output_fasta = "Genomad.faa"
rgi_output = "RGI_genomad.csv"

# Read the input CSV
logging.info("Reading input CSV file.")
df = pd.read_csv(query_csv)

# Write to a single fasta file
logging.info("Creating combined FASTA file.")
with open(output_fasta, "w") as fasta_file:
    for index, row in df.iterrows():
        header = f">{row['gene']}|{row['Source']}|{row['Sample']}"
        sequence = row['sequence'].replace("*", "")  # Remove '*' if present
        fasta_file.write(f"{header}\n{sequence}\n")

# Run RGI on the combined fasta file
try:
    command = f"rgi main -i {output_fasta} -o {rgi_output} -t protein -n 32 --clean --local --include_loose"
    logging.info("Running RGI command.")
    subprocess.run(command, shell=True, check=True)
    logging.info("RGI analysis completed.")
except subprocess.CalledProcessError as e:
    logging.error(f"Error running RGI: {e}")

print(f"RGI output saved to {rgi_output}")
