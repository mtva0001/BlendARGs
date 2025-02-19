import os
import glob
import shutil
import subprocess
import logging
import pandas as pd
import csv

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the source directory and target directory
source_dir = '.'
target_dir = '.'
db = "resfinderfg2.0"

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

# Locate the single .faa file generated earlier
faa_file = "Genomad.faa"  # Update if necessary
if not os.path.exists(faa_file):
    logging.error(f"Query file {faa_file} not found!")
    exit(1)

# Define output BLAST file
blast_output = os.path.join(target_dir, "RF_blast_genomad.txt")

# Run BLASTP on the .faa file
try:
    blast_command = (
        f'blastp -query {faa_file} -out {blast_output} -num_threads 12 '
        f'-db {db} -outfmt "6 qseqid sseqid qcovs pident length bitscore" -max_target_seqs 10'
    )
    subprocess.run(blast_command, shell=True, check=True)
    logging.info(f"BLASTP completed successfully, output saved to {blast_output}")
except subprocess.CalledProcessError as e:
    logging.error(f"Error running BLASTP: {e}")
    exit(1)

# Convert BLAST output to CSV
def convert_txt_to_csv(txt_file, csv_file):
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None)
        df.columns = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']
        df.to_csv(csv_file, index=False)
        logging.info(f"Converted {txt_file} to {csv_file}")
    except Exception as e:
        logging.error(f"Error converting {txt_file} to CSV: {e}")

csv_output = blast_output + ".csv"
convert_txt_to_csv(blast_output, csv_output)

# Filter the BLAST results
def filter_top_hits(csv_file):
    df = pd.read_csv(csv_file)
    df_sorted = df.sort_values(by=['query', 'bitscore'], ascending=[True, False])
    df_filtered = df_sorted.drop_duplicates(subset='query', keep='first')
    return df_filtered

filtered_df = filter_top_hits(csv_output)

# Apply identity and coverage filters
filtered_df = filtered_df[(filtered_df['identity%'] >= 90) & (filtered_df['coverage%'] >= 20)]

# Save final summary
final_summary = os.path.join(target_dir, "ResFinderFG_genomad.csv")
filtered_df.to_csv(final_summary, index=False)
logging.info(f"Filtered summary saved to {final_summary}")
