#This script maps all geNomad-identified MGEs to CARD db, using RGI.

import os
import subprocess
import pandas as pd
import csv
import logging

def run_rgi_on_genomad():
    # Load the genomad_summary.csv file
    input_csv = "/cephyr/NOBACKUP/groups/jbp/matev/OceanARM/BlendARGs/result_BlendARGs/VirusIdentification/geNomad/genomad_summary.csv"
    df = pd.read_csv(input_csv)
    
    # Ensure the output directory exists
    output_dir = "RGI_genomad"
    os.makedirs(output_dir, exist_ok=True)
    
    for index, row in df.iterrows():
        gene_id = row['gene']
        protein_sequence = row['sequence']
        source = row['Source']
        sample = row['Sample']
        
        # Define file paths
        fasta_file = os.path.join(output_dir, f"{sample}_{source}_{gene_id}.faa")
        output_file = os.path.join(output_dir, f"{sample}_{source}_{gene_id}_prediction")
        
        # Write the sequence to a FASTA file
        with open(fasta_file, 'w') as f:
            f.write(f">{sample}_{source}_{gene_id}\n{protein_sequence}\n")
        
        # Run RGI
        try:
            command = f"rgi main -i {fasta_file} -o {output_file} -t protein -n 32 --clean --local --include_loose"
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {gene_id}: {e}")
    
    # Collect results into a summary file
    merge_rgi_results(output_dir)

def merge_rgi_results(output_dir):
    dfs = []
    for file in os.listdir(output_dir):
        if file.endswith('_prediction.txt'):
            file_path = os.path.join(output_dir, file)
            csv_file = file_path.replace('.txt', '.csv')
            convert_txt_to_csv(file_path, csv_file)
            df = pd.read_csv(csv_file)
            df['Sample'] = file.split('_')[0]  # Extract Sample from filename
            df['Source'] = file.split('_')[1]  # Extract Source from filename
            df['Gene'] = '_'.join(file.split('_')[2:-1])  # Extract Gene ID
            dfs.append(df)
    
    if dfs:
        merged_df = pd.concat(dfs, ignore_index=True)
        merged_df.to_csv(os.path.join(output_dir, 'summary_RGI_genomad.csv'), index=False)
        print(f"Summary CSV file saved to {os.path.join(output_dir, 'summary_RGI_genomad.csv')}")

def convert_txt_to_csv(txt_file, csv_file):
    try:
        with open(txt_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter=',')
            for row in reader:
                writer.writerow(row)
    except Exception as e:
        logging.error(f"Error converting {txt_file} to CSV: {e}")

if __name__ == "__main__":
    run_rgi_on_genomad()
