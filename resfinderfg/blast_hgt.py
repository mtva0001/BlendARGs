import os
import shutil
import subprocess
import pandas as pd
import csv
import logging
import glob

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the source directory and the target directory
source_dir = '.'
target_dir = 'BLAST'

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

# Walk through the directory tree starting from the source directory
for root, dirs, files in os.walk(source_dir):
    if root.endswith('_wd'):
        for subdir in dirs:
            if '_HGT_' in subdir:
                hgt_file_pattern = os.path.join(root, subdir, '*_detected_HGTs.txt')
                hgt_files = glob.glob(hgt_file_pattern)

                for hgt_file in hgt_files:
                    with open(hgt_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 2:
                            src_path = os.path.join(root, subdir)
                            dst_path = os.path.join(target_dir, subdir)
                            try:
                                shutil.copytree(src_path, dst_path)
                                logging.info(f"Copied {src_path} to {dst_path}")
                            except Exception as e:
                                logging.error(f"Failed to copy {src_path} to {dst_path}: {e}")
                            break

# Function to remove '*' characters from the input .faa files
def clean_faa_file(file_path):
    cleaned_file_path = file_path + ".fasta"
    with open(file_path, 'r') as input_file, open(cleaned_file_path, 'w') as output_file:
        for line in input_file:
            if not line.startswith('>'):
                line = line.replace('*', '')
            output_file.write(line)
    return cleaned_file_path

# Function to process each donor or recipient file
def process_folder(folder_path):
    donor_suffix = "_donor_genes.faa"
    recipient_suffix = "_recipient_genes.faa"

    for file_name in os.listdir(folder_path):
        if file_name.endswith(donor_suffix) or file_name.endswith(recipient_suffix):
            file_path = os.path.join(folder_path, file_name)
            cleaned_file_path = clean_faa_file(file_path)
            suffix = "_donor_prediction.txt" if file_name.endswith(donor_suffix) else "_recipient_prediction.txt"
            output_file_path = os.path.join(folder_path, file_name.replace(donor_suffix if suffix == "_donor_prediction.txt" else recipient_suffix, suffix))
            command = f'blastp -query {cleaned_file_path} -out {output_file_path} -num_threads 32 -evalue 1e-10 -db resfinderfg2.0 -outfmt "6 qseqid sseqid qcovs pident length bitscore" -max_target_seqs 10'
            
            try:
                subprocess.run(command, shell=True, check=True)
                if os.path.exists(output_file_path):
                    # Convert to CSV and retain only the top hit for each query
                    convert_txt_to_csv(output_file_path, output_file_path + ".csv")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing file {cleaned_file_path}: {e}")

# Function to convert text output to CSV format, retaining only the top hit based on bitscore
def convert_txt_to_csv(txt_file, csv_file):
    try:
        # Read the BLAST output into a DataFrame
        columns = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']
        df = pd.read_csv(txt_file, sep='\t', names=columns, header=None)
        
        # Sort by 'query' and 'bitscore' to retain only the top hit per query
        top_hits = df.sort_values(by=['query', 'bitscore'], ascending=[True, False]).drop_duplicates(subset='query', keep='first')
        
        # Save the top hits to CSV
        top_hits.to_csv(csv_file, index=False)
        
        logging.info(f"Converted {txt_file} to {csv_file} with top hits only")
    except Exception as e:
        logging.error(f"Error converting {txt_file} to CSV: {e}")

# Main function to run processing on all folders in BLAST
def main():
    main_folder_path = "BLAST"
    for folder_name in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path)

if __name__ == "__main__":
    main()

# Concatenate all CSV files into a single DataFrame
folder_path = 'BLAST'
dfs = []
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith('_prediction.txt.csv'):
            file_path = os.path.join(root, file)
            if os.path.getsize(file_path) > 0:
                df = pd.read_csv(file_path, header=0)
                dfs.append(df)
merged_df = pd.concat(dfs, ignore_index=True)

# Filter by identity >= 95% and coverage >= 80%
filtered_blast_df = merged_df[(merged_df["identity%"] >= 95) & (merged_df["coverage%"] >= 80)]
output_file_path = os.path.join(folder_path, 'summary_BLAST.csv')
filtered_blast_df.to_csv(output_file_path, index=False)
print(f"Summary CSV file saved to {output_file_path}")

# Create a gene list
gene_list = []
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith('_detected_HGTs.txt'):
            filepath = os.path.join(root, file)
            with open(filepath, 'r') as txtfile:
                reader = csv.DictReader(txtfile, delimiter='\t')
                genes_in_file = {row['Gene_1'] for row in reader}.union({row['Gene_2'] for row in reader})
            gene_list.extend([gene for gene in genes_in_file if gene not in gene_list])

output_file_path = os.path.join('BLAST', 'gene_list.csv')
with open(output_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(['Gene'])
    writer.writerows([[gene] for gene in gene_list])
logging.info(f"Gene list saved to {output_file_path}")

# Function to determine the direction (donor/recipient) and retrieve Identity
def determine_direction_and_identity(gene_id, hgt_files_folder):
    for root, dirs, files in os.walk(hgt_files_folder):
        for filename in files:
            if filename.endswith('_detected_HGTs.txt'):
                file_path = os.path.join(root, filename)
                with open(file_path, 'r') as file:
                    for line in file:
                        columns = line.strip().split('\t')
                        direction_parts = columns[-1]
                        identity = columns[2]
                        direction_values = direction_parts.split('-->')
                        if len(direction_values) == 2:
                            donor_part = direction_values[0].strip()
                            recipient_part = direction_values[1].strip()
                            if gene_id.startswith(donor_part):
                                return 'donor', identity
                            elif gene_id.startswith(recipient_part):
                                return 'recipient', identity
    return None, None

# Function to update the CSV file with the new 'DorR' and 'Identity' columns
def update_csv(csv_file, hgt_files_folder):
    updated_rows = []
    with open(csv_file, 'r') as file:
        header = file.readline().strip()
        updated_header = f"{header},DorR,HGT_identity"
        updated_rows.append(updated_header)

        for line in file:
            gene_id = line.strip().split(',')[0].strip('"')
            direction, identity = determine_direction_and_identity(gene_id, hgt_files_folder)
            if direction and identity:
                updated_row = f"{line.strip()},{direction},{identity}"
            else:
                updated_row = f"{line.strip()},,"
            updated_rows.append(updated_row)

    output_csv_file = os.path.splitext(csv_file)[0] + '_DorR.csv'
    with open(output_csv_file, 'w') as outfile:
        for row in updated_rows:
            outfile.write(row + '\n')
    return output_csv_file

# Usage
csv_file = 'BLAST/summary_BLAST.csv'
hgt_files_folder = 'BLAST'
annotation_file = 'annotation_resfinder.csv'

updated_csv_file = update_csv(csv_file, hgt_files_folder)

if updated_csv_file:
    updated_df = pd.read_csv(updated_csv_file)
    annotation_df = pd.read_csv(annotation_file)
    merged_df = pd.merge(updated_df, annotation_df, left_on='target', right_on='ID', how='left')
    merged_df = merged_df.drop(merged_df.columns[[8]], axis=1)
    merged_df.to_csv('BLAST/summary_BLAST_DorR.csv', index=False)
else:
    print("Error: The updated CSV file path is None.")

print("The pipeline has completed, go home and sleep!")
