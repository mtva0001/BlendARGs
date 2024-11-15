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
target_dir = 'BLAST_allbins_latent'

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)


# Walk through the directory tree starting from the source directory
for root, dirs, files in os.walk(source_dir):
    if root.endswith('_wd'):
        for subdir in dirs:
            if '_HGT_' in subdir:
                # Search for any file that matches *_detected_HGTs.txt
                hgt_file_pattern = os.path.join(root, subdir, '*_detected_HGTs.txt')
                hgt_files = glob.glob(hgt_file_pattern)  # Get list of matching files

                # Check if any matching file exists and has at least two lines
                for hgt_file in hgt_files:
                    with open(hgt_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 1:
                            src_path = os.path.join(root, subdir)
                            dst_path = os.path.join(target_dir, subdir)

                            try:
                                shutil.copytree(src_path, dst_path)
                                logging.info(f"Copied {src_path} to {dst_path}")
                            except Exception as e:
                                logging.error(f"Failed to copy {src_path} to {dst_path}: {e}")
                            break 


def process_folder(folder_path):
    folder_name = os.path.basename(folder_path)
    folder_name_parts = folder_name.split('_MetaCHIP_wd')
    folder_name_sampleID = folder_name_parts[0]
    logging.info(f"Processing sample: {folder_name_sampleID}")

    # Define the file suffix
    combined_fasta_suffix = "_combined_faa.fasta"

    for file_name in os.listdir(folder_path):
        if file_name.endswith(combined_fasta_suffix):
            file_path = os.path.join(folder_path, file_name)
            logging.info(f"Processing file: {file_path}")
            
            try:
                output_file_name = file_name.replace(combined_fasta_suffix, '_prediction.txt')
                output_file_path = os.path.join(folder_path, output_file_name)
                command = f'blastp -query {file_path} -out {output_file_path} -num_threads 32 -db latent_ARGs_db -outfmt "6 qseqid sseqid qcovs pident length bitscore" -max_target_seqs 10'
                subprocess.run(command, shell=True, check=True)

                # Convert only the _prediction.txt files
                txt_output_file = output_file_path
                if txt_output_file.endswith('_prediction.txt') and os.path.exists(txt_output_file):
                    convert_txt_to_csv(txt_output_file, output_file_path + ".csv")

            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing file {file_path}: {e}")

        
def convert_txt_to_csv(txt_file, csv_file):
    try:
        with open(txt_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter=',')
            writer = csv.writer(outfile, delimiter=',')
            for row in reader:
                writer.writerow(row)
        logging.info(f"Converted {txt_file} to {csv_file}")

    except Exception as e:
        logging.error(f"Error converting {txt_file} to CSV: {e}")

def filter_top_hits_per_file(csv_file):
    df = pd.read_csv(csv_file, header=None, sep='\t')
    df.columns = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']

    # Sort by query and bitscore (descending)
    df_sorted = df.sort_values(by=['query', 'bitscore'], ascending=[True, False])

    # Drop duplicates based on 'query', keeping the row with the highest bitscore
    df_filtered = df_sorted.drop_duplicates(subset='query', keep='first')

    return df_filtered

def main():
    main_folder_path = "BLAST_allbins_latent"

    dfs = []

    # Process each folder
    for folder_name in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path)

    # Walk through the directory tree to find all prediction CSV files
    for root, dirs, files in os.walk(main_folder_path):
        for file in files:
            if file.endswith('_prediction.txt.csv'):
                file_path = os.path.join(root, file)
                if os.path.getsize(file_path) > 0:
                    # Filter each file to retain only the top BLAST hit based on bitscore
                    df_filtered = filter_top_hits_per_file(file_path)
                    dfs.append(df_filtered)
                else:
                    print(f"Skipping empty file: {file_path}")

    # Concatenate all DataFrames into a single DataFrame
    if dfs:
        merged_df = pd.concat(dfs, ignore_index=True)
        new_header = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']
        merged_df.columns = new_header

        # Apply additional filters for identity and coverage
        filtered_df = merged_df[(merged_df['identity%'] >= 90) & (merged_df['coverage%'] >= 20)]

        output_file_path = os.path.join(main_folder_path, 'summary_BLAST_allbins_latent.csv')
        filtered_df.to_csv(output_file_path, index=False)
        print(f"Filtered summary CSV file saved to {output_file_path}")

if __name__ == "__main__":
    main()

print(f'The pipeline has completed, go home and sleep!')
