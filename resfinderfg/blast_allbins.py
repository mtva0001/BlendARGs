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
target_dir = 'BLAST_allbins'

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

# Counters
total_wd_folders = 0
count_hgt_folders = 0

# Walk through the directory tree starting from the source directory
for root, dirs, files in os.walk(source_dir):
    if root.endswith('_wd'):
        total_wd_folders += 1
        has_hgt_subfolder = False

        for subdir in dirs:
            if '_HGT_' in subdir:
                # Search for any file that matches *_detected_HGTs.txt
                hgt_file_pattern = os.path.join(root, subdir, '*_detected_HGTs.txt')
                hgt_files = glob.glob(hgt_file_pattern)

                # Check if any matching file exists and has at least two lines
                for hgt_file in hgt_files:
                    with open(hgt_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 1:
                            has_hgt_subfolder = True
                            src_path = os.path.join(root, subdir)
                            dst_path = os.path.join(target_dir, subdir)

                            try:
                                shutil.copytree(src_path, dst_path)
                                logging.info(f"Copied {src_path} to {dst_path}")
                            except Exception as e:
                                logging.error(f"Failed to copy {src_path} to {dst_path}: {e}")
                            break

        if has_hgt_subfolder:
            count_hgt_folders += 1

logging.info(f"\nProcessing {count_hgt_folders} number of samples.")


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
            logging.info(f"Processing donor file: {file_path}")

            try:
                output_file_name = file_name.replace(combined_fasta_suffix, '_prediction.txt')
                output_file_path = os.path.join(folder_path, output_file_name)
                command = f'blastp -query {file_path} -out {output_file_path} -num_threads 32 -db resfinderfg2.0 -outfmt "6 qseqid sseqid qcovs pident length bitscore" -max_target_seqs 10'
                subprocess.run(command, shell=True, check=True)

                # Convert only the _prediction.txt files
                txt_output_file = output_file_path
                if txt_output_file.endswith('_prediction.txt') and os.path.exists(txt_output_file):
                    convert_txt_to_csv(txt_output_file, output_file_path + ".csv")

            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing donor file {file_path}: {e}")


def convert_txt_to_csv(txt_file, csv_file):
    """
    Converts a comma-separated .txt file to a comma-separated .csv file.
    Only processes files ending with '_prediction.txt'.
    """
    try:
        with open(txt_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter=',')
            writer = csv.writer(outfile, delimiter=',')
            for row in reader:
                writer.writerow(row)
        logging.info(f"Converted {txt_file} to {csv_file}")

    except Exception as e:
        logging.error(f"Error converting {txt_file} to CSV: {e}")


def main():
    main_folder_path = "BLAST_allbins"

    for folder_name in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path)

if __name__ == "__main__":
    main()


# Read all CSV files into a single DataFrame
folder_path = 'BLAST_allbins'
dfs = []

# Walk through the directory tree to find all CSV files
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith('_prediction.txt.csv'):
            file_path = os.path.join(root, file)
            if os.path.getsize(file_path) > 0:
                df = pd.read_csv(file_path, header=None)
                dfs.append(df)
            else:
                print(f"Skipping empty file: {file_path}")

# Concatenate all DataFrames into a single DataFrame
merged_df = pd.concat(dfs, ignore_index=True)
new_header = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']
merged_df.columns = new_header

# Apply filters
filtered_df = merged_df[(merged_df['identity%'] >= 90) & (merged_df['coverage%'] >= 20)]
output_file_path = os.path.join(folder_path, 'summary_BLAST_allbins_filtered.csv')
filtered_df.to_csv(output_file_path, index=False)
print(f"Filtered summary CSV file saved to {output_file_path}")


csv_file = 'BLAST_allbins/summary_BLAST_allbins_filtered.csv'
annotation_file = 'annotation_resfinder.csv'

if csv_file:
    updated_df = pd.read_csv(csv_file)
    annotation_df = pd.read_csv(annotation_file)

    # Merge the updated CSV with the annotation file based on 'target' and 'ID'
    merged_df = pd.merge(updated_df, annotation_df, left_on='target', right_on='ID', how='left')
    merged_df.to_csv('BLAST_allbins/summary_BLAST_allbins_RF.csv', index=False)
else:
    print("Error: The updated CSV file path is None.")

print(f'The pipeline has completed, go home and sleep!')
