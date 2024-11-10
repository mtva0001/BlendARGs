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
target_dir = 'DeepNOG_allbins'

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

# Counters
total_wd_folders = 0
count_hgt_folders = 0

# Walk through the directory tree starting from the source directory
for root, dirs, files in os.walk(source_dir):
    if root.endswith('_wd'):
        total_wd_folders += 1  # Increment the total _wd folders counter
        has_hgt_subfolder = False

        for subdir in dirs:
            if '_HGT_' in subdir:
                # Search for any file that matches *_detected_HGTs.txt
                hgt_file_pattern = os.path.join(root, subdir, '*_detected_HGTs.txt')
                hgt_files = glob.glob(hgt_file_pattern)  # Get list of matching files

                # Check if any matching file exists and has at least two lines
                for hgt_file in hgt_files:
                    with open(hgt_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 1:  # At least one header
                            has_hgt_subfolder = True
                            src_path = os.path.join(root, subdir)
                            dst_path = os.path.join(target_dir, subdir)

                            try:
                                shutil.copytree(src_path, dst_path)
                                logging.info(f"Copied {src_path} to {dst_path}")
                            except Exception as e:
                                logging.error(f"Failed to copy {src_path} to {dst_path}: {e}")
                            break  # Break after finding the first valid HGT file

        # Increment the count only if a valid HGT folder with matching *_detected_HGTs.txt was found
        if has_hgt_subfolder:
            count_hgt_folders += 1

logging.info(f"\nProcessing {count_hgt_folders} number of samples.")



def process_folder(folder_path):
    folder_name = os.path.basename(folder_path)
    folder_name_parts = folder_name.split('_MetaCHIP_wd')
    folder_name_sampleID = folder_name_parts[0]
    logging.info(f"Processing sample: {folder_name_sampleID}")

    # Define the file suffixes
    combined_fasta_suffix = "_combined_faa.fasta"

    for file_name in os.listdir(folder_path):
        if file_name.endswith(combined_fasta_suffix):
            file_path = os.path.join(folder_path, file_name)
            logging.info(f"Processing file: {file_path}")
            
            try:
                output_file_name = file_name.replace(combined_fasta_suffix, '_prediction.csv')
                output_file_path = os.path.join(folder_path, output_file_name)
                command = f"deepnog infer {file_path} -db eggNOG5 -d auto --out {output_file_path} -c 0.99 --verbose 3"
                subprocess.run(command, shell=True, check=True)

            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing donor file {file_path}: {e}")



def main():
    main_folder_path = "DeepNOG_allbins"

    for folder_name in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path)

if __name__ == "__main__":
    main()



# Read all CSV files into a single DataFrame
folder_path = 'DeepNOG_allbins'

dfs = []
# Walk through the directory tree to find all CSV files
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith('_prediction.csv'):
            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path, header=0, names=['ID', 'COG_ID', 'Value'])
            dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
merged_df = pd.concat(dfs, ignore_index=True)

output_file_path = os.path.join(folder_path, 'summary_DeepNOG_allbins.csv')
merged_df.to_csv(output_file_path, index=False)

print(f"Summary CSV file saved to {output_file_path}")


csv_file = "DeepNOG_allbins/summary_DeepNOG_allbins.csv"
hgt_files_folder = "DeepNOG_allbins"
annotation_file = "annotation.csv"


if csv_file:
    updated_df = pd.read_csv(csv_file)  # Load the CSV as a DataFrame
    annotation_df = pd.read_csv(annotation_file)  # Load the annotation file as a DataFrame

    # Merge the updated CSV with the annotation file based on 'COG_ID' and 'COG'
    merged_df = pd.merge(updated_df, annotation_df, left_on='COG_ID', right_on='COG', how='left')
    merged_df = merged_df.drop(merged_df.columns[[3,5]], axis=1)

    # Save the merged DataFrame to a new CSV file
    merged_df.to_csv('DeepNOG_allbins/summary_DeepNOG_allbins_COG.csv', index=False)
else:
    print("Error: The updated CSV file path is None.")



print(f'The pipeline has completed, go home and sleep!')
