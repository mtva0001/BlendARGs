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
                        if len(lines) >= 2:  # At least one header and one data row
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

logging.info(f"\nSummary: In {count_hgt_folders} out of {total_wd_folders} samples MetaCHIP have detected HGT events.")



# Function to remove '*' characters from the input .faa files
def clean_faa_file(file_path):
    cleaned_file_path = file_path + ".fasta"
    with open(file_path, 'r') as input_file, open(cleaned_file_path, 'w') as output_file:
        for line in input_file:
            if not line.startswith('>'):  # Only clean sequence lines, not headers
                line = line.replace('*', '')
            output_file.write(line)
    return cleaned_file_path


def process_folder(folder_path):
    folder_name = os.path.basename(folder_path)
    folder_name_parts = folder_name.split('_MetaCHIP_wd')
    folder_name_sampleID = folder_name_parts[0]
    logging.info(f"Processing sample: {folder_name_sampleID}")

    # Define the file suffixes
    donor_suffix = "_donor_genes.faa"
    recipient_suffix = "_recipient_genes.faa"

    for file_name in os.listdir(folder_path):
        # Check for donor files
        if file_name.endswith(donor_suffix):
            file_path = os.path.join(folder_path, file_name)
            logging.info(f"Processing donor file: {file_path}")

            # Clean the .faa file to remove '*' characters
            cleaned_file_path = clean_faa_file(file_path)
            try:
                output_file_name = file_name.replace(donor_suffix, '_donor_prediction.txt')
                output_file_path = os.path.join(folder_path, output_file_name)
                command = f'blastp -query {cleaned_file_path} -out {output_file_path} -num_threads 32 -db resfinderfg2.0 -outfmt "10 qseqid sseqid qcovs pident length" -max_target_seqs 1'
                subprocess.run(command, shell=True, check=True)

                # Convert only the _prediction.txt files
                txt_output_file = output_file_path
                if txt_output_file.endswith('_prediction.txt') and os.path.exists(txt_output_file):
                    convert_txt_to_csv(txt_output_file, output_file_path + ".csv")

            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing donor file {cleaned_file_path}: {e}")

        # Check for recipient files
        elif file_name.endswith(recipient_suffix):
            file_path = os.path.join(folder_path, file_name)
            logging.info(f"Processing recipient file: {file_path}")

            # Clean the .faa file to remove '*' characters
            cleaned_file_path = clean_faa_file(file_path)
            try:
                output_file_name = file_name.replace(recipient_suffix, '_recipient_prediction.txt')
                output_file_path = os.path.join(folder_path, output_file_name)
                command = f'blastp -query {cleaned_file_path} -out {output_file_path} -num_threads 32 -db resfinderfg2.0 -outfmt "10 qseqid sseqid qcovs pident length" -max_target_seqs 1'
                subprocess.run(command, shell=True, check=True)

                # Convert only the _prediction.txt files
                txt_output_file = output_file_path
                if txt_output_file.endswith('_prediction.txt') and os.path.exists(txt_output_file):
                    convert_txt_to_csv(txt_output_file, output_file_path + ".csv")

            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing recipient file {cleaned_file_path}: {e}")


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
    main_folder_path = "BLAST"

    for folder_name in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path)

if __name__ == "__main__":
    main()



# Read all CSV files into a single DataFrame
folder_path = 'BLAST'

dfs = []
# Walk through the directory tree to find all CSV files
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith('_prediction.txt.csv'):
            file_path = os.path.join(root, file)

            # Check if the file is empty before reading
            if os.path.getsize(file_path) > 0:
                df = pd.read_csv(file_path, header=None)
                dfs.append(df)
            else:
                print(f"Skipping empty file: {file_path}")



# Concatenate all DataFrames into a single DataFrame
merged_df = pd.concat(dfs, ignore_index=True)

new_header=['query','target','coverage%','identity%','align_length']
merged_df.columns = new_header

output_file_path = os.path.join(folder_path, 'summary_BLAST.csv')
merged_df.to_csv(output_file_path, index=False)

print(f"Summary CSV file saved to {output_file_path}")



#Create gene list
source_dir = 'BLAST'

gene_list = []

for root, dirs, files in os.walk(source_dir):
    for file in files:
        if file.endswith('_detected_HGTs.txt'):
            filepath = os.path.join(root, file)
            logging.info(f"Processing file: {filepath}")
            genes_in_file = set()
            try:
                with open(filepath, 'r') as txtfile:
                    reader = csv.DictReader(txtfile, delimiter='\t')
                    for row in reader:
                        genes_in_file.add(row['Gene_1'])
                        genes_in_file.add(row['Gene_2'])
            except Exception as e:
                logging.error(f"Error reading file '{filepath}': {e}")

            gene_list.extend([gene for gene in genes_in_file if gene not in gene_list])

logging.info(f"Gene list: {gene_list}")

output_file_path = os.path.join('BLAST', 'gene_list.csv')
with open(output_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(['Gene'])
    writer.writerows([[gene] for gene in gene_list])

logging.info(f"Gene list saved to {output_file_path} successfully in the BLAST folder.")
logging.info(f"Number of unique genes: {len(gene_list)}")


# Function to determine the direction (donor/recipient) and retrieve Identity
def determine_direction_and_identity(gene_id, hgt_files_folder):
    # Walk through the directory tree
    for root, dirs, files in os.walk(hgt_files_folder):
        for filename in files:
            if filename.endswith('_detected_HGTs.txt'):  # Check for specific HGT files
                file_path = os.path.join(root, filename)
                with open(file_path, 'r') as file:
                    for line in file:
                        columns = line.strip().split('\t')
                        direction_parts = columns[-1]  # Get the last part of the line (assuming it contains the direction)
                        identity = columns[2]  # Assuming the 3rd column contains the Identity
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
        # Read and process header
        header = file.readline().strip()
        updated_header = f"{header},DorR,HGT_identity"  # Add new columns 'DorR' and 'Identity'
        updated_rows.append(updated_header)

        # Process each row in the CSV file
        for line in file:
            gene_id = line.strip().split(',')[0].strip('"')  # Assuming gene_id is in the 1st column
            direction, identity = determine_direction_and_identity(gene_id, hgt_files_folder)
            if direction and identity:
                updated_row = f"{line.strip()},{direction},{identity}"  # Add direction and identity
            else:
                updated_row = f"{line.strip()},,"  # Default values if not found
            updated_rows.append(updated_row)

    # Write updated rows to a new CSV file
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
    updated_df = pd.read_csv(updated_csv_file)  # Load the updated CSV as a DataFrame
    annotation_df = pd.read_csv(annotation_file)  # Load the annotation file as a DataFrame

    # Merge the updated CSV with the annotation file based on 'COG_ID' and 'COG'
    merged_df = pd.merge(updated_df, annotation_df, left_on='target', right_on='ID', how='left')
    merged_df = merged_df.drop(merged_df.columns[[7]], axis=1)

    # Save the merged DataFrame to a new CSV file
    merged_df.to_csv('BLAST/summary_BLAST_DorR.csv', index=False)
else:
    print("Error: The updated CSV file path is None.")



print(f'The pipeline has completed, go home and sleep!')
