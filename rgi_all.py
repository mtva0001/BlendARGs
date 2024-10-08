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
target_dir = 'RGI_allbins'

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

# Counters
total_wd_folders = 0
count_hgt_folders = 0

# Walk through the directory tree starting from the source directory
for root, dirs, files in os.walk(source_dir):
    if 'RGI' in root:
        continue
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
    cleaned_file_path = file_path + ".cleaned"
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

    # Define the file suffix
    combined_fasta_suffix = "_combined_faa.fasta"

    for file_name in os.listdir(folder_path):
        # Check for the combined fasta file
        if file_name.endswith(combined_fasta_suffix):
            file_path = os.path.join(folder_path, file_name)
            logging.info(f"Processing combined file: {file_path}")

            # Clean the .fasta file to remove '*' characters
            cleaned_file_path = clean_faa_file(file_path)
            try:
                output_file_name = file_name.replace(combined_fasta_suffix, '_prediction')
                output_file_path = os.path.join(folder_path, output_file_name)
                command = f"rgi main -i {cleaned_file_path} -o {output_file_path} -t protein -n 32 --clean --local --include_loose"
                subprocess.run(command, shell=True, check=True)

                # Convert only the _prediction.txt files
                txt_output_file = output_file_path + ".txt"
                if txt_output_file.endswith('_prediction.txt') and os.path.exists(txt_output_file):
                    convert_txt_to_csv(txt_output_file, output_file_path + ".csv")

            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing combined file {cleaned_file_path}: {e}")



def convert_txt_to_csv(txt_file, csv_file):
    """
    Converts a tab-separated .txt file to a comma-separated .csv file.
    Only processes files ending with '_prediction.txt'.
    """
    try:
        with open(txt_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter=',')
            for row in reader:
                writer.writerow(row)
        logging.info(f"Converted {txt_file} to {csv_file}")

    except Exception as e:
        logging.error(f"Error converting {txt_file} to CSV: {e}")


def main():
    main_folder_path = "RGI_allbins"

    for folder_name in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path)

if __name__ == "__main__":
    main()



# Read all CSV files into a single DataFrame
folder_path = 'RGI_allbins'

dfs = []
# Walk through the directory tree to find all CSV files
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith('_prediction.csv'):
            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path)
            dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
merged_df = pd.concat(dfs, ignore_index=True)

output_file_path = os.path.join(folder_path, 'summary_RGI_allbins.csv')
merged_df.to_csv(output_file_path, index=False)

print(f"Summary CSV file saved to {output_file_path}")


