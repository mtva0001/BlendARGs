import pandas as pd
import os
import gzip
import shutil
import csv

directory = "./bins"

# Reformat filenames
print("Reformatting filenames...")
for filename in os.listdir(directory):
    if filename.startswith("MEGAHIT-MetaBAT2-"):
        new_filename = filename.replace("MEGAHIT-MetaBAT2-", "")
        os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))

# Replace dots with underscores in filenames
print("Replacing dots with underscores in filenames...")
for filename in os.listdir(directory):
    if filename.endswith(".fa.gz"):
        base = filename[:-6]  # Remove the last ".fa.gz"
        ext = ".fa.gz"  # The extension to be added back
        new_base = base.replace('.', '_')
        
        new_filename = new_base + ext
        
        os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))

print(f"Filenames have been reformatted and processed in {directory}")


#Reformatting taxonomy sample names
input_file = 'gtdbtk_summary.tsv'
output_file = 'gtdbtk_summary_metachip.tsv'

# Check if the output file already exists
if os.path.exists(output_file):
    print(f"{output_file} already exists. Skipping reformatting taxonomy.")
else:
    # Open input and output files
    try:
        with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')

            # Process rows
            for row in reader:
                print("Original Row:", row)  # Debugging: print original row
                if row[1] and not row[1].startswith('d__Archaea'):
                    # Renaming the value in the first column
                    row[0] = row[0].replace("MEGAHIT-MetaBAT2-", "").replace(".fa", "").replace('.', '_')
                    # Write the modified row to the output file
                    writer.writerow(row[:2])  # Keeping only the first two columns
                    print("Modified Row:", row[:2])  # Debugging: print modified row

    except Exception as e:
        print("An error occurred:", e)


#Creating groups
samplesheet_df = pd.read_csv("Samplesheet_input.csv") #Samplesheet that was used for the nf-core pipeline

gtdbtk_df = pd.read_csv("gtdbtk_summary_metachip.tsv", sep='\t')

gtdbtk_df['sample'] = gtdbtk_df['user_genome'].str.rsplit('_', n=1).str[0]

merged_df = pd.merge(samplesheet_df[['sample', 'group']], gtdbtk_df[['sample', 'user_genome']], on='sample')
merged_df = merged_df[['group', 'user_genome', 'sample']]

merged_df.to_csv("grouping_file.txt", sep=',', index=False, header=False)

for group, group_df in merged_df.groupby('group'):
    group_df = group_df[['sample', 'user_genome']]
    group_df.columns = ['sample', 'user_genome']
    group_df.to_csv(f"{group}_info.txt", sep=',', index=False, header=False)


#Create subfolder for bins
genome_folder_path = './bins'
group_file_path = 'grouping_file.txt'

genome_id_to_subfolder = {}
try:
    with open(group_file_path, 'r') as group_file:
        for line in group_file:
            fields = line.strip().split(',')
            if len(fields) >= 2:
                genome_id = fields[1].strip()  # Assuming genome ID is in the second column
                subfolder_name = fields[0].strip()  # Assuming subfolder name is in the first column
                genome_id_to_subfolder[genome_id] = subfolder_name
            else:
                print("Error: Line does not have enough fields:", line.strip())
except FileNotFoundError:
    print("Error: Grouping file not found:", group_file_path)
    exit()

for filename in os.listdir(genome_folder_path):
    if filename.endswith('.fa.gz'):
        # Get the genome ID from the filename (without extension)
        genome_id = os.path.splitext(os.path.splitext(filename)[0])[0]
        # Check if the genome ID exists in the dictionary
        if genome_id in genome_id_to_subfolder:
            # Create subfolder based on genome ID
            subfolder_path = os.path.join(genome_folder_path, genome_id_to_subfolder[genome_id])
            os.makedirs(subfolder_path, exist_ok=True)
            # Move the file into the subfolder
            shutil.move(os.path.join(genome_folder_path, filename), os.path.join(subfolder_path, filename))
            print("Moved", filename, "to", subfolder_path)

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith(".fa.gz"):
            # Full path to the gzipped file
            gz_file_path = os.path.join(root, file)
            
            # Path for the decompressed file
            uncompressed_file_path = os.path.splitext(gz_file_path)[0]
            
            # Print status
            print(f"Uncompressing: {gz_file_path}")

            # Uncompress the file
            with gzip.open(gz_file_path, 'rb') as f_in:
                with open(uncompressed_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            os.remove(gz_file_path)

print("All .fa.gz files have been uncompressed.")
