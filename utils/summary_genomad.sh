#This script outputs three summary files for plasmids, viruses and their protein sequences identified by geNomad.

import os
import pandas as pd

# Define the path to geNomad folder
base_path = "result_BlendARGs/VirusIdentification/geNomad"

# List to store the combined data
data = []
excluded_genes = set()  # Track excluded provirus genes

# Loop through each sample folder in the base directory
for sample_folder in os.listdir(base_path):
    sample_folder_path = os.path.join(base_path, sample_folder)
    
    # Check if it's a directory
    if os.path.isdir(sample_folder_path):
        print(f"Processing sample folder: {sample_folder}")
        
        # Extract sample name
        sample_name = sample_folder

        # Find the *_summary subfolder
        summary_folder = next((f for f in os.listdir(sample_folder_path) if f.endswith("_summary")), None)
        if not summary_folder:
            print(f"No summary folder found in {sample_folder}. Skipping...")
            continue
        
        summary_folder_path = os.path.join(sample_folder_path, summary_folder)

        # Process each *_proteins.faa file
        for filename in os.listdir(summary_folder_path):
            if filename.endswith("_proteins.faa"):
                print(f"Processing file: {filename}")
                file_path = os.path.join(summary_folder_path, filename)
                source_type = "plasmid" if "_plasmid_proteins.faa" in filename else "virus"
                
                with open(file_path, 'r') as faa_file:
                    gene, sequence = None, ""
                    for line in faa_file:
                        line = line.strip()
                        if line.startswith(">"):  # Header line
                            if gene:  # Save previous record if exists
                                if "provirus" in gene:  # Track excluded provirus genes
                                    excluded_genes.add(gene)
                                else:
                                    data.append([gene, sequence.rstrip('*'), source_type, sample_name, contig_name])
                            
                            # Process the new header
                            header_parts = line[1:].split(" ", 1)
                            gene = header_parts[0]
                            contig_name = "_".join(header_parts[0].split("_")[:2])
                            sequence = ""
                        else:
                            sequence += line
                    
                    # Append last entry if exists and not provirus
                    if gene:
                        if "provirus_" in gene:
                            excluded_genes.add(gene)
                        else:
                            data.append([gene, sequence.rstrip('*'), source_type, sample_name, contig_name])

# Check if data was collected
if not data:
    print("No data collected. Please check the file structure and content.")
else:
    # Create the DataFrame from collected data
    columns = ["gene", "sequence", "Source", "Sample", "Contig"]
    df = pd.DataFrame(data, columns=columns)
    
    # Save as genomad_summary.csv
    summary_path = os.path.join(base_path, "genomad_summary.csv")
    df.to_csv(summary_path, index=False)
    print(f"Data saved to {summary_path}")

# Initialize empty dataframes for virus and plasmid summaries
summary_virus = pd.DataFrame()
summary_plasmid = pd.DataFrame()

# Process virus and plasmid .tsv files separately to create initial summary files
for sample_folder in os.listdir(base_path):
    sample_folder_path = os.path.join(base_path, sample_folder)
    
    if os.path.isdir(sample_folder_path):
        print(f"Looking for gene files in {sample_folder}...")

        # Find the *_summary subfolder
        summary_folder = next((f for f in os.listdir(sample_folder_path) if f.endswith("_summary")), None)
        if not summary_folder:
            print(f"No summary folder found in {sample_folder}. Skipping...")
            continue

        summary_folder_path = os.path.join(sample_folder_path, summary_folder)

        # Process each *_plasmid_genes.tsv and *_virus_genes.tsv file
        for gene_file_type in ["plasmid", "virus"]:
            gene_file = next((f for f in os.listdir(summary_folder_path) if f.endswith(f"_{gene_file_type}_genes.tsv")), None)
            if gene_file:
                gene_file_path = os.path.join(summary_folder_path, gene_file)
                print(f"Reading {gene_file} from {sample_folder}...")
                
                try:
                    # Read the .tsv file and reset index to keep 'gene' as a column
                    gene_df = pd.read_csv(gene_file_path, sep='\t').reset_index()
                    gene_df.columns = gene_df.columns.str.strip()

                    # Exclude rows where 'gene' is in excluded_genes
                    gene_df = gene_df[~gene_df['gene'].isin(excluded_genes)]
                except Exception as e:
                    print(f"Error reading {gene_file_path}: {e}")
                    continue
                
                # Append gene data to either summary_virus or summary_plasmid DataFrame
                if gene_file_type == "virus":
                    summary_virus = pd.concat([summary_virus, gene_df], ignore_index=True)
                else:
                    summary_plasmid = pd.concat([summary_plasmid, gene_df], ignore_index=True)

# Save initial virus and plasmid summaries
summary_virus_path = os.path.join(base_path, "summary_virus.csv")
summary_plasmid_path = os.path.join(base_path, "summary_plasmid.csv")
summary_virus.to_csv(summary_virus_path, index=False)
summary_plasmid.to_csv(summary_plasmid_path, index=False)
print(f"Virus data saved to {summary_virus_path}")
print(f"Plasmid data saved to {summary_plasmid_path}")
