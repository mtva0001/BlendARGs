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

        # Process each *.fna file
        for filename in os.listdir(summary_folder_path):
            if filename.endswith(".fna"):
                print(f"Processing file: {filename}")
                file_path = os.path.join(summary_folder_path, filename)
                source_type = "plasmid" if "plasmid.fna" in filename else "virus"
                
                with open(file_path, 'r') as fna_file:
                    gene, sequence = None, ""
                    for line in fna_file:
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
    summary_path = os.path.join(base_path, "genomad_nucl_summary.csv")
    df.to_csv(summary_path, index=False)
    print(f"Data saved to {summary_path}")
