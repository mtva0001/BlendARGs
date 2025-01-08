#Taking the genomad_nucl_summary.csv to create a fasta file for PhaBOX v2.0 typing.
import csv
import sys

csv.field_size_limit(sys.maxsize)

# Input and output file paths
input_file = "result_LakesARGs/VirusIdentification/geNomad/genomad_nucl_summary.csv"  # Replace with your actual file name
output_fasta = "virus_sequences.fasta"

# Open the input CSV file and the output FASTA file
with open(input_file, "r") as csvfile, open(output_fasta, "w") as fasta_file:
    reader = csv.DictReader(csvfile)
    
    # Iterate through rows and filter for "virus" in the "Source" column
    for row in reader:
        if row["Source"].lower() == "virus":  # Case-insensitive match for "virus"
            # Create the header with gene and Sample values
            header = f">{row['gene']}_{row['Sample']}"
            sequence = row["sequence"]
            
            # Write the header and sequence to the FASTA file
            fasta_file.write(f"{header}\n{sequence}\n")

print(f"FASTA file created: {output_fasta}")
