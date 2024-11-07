#This script create two summary files for plasmids and viruses identified by geNomad
#!/bin/bash

# Set paths for your geNomad directory and output summary files
genomad_dir="/VirusIdentification/geNomad"
virus_summary="virus_summary.tsv"
plasmid_summary="plasmid_summary.tsv"

# Clear or create the summary files with headers
echo -e "Sample\t$(head -n 1 $(find $genomad_dir -type f -name "*_virus_summary.tsv" | head -n 1))" > $virus_summary
echo -e "Sample\t$(head -n 1 $(find $genomad_dir -type f -name "*_plasmid_summary.tsv" | head -n 1))" > $plasmid_summary

# Loop through each sample folder and process the .tsv files
for sample_folder in "$genomad_dir"/*; do
    sample_name=$(basename "$sample_folder")

    # Paths to the virus and plasmid summary files
    virus_file="$sample_folder/MEGAHIT-${sample_name}.contigs_summary/MEGAHIT-${sample_name}.contigs_virus_summary.tsv"
    plasmid_file="$sample_folder/MEGAHIT-${sample_name}.contigs_summary/MEGAHIT-${sample_name}.contigs_plasmid_summary.tsv"

    # Append data to virus summary file if it exists
    if [ -f "$virus_file" ]; then
        tail -n +2 "$virus_file" | awk -v sample="$sample_name" '{print sample"\t"$0}' >> $virus_summary
    fi

    # Append data to plasmid summary file if it exists
    if [ -f "$plasmid_file" ]; then
        tail -n +2 "$plasmid_file" | awk -v sample="$sample_name" '{print sample"\t"$0}' >> $plasmid_summary
    fi
done

echo "Summaries compiled into $virus_summary and $plasmid_summary."
