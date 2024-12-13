#!/bin/bash

# Define the path to the fastANI executable
FASTANI="/cephyr/NOBACKUP/groups/jbp/matev/OceanARM/FastANI-1.34/build/fastANI"

# Input GTDB-Tk summary file
gtdbtk_file="gtdbtk_summary_metachip.tsv"

# Output directory for fastANI results
output_dir="fastani_results"
mkdir -p "$output_dir"

# Output file for concatenated results
summary_file="fastani_summary.csv"
echo "bin_q,bin_t,ANI,Sample,source,target" > "$summary_file"

# Read user_genome values from the GTDB-Tk summary file
mapfile -t user_genomes < <(awk -F"\t" 'NR>1 {print $1}' "$gtdbtk_file")

# Function to check if a file matches the user_genomes list
file_matches_user_genome() {
    local file_basename=$(basename "$1" .fa)
    for genome in "${user_genomes[@]}"; do
        if [[ "$file_basename" == "$genome" ]]; then
            return 0  # Match found
        fi
    done
    return 1  # No match
}

# Loop through each sample-specific subdirectory in the bins directory
for sample_dir in bins/*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")  # Extract sample-specific string
        genome_list="${sample_name}_genome_list.txt"  # Temporary genome list file

        # Create genome list by matching files in the sample directory with user_genomes
        > "$genome_list"
        for file in "$sample_dir"/*.fa; do
            if file_matches_user_genome "$file"; then
                echo "$file" >> "$genome_list"
            fi
        done

        # Only run fastANI if the genome list is not empty
        if [ -s "$genome_list" ]; then
            "$FASTANI" --rl "$genome_list" --ql "$genome_list" -t 32 -o "${output_dir}/${sample_name}_fastani_output.txt"
        fi

        # Remove the temporary genome list file
        rm "$genome_list"

        # If fastANI output exists, process it and validate lines before appending
        output_file="${output_dir}/${sample_name}_fastani_output.txt"
        if [ -f "$output_file" ]; then
            awk -v sample="$sample_name" 'NF >= 3 && $3 ~ /^[0-9.]+$/ {
                split($1, q_parts, "/"); split($2, t_parts, "/");
                split(q_parts[length(q_parts)], q_id, "_"); split(t_parts[length(t_parts)], t_id, "_");
                print $1 "," $2 "," $3 "," sample "," q_id[1] "," t_id[1]
            }' "$output_file" >> "$summary_file"
        fi
    fi

done

echo "FastANI processing and concatenation complete. Results saved in $summary_file."
