#!/bin/bash

# Define the path to the fastANI executable
FASTANI="/cephyr/NOBACKUP/groups/jbp/matev/OceanARM/FastANI-1.34/build/fastANI"

# Input CSV file with sequences and metadata
input_file="result_LakesARGs/VirusIdentification/geNomad/genomad_nucl_summary.csv"

# Output directory for fastANI results
output_dir="fastani_results"
mkdir -p "$output_dir"

# Output files for combined results
virus_combined="$output_dir/virus_fastani_combined.csv"
plasmid_combined="$output_dir/plasmid_fastani_combined.csv"

# Initialize the combined files with headers
echo "bin_q,bin_t,ANI,Sample,Source" > "$virus_combined"
echo "bin_q,bin_t,ANI,Sample,Source" > "$plasmid_combined"

# Function to extract and prepare input for fastANI analysis
prepare_fastani_input() {
    local source_group="$1"
    local sample="$2"
    local temp_dir="$output_dir/${sample}_${source_group}_temp"
    mkdir -p "$temp_dir"

    # Extract sequences for the specific source group and sample
    awk -F"," -v source="$source_group" -v sample="$sample" \
        'NR > 1 && $5 == source && $4 == sample { \
            print ">" $1 "\n" $2 > "'"$temp_dir/"$1".fasta"'"}' "$input_file"

    # Create a genome list for fastANI
    find "$temp_dir" -name "*.fasta" > "$temp_dir/genome_list.txt"
    echo "$temp_dir"
}

# Function to run fastANI analysis
run_fastani() {
    local genome_list="$1/genome_list.txt"
    local output_file="$2"

    if [ -s "$genome_list" ]; then
        "$FASTANI" --rl "$genome_list" --ql "$genome_list" -t 32 -o "$output_file"
    else
        echo "No sequences found for this group/sample. Skipping."
    fi
}

# Process each combination of Source and Sample
sources=("virus" "plasmid")

# Extract unique Sample names
mapfile -t samples < <(awk -F"," 'NR > 1 {print $4}' "$input_file" | sort | uniq)

for sample in "${samples[@]}"; do
    for source in "${sources[@]}"; do
        echo "Processing Sample: $sample, Source: $source"

        # Prepare input files for fastANI
        temp_dir=$(prepare_fastani_input "$source" "$sample")

        # Define temporary fastANI output file
        output_file="$output_dir/${sample}_${source}_fastani_output.txt"

        # Run fastANI
        run_fastani "$temp_dir" "$output_file"

        # Append fastANI output to the combined file
        if [ -f "$output_file" ]; then
            combined_file="$output_dir/${source}_fastani_combined.csv"
            awk -F"\t" -v sample="$sample" -v source="$source" \
                'NF >= 3 && $3 ~ /^[0-9.]+$/ {
                    print $1 "," $2 "," $3 "," sample "," source
                 }' "$output_file" >> "$combined_file"
        fi

        # Clean up temporary directory
        rm -rf "$temp_dir"
    done
done

echo "FastANI processing complete. Combined results saved in $output_dir."
