#!/bin/bash
# Function to process each subfolder
output_dir="./MetaCHIP_results"
mkdir -p "$output_dir"

process_subfolder() {
    local subfolder="$1"
    local grouping_file="$2"
    local output_prefix="$3"
    local log_file="$output_dir/$(basename "$subfolder")_log.txt"  # Define path to log file for this subfolder

    # Check existing MetaCHIP output(s)
    if [ -f "$log_file" ]; then
        echo "Skipping $subfolder because $log_file already exists."
        return 0
    fi

    # Run the MetaCHIP PI command with the files from the subfolder
    echo "Running MetaCHIP PI for $subfolder..."
    MetaCHIP PI -i "$subfolder" -taxon gtdbtk_summary_metachip.tsv -g "$grouping_file" -x fa -t 32 -p "$output_prefix" > "$log_file" 2>&1

    # Check exit status of MetaCHIP PI command
    if [ $? -ne 0 ]; then
        echo "Error occurred during MetaCHIP PI execution. Check log file: $log_file"
        return 1
    fi

    # Run the MetaCHIP BP command for the same subfolder after MetaCHIP PI finishes
    MetaCHIP BP -g "$grouping_file" -t 32 -p "$output_prefix" -pfr > "$log_file.bp" 2>&1

    # Check exit status of MetaCHIP BP command
    if [ $? -ne 0 ]; then
        echo "Error occurred during MetaCHIP BP execution. Check log file: $log_file.bp"
        return 1
    fi
}

# Path to the directory containing the grouping files
grouping_dir="."

# Path to the directory containing the .fa files and subfolders
genome_folder="./bins"

# Initialize counter for processed subfolders
processed_subfolders=0

# Count total number of subfolders
total_subfolders=$(find "$genome_folder" -mindepth 1 -maxdepth 1 -type d | wc -l)

# Iterate over the grouping files
for grouping_file in "$grouping_dir"/*_info.txt; do
    if [ -f "$grouping_file" ]; then
        # Extract the base filename (without extension and _info suffix)
        base_filename=$(basename "${grouping_file%.*}" | sed 's/_info$//')

        # Modify the -p parameter with the corresponding name
        output_prefix="$base_filename"

        # Find subfolders within the genome folder that match the base filename
        matching_subfolders=()
        while IFS= read -r -d '' subfolder; do
            matching_subfolders+=("$subfolder")
        done < <(find "$genome_folder" -mindepth 1 -maxdepth 1 -type d -name "$base_filename" -print0)

        # Check if matching subfolders were found
        if [ ${#matching_subfolders[@]} -eq 0 ]; then
            echo "No matching subfolders found for grouping file: $grouping_file"
            continue  # Move to the next grouping file
        fi

        # Process each subfolder in parallel
        for subfolder in "${matching_subfolders[@]}"; do
            process_subfolder "$subfolder" "$grouping_file" "$output_prefix" &
            ((processed_subfolders++))
            echo "Processing $processed_subfolders out of $total_subfolders subfolders."
        done
        # Wait for all background processes to finish before moving to the next grouping file
        wait
    fi
done


