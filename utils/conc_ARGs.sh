# This bash script finds the ARG-related output files and detected protein sequences and compiles them together for the HMMER analysis.

#!/bin/bash

# Define file paths
RF_CSV="./BLAST_allbins/summary_BLAST_allbins_RF_ARG.csv"
Latent_CSV="./BLAST_allbins_latent/summary_BLAST_allbins_Latent_ARG.csv"
RGI_CSV="./RGI_allbins/summary_RGI_allbins_ARG.csv"
RF_FAA="ResFinder_FG_AA_fixed.faa"
Latent_FAA="Latent_ARGs_AA.fasta"
OUTPUT="ARG_AA_summary.csv"

# Initialize output file with a header (modify if needed)
echo "Query,Protein_Sequence,Source" > "$OUTPUT"

# Function to extract protein sequence from .faa file
extract_protein_sequence() {
    local faa_file=$1
    local header=$2
    awk -v header="$header" -v seq="" 'BEGIN {found=0}
        # Locate the matching header
        $1 == (">"header) { found=1; next }
        # Collect sequence if the header was found
        found == 1 && $1 ~ /^>/ { exit }
        found == 1 { seq = seq $0 }
        END { print seq }' "$faa_file"
}

# Process RF CSV file
while IFS=',' read -r query target rest; do
    protein_sequence=$(extract_protein_sequence "$RF_FAA" "$target")
    if [[ -n "$protein_sequence" ]]; then
        echo "$query,$protein_sequence,ResFinder" >> "$OUTPUT"
    fi
done < <(tail -n +2 "$RF_CSV")  # Skip header row

# Process Latent CSV file
while IFS=',' read -r query target rest; do
    protein_sequence=$(extract_protein_sequence "$Latent_FAA" "$target")
    if [[ -n "$protein_sequence" ]]; then
        echo "$query,$protein_sequence,Latent" >> "$OUTPUT"
    fi
done < <(tail -n +2 "$Latent_CSV")  # Skip header row

# Process RGI CSV file (using ORF_ID as the query)
while IFS=',' read -r orf_id predicted_protein; do
    echo "$orf_id,$predicted_protein,RGI" >> "$OUTPUT"
done < <(awk -F',' 'NR==1 {next} {print $1","$19}' "$RGI_CSV")
