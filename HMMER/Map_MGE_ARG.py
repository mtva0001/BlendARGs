# This script uses MGEs identified by geNomad and maps to the detected ARGs by the three applied databases.

import os
import subprocess
import pandas as pd

# Define paths
query_file = "../BlendARGs/result_BlendARGs/VirusIdentification/geNomad/genomad_summary.csv"
target_file = "./ARG_AA_summary.csv"
hmmer_dir = "temp_hmmer"  # Directory to store HMMER results
alignment_tool = "mafft"  # Use MAFFT for alignment
final_output_file = "compiled_hmmer_results.csv"  # Final compiled output file

# Create directories if they don't exist
os.makedirs(hmmer_dir, exist_ok=True)

# Read the query and target dataframes
query_df = pd.read_csv(query_file)
target_df = pd.read_csv(target_file)

# Prepare the target dataframe by extracting the sample group from the Query column
def get_sample_from_query(query):
    if isinstance(query, str) and query.strip():
        return query.split('_')[0]
    else:
        return "Unknown"

target_df['Sample'] = target_df['Query'].apply(get_sample_from_query)

# Initialize a list to collect all HMMER results
compiled_results = []

# Iterate over samples in the query and target data
for sample in target_df['Sample'].unique():
    # Filter the query sequences based on the sample
    sample_query_df = query_df[query_df['Sample'] == sample]
    
    # Prepare FASTA file for this sample
    temp_fasta = os.path.join(hmmer_dir, f"{sample}_target.fasta")
    with open(temp_fasta, 'w') as fasta_file:
        for index, row in sample_query_df.iterrows():
            fasta_file.write(f">{row['gene']}\n{row['sequence']}\n")

    # Align the sequences using MAFFT
    aligned_fasta = os.path.join(hmmer_dir, f"{sample}_aligned.fasta")
    subprocess.run([alignment_tool, "--auto", temp_fasta], stdout=open(aligned_fasta, 'w'))
    
    # Build the HMM profile from the aligned sequences
    profile_file = os.path.join(hmmer_dir, f"{sample}_profile.hmm")
    subprocess.run(["hmmbuild", profile_file, aligned_fasta])

    # Now perform hmmsearch with the constructed HMM profile
    hmmer_output = os.path.join(hmmer_dir, f"{sample}_output.tbl")
    subprocess.run(["hmmsearch", "--domtblout", hmmer_output, profile_file, temp_fasta, "--cpu 32"])

    # Process the HMMER output (e.g., filtering by E-value)
    hmmer_df = pd.read_csv(hmmer_output, sep="\s+", comment="#", header=None)
    
    # Read HMMER output dynamically based on number of columns
    hmmer_df = pd.read_csv(hmmer_output, sep="\s+", comment="#", header=None)

    # Define column names for the domtblout format
    hmmer_df.columns = [
        "target_name", "t_accession", "tlen", "query_name", "q_accession", "qlen",
        "E_value", "score_sequence", "bias_sequence", "#",
        "of", "c-Evalue", "i-Evalue", "score", "bias",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc", "description"
    ]

    # Filter by E-value < 1e-3
    hmmer_df_filtered = hmmer_df[(hmmer_df["E_value_domain"] < 1e-3) & (hmmer_df["acc"] >= 0.95)]

    # Add the sample column for identification
    hmmer_df_filtered['Sample'] = sample

    # Append the filtered results to the compiled list
    compiled_results.append(hmmer_df_filtered)

    print(f"Processed sample {sample} - HMMER results saved.")

# Concatenate all results into a single dataframe
compiled_df = pd.concat(compiled_results, ignore_index=True)

# Save the final compiled dataframe to a CSV file
compiled_df.to_csv(final_output_file, index=False)

print(f"All HMMER results have been compiled and saved to {final_output_file}")
