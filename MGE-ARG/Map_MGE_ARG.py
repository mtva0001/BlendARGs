# This script uses MGEs identified by geNomad and maps to the newly compailed db of detected ARGs.

import pandas as pd
import subprocess
import os

# Paths to the input files
query_file = "../BlendARGs/result_BlendARGs/VirusIdentification/geNomad/genomad_summary.csv"
target_file = "./ARG_AA_summary.csv"
final_output_file = "blastp_MGE_ARG_results.csv"


# Read in query and target data
query_df = pd.read_csv(query_file)
target_df = pd.read_csv(target_file)

# Function to extract sample name from the 'Query' column
def get_sample_from_query(query):
    if isinstance(query, str) and query.strip():
        return query.split('_')[0]
    else:
        return "Unknown"

# Apply the sample extraction to the 'Sample' column in both dataframes
query_df['Sample'] = query_df['Sample']  # Sample column already present in query_df
target_df['Sample'] = target_df['Query'].apply(get_sample_from_query)

# Store results from each sample for final combination
all_results = []

# Process each sample group individually
for sample in query_df['Sample'].unique():
    print(f"Processing sample: {sample}")
    
    # Filter sequences for the current sample
    sample_query_df = query_df[query_df['Sample'] == sample]
    sample_target_df = target_df[target_df['Sample'] == sample]
    
    # Skip if there are no corresponding target sequences for the sample
    if sample_target_df.empty:
        print(f"No target sequences found for sample: {sample}")
        continue

    # Save target sequences to a FASTA file for makeblastdb
    target_fasta = f"{sample}_target.fasta"
    with open(target_fasta, "w") as f:
        for _, row in sample_target_df.iterrows():
            f.write(f">{row['Query']}\n{row['Protein_Sequence']}\n")

    # Run makeblastdb to create a database from target sequences
    target_db = f"{sample}_target_db"
    subprocess.run(["makeblastdb", "-in", target_fasta, "-dbtype", "prot", "-out", target_db])

    # Save query sequences to a FASTA file for blastp
    query_fasta = f"{sample}_query.fasta"
    with open(query_fasta, "w") as f:
        for _, row in sample_query_df.iterrows():
            f.write(f">{row['gene']}\n{row['sequence']}\n")

    # Run blastp to get top 10 hits per query
    blast_output = f"{sample}_blast_output.tsv"
    subprocess.run([
        "blastp", "-query", query_fasta, "-db", target_db, "-evalue", "1e-10",
        "-out", blast_output, "-outfmt", "6 qseqid sseqid qcovs pident length bitscore", "-max_target_seqs", "10"
    ])

    # Check if the BLAST output file is empty or nonexistent
    if not os.path.exists(blast_output) or os.stat(blast_output).st_size == 0:
        print(f"No BLAST hits for sample {sample}. Skipping.")
        os.remove(target_fasta)
        os.remove(query_fasta)
        # Optionally, remove the BLAST database files
        for ext in [".phr", ".pin", ".psq"]:
            db_file = f"{target_db}{ext}"
            if os.path.exists(db_file):
                os.remove(db_file)
        continue

    # Load the BLAST output
    blast_df = pd.read_csv(blast_output, sep="\t", header=None)
    blast_df.columns = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']

    # Sort by query and bitscore (descending) to keep only the highest bitscore hit per query
    blast_df = blast_df.sort_values(by=['query', 'bitscore'], ascending=[True, False])
    best_hits_df = blast_df.drop_duplicates(subset='query', keep='first')

    # Apply the filter for identity >= 95% and coverage >= 70%
    filtered_blast_df = best_hits_df[(best_hits_df["identity%"] >= 95) & (best_hits_df["coverage%"] >= 70)]
    
    # Add sample information and store in all_results
    filtered_blast_df["Sample"] = sample
    all_results.append(filtered_blast_df)
    
    # Clean up temporary files
    os.remove(target_fasta)
    os.remove(query_fasta)
    os.remove(blast_output)
    # Optionally, remove the BLAST database files
    for ext in [".phr", ".pin", ".psq"]:
        db_file = f"{target_db}{ext}"
        if os.path.exists(db_file):
            os.remove(db_file)

# Combine all results into a single DataFrame and save to CSV
final_df = pd.concat(all_results, ignore_index=True)
final_df.to_csv(final_output_file, index=False)
print(f"Final combined BLAST results saved to {final_output_file}.")
