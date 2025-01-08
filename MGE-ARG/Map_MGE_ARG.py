# This script uses MGEs identified by geNomad and maps to the newly compiled db of detected ARGs.

import pandas as pd
import subprocess
import os

# Paths to the input files
query_file = "../BlendARGs/result_BlendARGs/VirusIdentification/geNomad/genomad_summary.csv"
target_file = "./ARG_AA_summary.csv"
group_file = "sample_groups.txt"  # File specifying group-sample relationships (e.g., grouping vertical samples within sampling location)
final_output_file = "blastp_MGE_ARG_results.csv"

# Read in query, target, and group data
query_df = pd.read_csv(query_file)
target_df = pd.read_csv(target_file)
group_df = pd.read_csv(group_file)  # Assuming it has headers: Group,Sample

# Function to extract sample name from the 'Query' column
def get_sample_from_query(query):
    if isinstance(query, str) and query.strip():
        return query.split('_')[0]
    else:
        return "Unknown"

# Apply the sample extraction to the 'Sample' column in both dataframes
query_df['Sample'] = query_df['Sample']  # Sample column already present in query_df
target_df['Sample'] = target_df['Query'].apply(get_sample_from_query)

# Add Group information to query and target DataFrames
query_df = query_df.merge(group_df, on='Sample', how='left')
target_df = target_df.merge(group_df, on='Sample', how='left')

# Store results from each sample for final combination
all_results = []

# Process each sample group individually
for group in query_df['Group'].unique():
    print(f"Processing group: {group}")
    
    # Get samples in the current group
    group_samples = group_df[group_df['Group'] == group]['Sample'].tolist()
    
    # Filter sequences for the current group
    group_query_df = query_df[query_df['Sample'].isin(group_samples)]
    group_target_df = target_df[target_df['Sample'].isin(group_samples)]
    
    # Skip if there are no corresponding target sequences for the group
    if group_target_df.empty:
        print(f"No target sequences found for group: {group}")
        continue

    # Save target sequences to a FASTA file for makeblastdb
    target_fasta = f"{group}_target.fasta"
    with open(target_fasta, "w") as f:
        for _, row in group_target_df.iterrows():
            f.write(f">{row['Query']}\n{row['Protein_Sequence']}\n")

    # Run makeblastdb to create a database from target sequences
    target_db = f"{group}_target_db"
    subprocess.run(["makeblastdb", "-in", target_fasta, "-dbtype", "prot", "-out", target_db])

    # Process each sample in the group as a separate query
    for sample in group_samples:
        print(f"Processing sample: {sample} in group: {group}")
        
        # Filter query sequences for the current sample
        sample_query_df = group_query_df[group_query_df['Sample'] == sample]
        
        # Skip if no queries exist for the current sample
        if sample_query_df.empty:
            print(f"No query sequences found for sample: {sample}")
            continue

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

        # Load the BLAST output
        if os.path.exists(blast_output) and os.path.getsize(blast_output) > 0:
            blast_df = pd.read_csv(blast_output, sep="\t", header=None)
            blast_df.columns = ['query', 'target', 'coverage%', 'identity%', 'align_length', 'bitscore']

            # Sort by query and bitscore (descending) to keep only the highest bitscore hit per query
            blast_df = blast_df.sort_values(by=['query', 'bitscore'], ascending=[True, False])
            best_hits_df = blast_df.drop_duplicates(subset='query', keep='first')

            # Apply the filter for identity >= 95% and coverage >= 70%
            filtered_blast_df = best_hits_df[(best_hits_df["identity%"] >= 95) & (best_hits_df["coverage%"] >= 70)]
            
            # Add sample and group information and store in all_results
            filtered_blast_df["Sample"] = sample
            filtered_blast_df["Group"] = group
            all_results.append(filtered_blast_df)
        else:
            print(f"No BLAST results for sample: {sample}")

        # Clean up temporary files
        os.remove(query_fasta)
        os.remove(blast_output)

    # Optionally, remove the BLAST database files
    for ext in [".phr", ".pin", ".psq"]:
        db_file = f"{target_db}{ext}"
        if os.path.exists(db_file):
            os.remove(db_file)
    os.remove(target_fasta)

# Combine all results into a single DataFrame and save to CSV
if all_results:
    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(final_output_file, index=False)
    print(f"Final combined BLAST results saved to {final_output_file}.")
else:
    print("No results to save.")
