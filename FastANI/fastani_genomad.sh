#This script takes the nucleotide sequences of each plasmid and virus contig and run FastANI across them within each sampling location.
import os
import csv
import sys
import subprocess
from pathlib import Path

# Define the path to the fastANI executable
FASTANI = "/cephyr/NOBACKUP/groups/jbp/matev/OceanARM/FastANI-1.34/build/fastANI"

# Input CSV file with sequences and metadata
input_file = "result_BlendARGs/VirusIdentification/geNomad/genomad_nucl_summary.csv"

# Output directory for fastANI results
csv.field_size_limit(sys.maxsize)
output_dir = "fastani_results"
os.makedirs(output_dir, exist_ok=True)


def prepare_fastani_input(source_group, sample):
    """
    Prepare input FASTA files and genome list for fastANI analysis.

    Args:
    - source_group (str): The source group ('virus' or 'plasmid').
    - sample (str): The sample name.

    Returns:
    - temp_dir (Path): Path to the temporary directory for this source group and sample.
    """
    temp_dir = Path(output_dir) / f"{sample}_{source_group}_temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing Source Group: {source_group}, Sample: {sample}")

    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Source'] == source_group and row['Sample'] == sample:
                fasta_path = temp_dir / f"{row['gene']}.fasta"
                with open(fasta_path, 'w') as fasta_file:
                    fasta_file.write(f">{row['gene']}\n")
                    fasta_file.write(f"{row['sequence'].strip()}\n")

    genome_list_path = temp_dir / "genome_list.txt"
    with open(genome_list_path, 'w') as genome_list:
        for fasta_file in temp_dir.glob("*.fasta"):
            genome_list.write(str(fasta_file.resolve()) + "\n")

    print(f"Genome list created: {genome_list_path}")
    return temp_dir


def run_fastani(temp_dir, output_file):
    """
    Run fastANI analysis.

    Args:
    - temp_dir (Path): Path to the temporary directory containing the genome list.
    - output_file (Path): Path to the output file for fastANI results.
    """
    genome_list = temp_dir / "genome_list.txt"

    if genome_list.exists() and genome_list.stat().st_size > 0:
        print(f"Running fastANI for genome list: {genome_list}")
        subprocess.run(
            [FASTANI, "--rl", str(genome_list), "--ql", str(genome_list), "-t", "64", "-o", str(output_file)],
            check=True
        )
    else:
        print(f"No sequences found for this group/sample. Skipping: {genome_list}")


def process_fastani_output(output_file, sample, source, combined_output_file):
    """
    Process fastANI output and append results to a combined output file.

    Args:
    - output_file (Path): Path to the fastANI result file.
    - sample (str): Sample name.
    - source (str): Source group.
    - combined_output_file (Path): Path to the combined output file.
    """
    if output_file.exists():
        with open(output_file, 'r') as infile, open(combined_output_file, 'a') as outfile:
            for line in infile:
                columns = line.strip().split('\t')
                if len(columns) >= 3 and columns[2].replace('.', '', 1).isdigit():
                    outfile.write(f"{columns[0]},{columns[1]},{columns[2]},{sample},{source}\n")


def main():
    # Define sources and initialize combined output files
    sources = ["virus", "plasmid"]
    combined_output_files = {
        "virus": Path(output_dir) / "fastani_combined_virus.csv",
        "plasmid": Path(output_dir) / "fastani_combined_plasmid.csv",
    }

    for source in sources:
        with open(combined_output_files[source], 'w') as f:
            f.write("bin_q,bin_t,ANI,Sample,Source\n")

    # Extract unique sample names from the input file
    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        samples = sorted({row['Sample'] for row in reader})

    # Process each combination of source and sample
    for sample in samples:
        for source in sources:
            print(f"Processing Sample: {sample}, Source: {source}")

            # Prepare input files for fastANI
            temp_dir = prepare_fastani_input(source, sample)

            # Define fastANI output file
            output_file = Path(output_dir) / f"{sample}_{source}_fastani_output.csv"

            # Run fastANI
            run_fastani(temp_dir, output_file)

            # Process fastANI output
            process_fastani_output(output_file, sample, source, combined_output_files[source])

            # Clean up temporary directory
            for file in temp_dir.iterdir():
                file.unlink()
            temp_dir.rmdir()

    print(f"FastANI processing complete. Combined outputs saved in {output_dir}.")


if __name__ == "__main__":
    main()
