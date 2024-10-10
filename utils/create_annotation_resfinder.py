#This script creates an annotation csv file that can be used to link the ARG IDs (reported in the blastp output files) to ARG classes, sources and antibiotics.

import os
import pandas as pd

# Path to the ResFinderFG2.0 AA file
faa_file_path = 'ResFinder_FG_AA.faa'
# Path to the output CSV file
output_csv_path = 'annotation_resfinder.csv'

headers_data = []

with open(faa_file_path, 'r') as file:
    for line in file:
        if line.startswith('>'):
            header = line[1:].strip().split('|')
            headers_data.append(header)


headers_df = pd.DataFrame(headers_data, columns=['ARG', 'ID', 'Source', 'AB'])

headers_df.to_csv(output_csv_path, index=False)

print(f"Annotation CSV file saved to {output_csv_path}")
