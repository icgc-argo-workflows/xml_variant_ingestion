#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  Ontario Institute for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Guanqiao Feng
"""

import xml.etree.ElementTree as ET
import pandas as pd
import argparse
import re
from datetime import date
from pyfaidx import Fasta

def create_vcf_header(date_str, chrs, chr_dic, input_file_name):

    headers = [
        '##fileformat=VCFv4.2',
        f'##fileDate={date_str}',
        f'##source={input_file_name}',
        '##reference=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    ]

    # Remove duplicates while preserving order
    unique_chrs = list(dict.fromkeys(chrs))

    # Add contig lines based on chrs and chr_dic
    for chr_id in unique_chrs:
        if chr_id in chr_dic:
            length = chr_dic[chr_id]
            headers.append(f'##contig=<ID={chr_id},length={length}>')

    # Add INFO headers
    headers.extend([
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END_POS,Number=1,Type=Integer,Description="End position of this structural variant">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of copy number alteration">',
        '##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number">',
        '##INFO=<ID=RATIO,Number=1,Type=Float,Description="Copy number ratio">'
    ])

    return headers

def extract_data_from_cna(cna):
    return {
        # Extract useful information from the xml file
        'copy-number': cna.get('copy-number'),
        'position': cna.get('position'),
        'ratio': cna.get('ratio'),
        'type': cna.get('type')
    }

def chr_pos_order(chr_list): # this code is specific for hg19 genome assembly it sort chr based on "chr1, chr2, chrX, chrY, chrUn_, chrM"
    chr_sort = {}
    chr_order = []
    for i in chr_list:
        match = re.search(r'chr(\d)(?=\D|$)', i)
        if match:
            digit = match.group(1)
            new_i = re.sub(r'chr(\d)(?=\D|$)', rf'chr0{digit}', i)
            chr_sort[i] = new_i
        elif re.match(r'chrX', i):
            chr_sort[i] = 'chr66X'
        elif re.match(r'chrY', i):
            chr_sort[i] = 'chr77Y'
        elif re.search(r'chrUn', i):
            chr_sort[i] = re.sub(r'chrUn', r'chr88Un', i)
        elif re.match(r'chrM', i):
            chr_sort[i] = 'chr99M'
        else:
            chr_sort[i] = i

    # Sort the dictionary based on its values and then convert the keys to a list
    chr_order = sorted(chr_sort, key=lambda k: chr_sort[k])

    return chr_order

# extract ref nucletide
def extract_nucleotide(fasta_file, pos):
    fasta = Fasta(fasta_file, as_raw=True)
    sequence = fasta[pos.split(':')[0]][int(pos.split(':')[1])-1].upper()
    fasta.close()
    return sequence

# Process dataframe by adding "ID", "QUAL", "FILTER" with default value "." and "INFO" with "DP=#;AF=#"
# where DP means depth and AF means allele frequency
def process_dataframe(df, chr_list, fasta_file):

    df['#CHROM'] = df['position'].apply(lambda x: x.split(':')[0])
    df['POS'] = df['position'].apply(lambda x: x.split(':')[1].split('-')[0]).astype(int)
    df['ID'] = '.'
    df['ref_pos'] = df['position'].apply(lambda x: x.split('-')[0])
    df['REF'] = df['ref_pos'].apply(lambda x: extract_nucleotide(fasta_file, x))
    df['ALT'] = '<CNV>'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = 'SVTYPE=CNV;END=' + df['position'].apply(lambda x: x.split('-')[1]) + ';TYPE=' + df['type'] + ';CN=' + df['copy-number'] + ';RATIO=' + df['ratio']
    extract_df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

    sorted_df = chr_pos_sort(extract_df.copy(), chr_list)

    return sorted_df,sorted_df['#CHROM']

def chr_pos_sort(df, chr_ordered): # sort dataframe based on chr and pos

    cat_type = pd.CategoricalDtype(categories=chr_ordered, ordered=True)

    # Convert the #CHROM column to the categorical data type
    df['#CHROM'] = df['#CHROM'].astype(cat_type)

    # Sort the DataFrame first by '#CHROM' and then by 'POS'
    sorted_df = df.sort_values(by=['#CHROM', 'POS'])

    # Convert the chromosome numbers back to their original format
    sorted_df['#CHROM'] = sorted_df['#CHROM'].astype(str)

    return sorted_df

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some files.")

    # Add input file argument
    parser.add_argument('-i', '--input', type=str, required=True, help="Input file name")

    # Add hg19 fasta reference file
    parser.add_argument('-r', '--reference', type=str, required=True, help="hg19 reference fasta file name")

    # Add hg19 fasta reference file
    parser.add_argument('-r2', '--reference2', type=str, required=True, help="hg19 reference fai file name")

    # Add output file argument
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file name")

    # Add output file argument
    parser.add_argument('-o2', '--output2', type=str, required=True, help="Output file name for end of CNV")

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed input and output file names
    input_file_name = args.input
    reference_file_name = args.reference
    reference_file_name2 = args.reference2
    output_file_name = args.output
    output_file_end_name = args.output2
    date_str = date.today().strftime("%Y%m%d")

    chr_dic = {}
    with open(reference_file_name2, 'r') as f:
        for line in f:
            chr_name = line.split('\t')[0].strip()
            chr_len = int(line.split('\t')[1].strip())
            chr_dic[chr_name] = chr_len

    tree = ET.parse(input_file_name)
    root = tree.getroot()

    # Extract information from xml
    data = []

    for cna in root.findall('.//copy-number-alteration'):
        cna_data = extract_data_from_cna(cna)
        data.append(cna_data)
    df = pd.DataFrame(data)

    # Get chromosome order information
    chr_ordered = chr_pos_order(list(chr_dic))

    # Process the DataFrame
    df_processed,chrs = process_dataframe(df, chr_ordered, reference_file_name)

    # Generate sequential IDs which will be used later to match the CNV record with END position
    df_processed['ID'] = [f'r_{i+1}' for i in range(len(df_processed))]
    # print(df_processed)

    # used to generate END vcf
    df_end = df_processed.copy()

    # Apply lambda function to remove END=... part for the output vcf
    df_processed['INFO'] = df_processed['INFO'].apply(lambda x: re.sub(r'END=[^;]+;', '', x))

    # Extract the END values and update POS column for the END vcf
    df_end['POS'] = df_end['INFO'].apply(lambda x: re.search(r'END=([0-9]+);', x).group(1) if re.search(r'END=([0-9]+);', x) else df_end['POS'])

    # Remove END=... from INFO column for END vcf
    df_end['INFO'] = df_end['INFO'].apply(lambda x: re.sub(r'END=[^;]+;', '', x))

    # Update Ref with the chr and pos using hg19 for END vcf
    df_end['CHROM_POS'] = df_end.apply(lambda row: f"{row['#CHROM']}:{row['POS']}", axis=1)
    df_end['REF'] = df_end['CHROM_POS'].apply(lambda x: extract_nucleotide(reference_file_name, x))
    df_end.drop(columns=['CHROM_POS'], inplace=True)

    # Create VCF headers
    vcf_headers = create_vcf_header(date_str, chrs, chr_dic, input_file_name)

    # Write headers and data to VCF file
    with open(output_file_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        df_processed.to_csv(f, sep='\t', index=False)

    # Write headers and data to VCF file
    with open(output_file_end_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        df_end.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()
