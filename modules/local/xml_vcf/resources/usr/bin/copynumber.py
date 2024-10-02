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
        '##fileformat=VCFv4.3',
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
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural Variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of Copy Number Alteration">',
        '##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy Number">',
        '##INFO=<ID=RATIO,Number=1,Type=Float,Description="Copy Number Ratio">',
        '##INFO=<ID=STATUS,Number=1,Type=String,Description="Status of Variant">',
        '##INFO=<ID=EQUIVOCAL,Number=0,Type=Flag,Description="Equivocal of Variant">',
        '##INFO=<ID=AO,Number=0,Type=Flag,Description="Analytical Only">'
    ])

    return headers

def extract_data_from_cna(cna):
    return {
        # Extract useful information from the xml file
        'copy-number': cna.get('copy-number'),
        'position': cna.get('position'),
        'ratio': cna.get('ratio'),
        'type': cna.get('type'),
        'status': cna.get('status'),
        'equivocal': cna.get('equivocal'),
        'analytical-only': cna.get('analytical-only')
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

# extract ref nucleotide
def extract_nucleotide(fasta_file, pos):
    fasta = Fasta(fasta_file, as_raw=True)
    sequence = fasta[pos.split(':')[0]][int(pos.split(':')[1])-2].upper() # extract the position before this position's nucleotide, so -2
    fasta.close()
    return sequence

# Function to calculate the interval length
def calculate_length(position):
    # Split the position to get start and end as integers
    chrom, positions = position.split(':')
    start, end = map(int, positions.split('-'))
    return end - start + 1

# Process dataframe by adding "ID", "QUAL", "FILTER" with default value "." and "INFO"
def process_dataframe(df, chr_list, fasta_file):

    df = df.copy()

    df['#CHROM'] = df['position'].apply(lambda x: x.split(':')[0])
    df['POS'] = df['position'].apply(lambda x: x.split(':')[1].split('-')[0]).astype(int) - 1
    df['ID'] = '.'
    df['ref_pos'] = df['position'].apply(lambda x: x.split('-')[0])
    df['REF'] = df['ref_pos'].apply(lambda x: extract_nucleotide(fasta_file, x))
    df['ALT'] = '<CNV>'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = 'SVTYPE=CNV;SVLEN=' + df['position'].apply(lambda x: str(calculate_length(x))) + \
                    ';TYPE=' + df['type'] + \
                    ';CN=' + df['copy-number'] + \
                    ';RATIO=' + df['ratio'] + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

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

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed input and output file names
    input_file_name = args.input
    reference_file_name = args.reference
    reference_file_name2 = args.reference2
    output_file_name = args.output
    # output_file_end_name = args.output2
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
    print(df)

    # Get chromosome order information
    chr_ordered = chr_pos_order(list(chr_dic))

    # Process the DataFrame
    df_processed,chrs = process_dataframe(df, chr_ordered, reference_file_name)

    # Create VCF headers
    vcf_headers = create_vcf_header(date_str, chrs, chr_dic, input_file_name)

    # Write headers and data to VCF file
    with open(output_file_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        df_processed.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()
