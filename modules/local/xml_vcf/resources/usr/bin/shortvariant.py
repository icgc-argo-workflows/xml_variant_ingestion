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

def create_vcf_header(root, chrs, chr_dic):
    # Extract attributes from variant-report element
    curation_version_id = root.get('curation-version-id')
    disease = root.get('disease')
    flowcell_analysis = root.get('flowcell-analysis')
    gender = root.get('gender')
    human_genome_assembly = root.get('human-genome-assembly')
    specimen_id = root.get('specimen-id')
    study = root.get('study')
    test_request = root.get('test-request')
    test_type = root.get('test-type')

    # Extract attributes from sample element (assuming there is only one sample)
    sample = root.find('.//sample')
    name = sample.get('name')
    bait_set = sample.get('bait-set')

    # Construct the headers
    headers = [
        '##fileformat=VCFv4.3',
        f'##curation-version-id=\"{curation_version_id}\"',
        f'##disease="{disease}"',
        f'##flowcell-analysis="{flowcell_analysis}"',
        f'##gender="{gender}"',
        f'##human-genome-assembly="{human_genome_assembly}"',
        f'##specimen-id="{specimen_id}"',
        f'##study="{study}"',
        f'##test-request="{test_request}"',
        f'##test-type="{test_type}"',
        f'##name="{name}"',
        f'##bait-set="{bait_set}"',
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
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=SOMATIC-OR-GERMLINE,Number=1,Type=String,Description="somatic or germline">',
        '##INFO=<ID=ZYGOSITY,Number=1,Type=String,Description="sample zygosity">'
    ])

    return headers

def extract_data_from_short_variant(short_variant):
    return {
        # Extract useful information from the xml file
        'AF': short_variant.get('allele-fraction'),
        'ALT': short_variant.get('alternate-sequence'),
        '#CHROM': short_variant.get('chromosome'),
        'DP': short_variant.get('depth'),
        'POS': short_variant.get('genomic-start'),
        'REF': short_variant.get('reference-sequence'),
        'SOMATIC-or-GERMLINE': short_variant.get('germline-status'),
        'ZYGOSITY': short_variant.get('tumor-zygosity')
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


def chr_pos_sort(df, chr_ordered): # sort dataframe based on chr and pos

    cat_type = pd.CategoricalDtype(categories=chr_ordered, ordered=True)

    # Convert the #CHROM column to the categorical data type
    df['#CHROM'] = df['#CHROM'].astype(cat_type)

    # Sort the DataFrame first by '#CHROM' and then by 'POS'
    sorted_df = df.sort_values(by=['#CHROM', 'POS'])

    # Convert the chromosome numbers back to their original format
    sorted_df['#CHROM'] = sorted_df['#CHROM'].astype(str)

    return sorted_df

# Process dataframe by adding "ID", "QUAL", "FILTER" with default value "." and "INFO" with "DP=#;AF=#"
# where DP means depth and AF means allele frequency
def process_dataframe(df, chr_list):

    df['ID'] = '.'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = 'DP=' + df['DP'].astype(str) + ';AF=' + df['AF'].astype(str) + ';SOMATIC-OR-GERMLINE=' + df['SOMATIC-or-GERMLINE'] + ';ZYGOSITY=' + df['ZYGOSITY']
    df.drop(['DP', 'AF'], axis=1, inplace=True)
    df['POS'] = df['POS'].astype(int)
    desired_order = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df = df[desired_order]
    sorted_df = chr_pos_sort(df.copy(), chr_list)

    return sorted_df,sorted_df['#CHROM']

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some files.")

    # Add input file argument
    parser.add_argument('-i', '--input', type=str, required=True, help="Input file name")

    # Add hg19 fasta reference file
    parser.add_argument('-r', '--reference', type=str, required=True, help="hg19 reference fai file name")

    # Add output file argument
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file name")

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed input and output file names
    input_file_name = args.input
    reference_file_name = args.reference
    output_file_name = args.output

    tree = ET.parse(input_file_name)
    root = tree.getroot()

    chr_dic = {}
    with open(reference_file_name, 'r') as f:
        for line in f:
            chr_name = line.split('\t')[0].strip()
            chr_len = int(line.split('\t')[1].strip())
            chr_dic[chr_name] = chr_len

    # Extract information from xml
    data = []

    for short_variant in root.findall('.//short-variant'):
        short_variant_data = extract_data_from_short_variant(short_variant)
        data.append(short_variant_data)
    df = pd.DataFrame(data)

    # Get chromosome order information
    chr_ordered = chr_pos_order(list(chr_dic))

    # Process the DataFrame
    df_processed,chrs = process_dataframe(df, chr_ordered)

    # Create VCF headers
    vcf_headers = create_vcf_header(root, chrs, chr_dic)

    # Write headers and data to VCF file
    with open(output_file_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        df_processed.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()
