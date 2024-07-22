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
from pyfaidx import Fasta
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
        f'##bait-set="{bait_set}"'
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
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">'
    ])

    return headers

def custom_sort_key(item):
    parts = item.split(':')
    try:
        chr_num = int(parts[0].replace('chr', ''))
        coord = int(parts[1])
        return (chr_num, coord)
    except (ValueError, IndexError):
        return item

# extract ref nucletide
def extract_nucleotide(fasta_file, pos):
    fasta = Fasta(fasta_file)
    sequence = fasta[pos.split(':')[0]][int(pos.split(':')[1])-1].seq
    fasta.close()
    return sequence

# create alt
def alternative(pos1, pos2):
    dir1 = dir[pos1]
    dir2 = dir[pos2]
    ref_nt = ref[pos1]
    if dir1 == '+':
        if dir2 == '+':
            return ref_nt + ']' + pos2 + ']'
        elif dir2 == '-':
            return ref_nt + '[' + pos2 + '['
    elif dir1 == '-':
        if dir2 == '+':
            return ']' + pos2 + ']' + ref_nt
        elif dir2 == '-':
            return '[' + pos2 + '[' + ref_nt
    return None

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

# Sort the final output with chr and pos
def chr_pos_sort(df, chr_ordered):

    cat_type = pd.CategoricalDtype(categories=chr_ordered, ordered=True)

    # Convert the #CHROM column to the categorical data type
    df['#CHROM'] = df['#CHROM'].astype(cat_type)
    # df['POS'] = df['POS'].astype(int)

    # Sort the DataFrame first by '#CHROM' and then by 'POS'
    sorted_df = df.sort_values(by=['#CHROM', 'POS'])

    # Convert the chromosome numbers back to their original format
    sorted_df['#CHROM'] = sorted_df['#CHROM'].astype(str)

    return sorted_df,sorted_df['#CHROM']

# define these two so they can be used both in main and alternative function
ref = {} # contian reference nt
dir = {} # contian direction info

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

    # Parse the XML file
    tree = ET.parse(input_file_name)
    root = tree.getroot()
    fasta_file = reference_file_name
    fai_file = reference_file_name2

    data = [] # Create a list to store the extracted data
    junc = [] # contain junction

    # Iterate over each "short-variant"
    for rearrangement in root.findall('.//rearrangement'):
        # Extract data from "short-variant"
        rearrangement_data = {
            'AF': rearrangement.get('allele-fraction'),
            'Pos1': rearrangement.get('pos1'),
            'Pos2': rearrangement.get('pos2'),
            'DP': rearrangement.get('supporting-read-pairs'),
            'type': rearrangement.get('type'),
            'genomic-type': rearrangement.find('.//chimeric-junctions').get('genomic-type'),
            'description': rearrangement.find('.//chimeric-junction').get('description'),
        }
        junc.append(rearrangement.get('pos1'))
        junc.append(rearrangement.get('pos2'))

        ref[rearrangement.get('pos1')] = extract_nucleotide(fasta_file, rearrangement.get('pos1'))
        ref[rearrangement.get('pos2')] = extract_nucleotide(fasta_file, rearrangement.get('pos2'))
        print(rearrangement.get('pos1'), ref[rearrangement.get('pos1')])
        # Check if there are any "short-variant-intragenic-annotation" elements

        partner_breakpoints = rearrangement.findall('.//partner-breakpoint')
        for partner_breakpoint in partner_breakpoints:
            dir[partner_breakpoint.get('chromosome') + ":" + partner_breakpoint.get('annotation-position')] = partner_breakpoint.get('genomic-disruption-direction')
        #    print(partner_breakpoint.get('chromosome') + ":" + partner_breakpoint.get('annotation-position'), dir[partner_breakpoint.get('chromosome') + ":" + partner_breakpoint.get('annotation-position')])

        data.append(rearrangement_data)

    sorted_junc = sorted(junc, key=custom_sort_key) # sort junction based on first chr, than pos

    # create a dictionary, naming bnd_# for all junctions based on sorted order
    junction = {}
    num = 0
    for item in sorted_junc:
        num +=1
        junction[item] = f'bnd_{num}'

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data)

    # Create a DataFrame for vcf file
    df_a = pd.DataFrame()

    df_a[['#CHROM', 'POS']] = df['Pos1'].str.split(':', expand=True)
    df_a['POS'] = df_a['POS'].astype(int)
    df_a['ID'] = df['Pos1'].map(junction)
    df_a['REF'] = df['Pos1'].map(ref)
    df_a['ALT'] = df.apply(lambda row: alternative(row['Pos1'], row['Pos2']), axis=1)
    df_a['QUAL'] = '.'
    df_a['FILTER'] = '.'
    df_a['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos2'].map(junction) + ';DP=' + df['DP'].astype(str) + ';AF=' + df['AF'].astype(str)

    df_b = pd.DataFrame()

    df_b[['#CHROM', 'POS']] = df['Pos2'].str.split(':', expand=True)
    df_b['POS'] = df_b['POS'].astype(int)
    df_b['ID'] = df['Pos2'].map(junction)
    df_b['REF'] = df['Pos2'].map(ref)
    df_b['ALT'] = df.apply(lambda row: alternative(row['Pos2'], row['Pos1']), axis=1)
    df_b['QUAL'] = '.'
    df_b['FILTER'] = '.'
    df_b['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos1'].map(junction) + ';DP=' + df['DP'].astype(str) + ';AF=' + df['AF'].astype(str)

    result = pd.concat([df_a, df_b], ignore_index=True)

    chr_dic = {}
    with open(fai_file, 'r') as f:
        for line in f:
            chr_name = line.split('\t')[0].strip()
            chr_len = int(line.split('\t')[1].strip())
            chr_dic[chr_name] = chr_len

    # Get chromosome order information
    chr_ordered = chr_pos_order(list(chr_dic))

    # Process the DataFrame
    sorted_df,chrs = chr_pos_sort(result.copy(), chr_ordered)

    # Create VCF headers
    vcf_headers = create_vcf_header(root, chrs, chr_dic)

    # Write headers and data to VCF file
    with open(output_file_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        sorted_df.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()
