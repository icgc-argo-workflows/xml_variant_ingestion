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
from Bio import SeqIO
from Bio.Seq import Seq
import requests

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
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=SOMATIC-OR-GERMLINE,Number=1,Type=String,Description="Somatic or Germline">',
        '##INFO=<ID=ZYGOSITY,Number=1,Type=String,Description="Sample Zygosity">',
        '##INFO=<ID=STATUS,Number=1,Type=String,Description="Status of Variant">',
        '##INFO=<ID=EQUIVOCAL,Number=0,Type=Flag,Description="Equivocal of Variant">',
        '##INFO=<ID=AO,Number=0,Type=Flag,Description="Analytical Only">'
    ])

    return headers

def extract_data_from_short_variant(short_variant):
    return {
        # Extract useful information from the xml file
        'AF': short_variant.get('allele-fraction'),
        'ALT': short_variant.get('alternate-sequence'),
        '#CHROM': short_variant.get('chromosome') if short_variant.get('chromosome') else short_variant.get('position').split(":")[0],
        'DP': short_variant.get('depth'),
        'POS': short_variant.get('genomic-start') if short_variant.get('genomic-start') else short_variant.get('position').split(":")[-1],
        'REF': short_variant.get('reference-sequence'),
        'SOMATIC-or-GERMLINE': short_variant.get('germline-status'),
        'ZYGOSITY': short_variant.get('tumor-zygosity'),
        'status': short_variant.get('status'),
        'equivocal': short_variant.get('equivocal'),
        'analytical-only': short_variant.get('analytical-only') if short_variant.get('analytical-only') else "false",
        'variant-type': short_variant.get('variant-type'),
        'cds-effect': short_variant.get('cds-effect'),
        'gene': short_variant.get('gene'),
    }

def update_missing_info_from_cds(short_variant,reference_file):
    ###Find matching reference seqeucne
    for record in SeqIO.parse(reference_file, "fasta"):
        if record.id==short_variant.get('#CHROM'):
            break

    ### Check for type of SNV
    if "del" in short_variant.get('cds-effect'):
        ###Populate initial variables
        start=int(short_variant.get("POS"))
        del_sequence=re.findall("[A-Z]+",short_variant.get("cds-effect"))[0]
        end=start+len(del_sequence)
        gene=short_variant.get('gene').replace("MLL2","KMT2D").replace("itfg3","FAM234A").replace("kiaa1377","CEP126").replace("fgf14","FGF14").replace("emsy","C11ORF30")

        ###Query for strand information from API
        url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
        response=requests.get(url)

        if response.status_code!=200:
            print("ERRORD1: Unable to ping https://www.genenetwork.nl/api/v1/gene/%s" % gene.lower())
            exit(1)
            
        strand=response.json()['gene']['strand']

        if strand==-1:
            ###If strand is reverse, find matching sequence and report as forward strand
            if record.seq[start:end].upper()==Seq(del_sequence).upper().reverse_complement():
                updated=start-1
                short_variant['REF']=record.seq[updated:end].upper()
                short_variant['ALT']=record.seq[updated:end-len(del_sequence)].upper()
            else:
                print("ErrorD2: Expected sequence does not match reverse complement of provided sequence: Original-%s vs Provided-%s" % (record.seq[start:end].upper(),Seq(del_sequence).upper().reverse_complement()))
        elif strand==1:
            ### If strand is forward, sanity check by finding matching sequence
            if record.seq[start:end].upper()==Seq(del_sequence).upper():
                updated=start-1
                short_variant['REF']=record.seq[updated:end].upper()
                short_variant['ALT']=record.seq[updated:end-len(del_sequence)].upper()
            else:
                print("ErrorD2: Expected sequence does not match provided sequence: Original-%s vs Provided-%s" % (record.seq[start:end].upper(),Seq(del_sequence).upper()))
        else:
            print("ErrorD4: RETURN STRAND DOES NOT COMPUTE")

        print("######",del_sequence,short_variant['REF'],short_variant['ALT'])

    elif "ins" in short_variant.get('cds-effect'):
        ###Populate initial variables
        start=int(short_variant.get("POS"))
        added_sequence=re.findall("[A-Z]+",short_variant.get("cds-effect"))[0]
        end=start+1
        gene=short_variant.get('gene').replace("MLL2","KMT2D").replace("itfg3","FAM234A").replace("kiaa1377","CEP126").replace("fgf14","FGF14").replace("emsy","C11ORF30")

        ###Query for strand information from API
        url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
        response=requests.get(url)

        if response.status_code!=200:
            print("ERRORI1: Unable to ping https://www.genenetwork.nl/api/v1/gene/%s" % gene.lower())
            exit(1)
            
        strand=response.json()['gene']['strand']

        if strand==-1:
            updated=start
            short_variant['REF']=record.seq[updated:end].upper()
            short_variant['ALT']=record.seq[updated:end].upper()+added_sequence
        elif strand==1:
            updated=start
            short_variant['REF']=record.seq[updated:end].upper()
            short_variant['ALT']=record.seq[updated:end].upper()+added_sequence
        else:
            print("ERRORI2: UNKNOWN STRAND FOUND IN %s" % (gene) )
            exit(1)
    elif ">" in short_variant.get('cds-effect'):
        ###Query for strand information from API
        gene=short_variant.get('gene').replace("MLL2","KMT2D").replace("itfg3","FAM234A").replace("kiaa1377","CEP126").replace("fgf14","FGF14").replace("emsy","C11ORF30")
        url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
        response=requests.get(url)

        if response.status_code!=200:
            print("ERRORD1: Unable to ping https://www.genenetwork.nl/api/v1/gene/%s" % gene.lower())
            exit(1)
            
        strand=response.json()['gene']['strand']

        if strand==1:
            short_variant['REF']=re.findall("[a-zA-Z]+",short_variant.get('cds-effect'))[0]
            short_variant['ALT']=re.findall("[a-zA-Z]+",short_variant.get('cds-effect'))[1]
        else:
            short_variant['REF']=Seq(re.findall("[a-zA-Z]+",short_variant.get('cds-effect'))[0]).reverse_complement()
            short_variant['ALT']=Seq(re.findall("[a-zA-Z]+",short_variant.get('cds-effect'))[1]).reverse_complement()
        ### Label single or multi based on length of affected REF sequence
        short_variant['variant-type']= 'multiple-nucleotide-substitution' if len(re.findall("[a-zA-Z]+",short_variant.get('cds-effect'))[0])>1 else 'single-nucleotide-substitution'
    else:
        print("ERROR: UNKNOWN VARIANT TYPE DETECTED" )


def split_MNV_SNV(df):
    # Create an empty list to store the expanded rows
    expanded_rows = []

    # Loop through each row
    for _, row in df.iterrows():
        ref = row['REF']
        alt = row['ALT']
        pos = int(row['POS'])

        # Ensure REF and ALT are the same length
        if len(ref) == len(alt):
            for i in range(len(ref)):
                new_row = row.copy()
                new_row['REF'] = ref[i]
                new_row['ALT'] = alt[i]
                new_row['POS'] = pos + i
                expanded_rows.append(new_row)

    # Create a new dataframe from the expanded rows
    df_expanded = pd.DataFrame(expanded_rows)

    # Reset the index
    df_expanded.reset_index(drop=True, inplace=True)

    return df_expanded

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

    df = df.copy()
    df['ID'] = '.'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = 'DP=' + df['DP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    ';SOMATIC-OR-GERMLINE=' + df['SOMATIC-or-GERMLINE'].fillna('NA') + \
                    ';ZYGOSITY=' + df['ZYGOSITY'].replace("not applicable","NA").fillna('NA') + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')
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
    parser.add_argument('-r', '--reference', type=str, required=True, help="hg19 reference file name")
    parser.add_argument('-r2', '--reference-fai', type=str, required=True, help="hg19 reference fai file name")

    # Add output file argument
    parser.add_argument('-o1', '--output1', type=str, required=True, help="Output file name for snv")
    parser.add_argument('-o2', '--output2', type=str, required=True, help="Output file name for indel")

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed input and output file names
    input_file_name = args.input
    reference_file = args.reference
    reference_fai_file_name = args.reference_fai
    output_file_snv = args.output1
    output_file_indel = args.output2
    date_str = date.today().strftime("%Y%m%d")

    tree = ET.parse(input_file_name)
    root = tree.getroot()

    chr_dic = {}
    with open(reference_fai_file_name, 'r') as f:
        for line in f:
            chr_name = line.split('\t')[0].strip()
            chr_len = int(line.split('\t')[1].strip())
            chr_dic[chr_name] = chr_len

    # Extract information from xml
    data = []
    total=root.findall('.//short-variant')
    for count,short_variant in enumerate(total):
        print("Processing Short Variants # %s/%s" % (str(count+1),len(total)))
        #print(short_variant)
        short_variant_data = extract_data_from_short_variant(short_variant)
        if short_variant_data['ALT'] is None and short_variant_data['REF'] is None and short_variant_data['cds-effect'] is not None and short_variant_data['gene'] is not None:
            update_missing_info_from_cds(short_variant_data,reference_file)
        data.append(short_variant_data)
        print(short_variant.get('REF'),short_variant.get('ALT'),short_variant.get('cds-effect'))
        print(short_variant_data)
    #print(data)

    df = pd.DataFrame(data)
    #print(df)
    df['SOMATIC-or-GERMLINE']=df['SOMATIC-or-GERMLINE'].replace("somatic","Somatic").replace("germline","Germline").replace("unknown","Unclassified")
    df_mnv = df[df['variant-type'] == 'multiple-nucleotide-substitution']
    df_snv = df[df['variant-type'] == 'single-nucleotide-substitution']
    df_indel = df[~df['variant-type'].isin(['multiple-nucleotide-substitution', 'single-nucleotide-substitution'])]

    # Get chromosome order information
    chr_ordered = chr_pos_order(list(chr_dic))

    # Split MNVs to multiple SNVs
    df_snvs = split_MNV_SNV(df_mnv)

    # Merge df_snvs with df_snv
    df_snv_total = pd.concat([df_snv, df_snvs]).reset_index(drop=True).query("REF!=ALT")

    # Process the DataFrame
    df_indel_processed,chrs_indel = process_dataframe(df_indel, chr_ordered)    # indel
    df_snv_total_processed,chrs_snv = process_dataframe(df_snv_total, chr_ordered)  # total snv

    # Create VCF headers
    vcf_headers_indel = create_vcf_header(date_str, chrs_indel, chr_dic, input_file_name)   # indel
    vcf_headers_snv_total = create_vcf_header(date_str, chrs_snv, chr_dic, input_file_name)

    # Write headers and data to VCF file
    if len(df_snv_total_processed)>0:
        with open(output_file_snv, 'w') as f:
            for header_line in vcf_headers_snv_total:
                f.write(f"{header_line}\n")
            df_snv_total_processed.to_csv(f, sep='\t', index=False)

    # Write headers and data to VCF file
    if len(df_indel_processed)>0:
        with open(output_file_indel, 'w') as f:
            for header_line in vcf_headers_indel:
                f.write(f"{header_line}\n")
            df_indel_processed.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()
