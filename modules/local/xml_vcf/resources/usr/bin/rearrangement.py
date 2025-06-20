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
from datetime import date
import requests
import string

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
        '##INFO=<ID=SRP,Number=1,Type=Integer,Description="Supporting Read Pairs">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural Variant">',
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of the Mate Breakend">',
        '##INFO=<ID=STATUS,Number=1,Type=String,Description="Status of Variant">',
        '##INFO=<ID=EQUIVOCAL,Number=0,Type=Flag,Description="Equivocal of Variant">',
        '##INFO=<ID=AO,Number=0,Type=Flag,Description="Analytical Only">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
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
    fasta = Fasta(fasta_file, as_raw=True)
    sequence = fasta[pos.split(':')[0]][int(pos.split(':')[1])-1]
    fasta.close()
    return sequence.upper()

# create alt
def alternative(pos1, pos2,dir,ref):
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

def generate_deletion(rearrangement_data,rearrangement,fasta_file,count):
    ###For whatever reason coordinates are reported as "chrN:N-N+1" or "chrN", if the prior....which set of coordinates to use?
    if "-" in rearrangement_data.get('Pos1'):
        ### Check if we have gene annotation and chimeric description.
        if rearrangement_data.get('targeted-gene') and rearrangement_data.get('other-gene'):
            if len(rearrangement)>0:
                description=rearrangement[0][0].get('description')
                ###Use chimeric description where if 3'-FGFR1(x18-3*)-5':FGFR1-upstream(20kB)
                if re.findall("3'|5'",description)[0]=="5'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                elif re.findall("3'|5'",description)[0]=="3'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]

                if "upstream" in description or "downstream" in description:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="5'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="3'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
            else:
                ### If we don't have annotations to work with query external API to determine strand
                gene=rearrangement.get('targeted-gene')
                url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                responseA=requests.get(url)

                if responseA.status_code!=200:
                    print("ERROR_DE1: Unable to ping %s" % url)
                    exit(1)
                rearrangement_data['pos1_strand']="+" if responseA.json()['gene']['strand']==1 else "-"

                ###Assign arbitary strand if partner gene is missing
                if rearrangement.get('other_gene')=='N/A':
                    rearrangement_data['pos2_strand']="+"
                else:
                    gene=rearrangement.get('other-gene')
                    url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                    responseB=requests.get(url)

                    if responseB.status_code!=200:
                        print("ERROR_DE2: Unable to ping %s" % url)
                        exit(1)
                    rearrangement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"
 
                ###Using strand info return 5' or 3' positions
                if re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]<re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]:
                    if rearrangement_data['pos1_strand']=="+":
                        rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                    else:
                        rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]


                    if rearrangement_data['pos2_strand']=="+":
                        rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                    else:
                        rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                else:
                    if rearrangement_data['pos1_strand']=="+":
                        rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]
                    else:
                        rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]


                    if rearrangement_data['pos2_strand']=="+":
                        rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
                    else:
                        rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]

    ###BECAUSE WE CANT EVEN EXPECT THEM TO PROVIDE COORDINATES WHERE POS1<POS2 CONSISTENTLY - SWAP THEM SO WE ALWAYS WRITE SMALLEST COORD FIRST
    if int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1]) > int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1]):
        tmpA=rearrangement_data['Pos1']
        tmpB=rearrangement_data['Pos2']
        rearrangement_data['Pos1']=tmpB
        rearrangement_data['Pos2']=tmpA

    data = [] # Create a list to store the extracted data
    junc = [] # contain junction
    dir = {} # contian direction info

    junc.append(rearrangement_data.get('Pos1'))
    junc.append(rearrangement_data.get('Pos2'))

    ref[rearrangement_data.get('Pos1')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos1'))
    ref[rearrangement_data.get('Pos2')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos2'))

    ####We're checking the partner breakpoint for additional context
    partner_breakpoints = rearrangement.findall('.//partner-breakpoint')
    if len(partner_breakpoints)>0:
        for partner_breakpoint in partner_breakpoints:
            chrom=partner_breakpoint.get('chromosome')
            position=int(partner_breakpoint.get('annotation-position'))
            pos1=int(rearrangement_data['Pos1'].split(":")[-1])
            pos2=int(rearrangement_data['Pos2'].split(":")[-1])
            ### Again things are weird where the annotation does not match the approximate position reported...this is the closest way I could attempt to match the two via distance
            if abs(position-pos1) < abs(position-pos2):
                dir[chrom + ":" + str(pos1)] = partner_breakpoint.get('genomic-disruption-direction')
            elif abs(position-pos2) < abs(position-pos1):
                dir[chrom + ":" + str(pos2)] = partner_breakpoint.get('genomic-disruption-direction')
            else : 
                ("PANIC nothing makes sense")
                exit(1)
    else:
        ### If we don't have breakpoint annotation, auto assign the strand
        dir[rearrangement_data['Pos1']] = rearrangement_data['pos1_strand']
        dir[rearrangement_data['Pos2']] = rearrangement_data['pos2_strand']

    data.append(rearrangement_data)

    sorted_junc = sorted(junc, key=custom_sort_key) # sort junction based on first chr, than pos

    # create a dictionary, naming bnd_# for all junctions based on sorted order
    junction = {}
    num = 0
    for item in sorted_junc:
        num +=1
        junction[item] = f'DEL_{count}'

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data)

    # Create a DataFrame for vcf file
    df_a = pd.DataFrame()


    df_a[['#CHROM', 'POS']] = df['Pos1'].str.split(':', expand=True)
    df_a['POS'] = df_a['POS'].astype(int)
    df_a['ID'] = df['Pos1'].map(junction)
    df_a['REF'] = df['Pos1'].map(ref)
    df_a['ALT'] = "<DEL>"
    df_a['QUAL'] = '.'
    df_a['FILTER'] = '.'
    df_a['INFO'] = 'SVTYPE=DEL' + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';END=' + (df['Pos2'].str.split(':', expand=True)[1].astype(int)).astype(str) + \
                    ";SVLEN=" + (-1*abs(df['Pos2'].str.split(':', expand=True)[1].astype(int)-df['Pos1'].str.split(':', expand=True)[1].astype(int))).astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')
    return(df_a)

def generate_duplication(rearrangement_data,rearrangement,fasta_file,count):
    ###For whatever reason coordinates are reported as "chrN:N-N+1" or "chrN", if the prior....which set of coordinates to use?
    if "-" in rearrangement_data.get('Pos1'):
        ### Check if we have gene annotation and chimeric description.
        if rearrangement_data.get('targeted-gene') and rearrangement_data.get('other-gene'):
            if len(rearrangement)>0:
                description=rearrangement[0][0].get('description')
                ###Use chimeric description where if 3'-FGFR1(x18-3*)-5':FGFR1-upstream(20kB)
                if re.findall("3'|5'",description)[1]=="5'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]
                elif re.findall("3'|5'",description)[1]=="3'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]

                if "upstream" in description or "downstream" in description:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="5'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="3'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
            else:
                ### If we don't have annotations to work with query external API to determine strand
                gene=rearrangement.get('targeted-gene')
                url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                responseA=requests.get(url)

                if responseA.status_code!=200:
                    print("ERROR_DU1: Unable to ping %s" % url)
                    exit(1)
                rearrangement_data['pos1_strand']="+" if responseA.json()['gene']['strand']==1 else "-"

                ###Assign arbitary strand if partner gene is missing
                if rearrangement.get('other_gene')=='N/A':
                    rearrangement_data['pos2_strand']="+"
                else:
                    gene=rearrangement.get('other-gene')
                    url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                    responseB=requests.get(url)

                    if responseB.status_code!=200:
                        print("ERROR_DU2: Unable to ping %s" % url)
                        exit(1)
                    rearrangement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"

                ###Using strand info return 5' or 3' positions
                if rearrangement_data['pos1_strand']=="+":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                else:
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]


                if rearrangement_data['pos2_strand']=="+":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                else:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]

    ###BECAUSE WE CANT EVEN EXPECT THEM TO PROVIDE COORDINATES WHERE POS1<POS2 CONSISTENTLY:
    if int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1]) > int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1]):
        tmpA=rearrangement_data['Pos1']
        tmpB=rearrangement_data['Pos2']
        rearrangement_data['Pos1']=tmpB
        rearrangement_data['Pos2']=tmpA

    data = [] # Create a list to store the extracted data
    junc = [] # contain junction
    dir = {} # contian direction info

    junc.append(rearrangement_data.get('Pos1'))
    junc.append(rearrangement_data.get('Pos2'))

    ref[rearrangement_data.get('Pos1')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos1'))
    ref[rearrangement_data.get('Pos2')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos2'))

    ###We're checking the partner breakpoint for additional context
    partner_breakpoints = rearrangement.findall('.//partner-breakpoint')
    if len(partner_breakpoints)>0:
        for partner_breakpoint in partner_breakpoints:
            chrom=partner_breakpoint.get('chromosome')
            position=int(partner_breakpoint.get('annotation-position'))
            pos1=int(rearrangement_data['Pos1'].split(":")[-1])
            pos2=int(rearrangement_data['Pos2'].split(":")[-1])
            ### Again things are weird where the annotation does not match the approximate position reported...this is the closest way I could attempt to match the two via distance
            if abs(position-pos1) < abs(position-pos2):
                dir[chrom + ":" + str(pos1)] = partner_breakpoint.get('genomic-disruption-direction')
            elif abs(position-pos2) < abs(position-pos1):
                dir[chrom + ":" + str(pos2)] = partner_breakpoint.get('genomic-disruption-direction')
            else : 
                ("PANIC nothing makes sense")
                exit(1)
    else:
        ### If we don't have breakpoint annotation, auto assign the strand
        dir[rearrangement_data['Pos1']] = rearrangement_data['pos1_strand']
        dir[rearrangement_data['Pos2']] = rearrangement_data['pos2_strand']

    data.append(rearrangement_data)

    sorted_junc = sorted(junc, key=custom_sort_key) # sort junction based on first chr, than pos

    # create a dictionary, naming bnd_# for all junctions based on sorted order
    junction = {}
    num = 0
    for item in sorted_junc:
        num +=1
        junction[item] = f'DUP_{count}'

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data)

    # Create a DataFrame for vcf file
    df_a = pd.DataFrame()


    df_a[['#CHROM', 'POS']] = df['Pos1'].str.split(':', expand=True)
    df_a['POS'] = df_a['POS'].astype(int)
    df_a['ID'] = df['Pos1'].map(junction)
    df_a['REF'] = df['Pos1'].map(ref)
    df_a['ALT'] = "<DUP>"
    df_a['QUAL'] = '.'
    df_a['FILTER'] = '.'
    df_a['INFO'] = 'SVTYPE=DUP' + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';END=' + (df['Pos2'].str.split(':', expand=True)[1].astype(int)).astype(str) + \
                    ";SVLEN=" + (abs(df['Pos2'].str.split(':', expand=True)[1].astype(int)-df['Pos1'].str.split(':', expand=True)[1].astype(int))).astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    return(df_a)

def generate_fusion(rearrangement_data,rearrangement,fasta_file,count):
    ###For whatever reason coordinates are reported as "chrN:N-N+1" or "chrN", if the prior....which set of coordinates to use?
    if "-" in rearrangement_data.get('Pos1'):
        if rearrangement_data.get('targeted-gene') and rearrangement_data.get('other-gene'):
            if len(rearrangement)>0:
                description=rearrangement[0][0].get('description')
                ###Use chimeric description where if 3'-FGFR1(x18-3*)-5':FGFR1-upstream(20kB)
                if re.findall("3'|5'",description)[1]=="5'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]
                elif re.findall("3'|5'",description)[1]=="3'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]

                if "upstream" in description or "downstream" in description:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="5'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="3'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
            else:
                ### If we don't have annotations to work with query external API to determine strand
                gene=rearrangement.get('targeted-gene')
                url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                responseA=requests.get(url)

                if responseA.status_code!=200:
                    print("ERROR_GF1: Unable to ping %s" % url)
                    exit(1)
                rearrangement_data['pos1_strand']="+" if responseA.json()['gene']['strand']==1 else "-"

                ###Assign arbitary strand if partner gene is missing
                if rearrangement.get('other_gene')=='N/A':
                    rearrangement_data['pos2_strand']="+"
                else:
                    gene=rearrangement.get('other-gene')
                    url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                    responseB=requests.get(url)

                    if responseB.status_code!=200:
                        print("ERROR_GF2: Unable to ping %s" % url)
                        exit(1)
                    rearrangement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"

                ###Using strand info return 5' or 3' positions
                if rearrangement_data['pos1_strand']=="+":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                else:
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]


                if rearrangement_data['pos2_strand']=="+":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                else:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]

    ###BECAUSE WE CANT EVEN EXPECT THEM TO PROVIDE COORDINATES WHERE POS1<POS2 CONSISTENTLY:
    if int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1]) > int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1]):
        tmpA=rearrangement_data['Pos1']
        tmpB=rearrangement_data['Pos2']
        rearrangement_data['Pos1']=tmpB
        rearrangement_data['Pos2']=tmpA

    data = [] # Create a list to store the extracted data
    junc = [] # contain junction
    dir = {} # contian direction info

    junc.append(rearrangement_data.get('Pos1'))
    junc.append(rearrangement_data.get('Pos2'))

    ref[rearrangement_data.get('Pos1')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos1'))
    ref[rearrangement_data.get('Pos2')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos2'))

    ####We're checking the partner breakpoint for additional context
    partner_breakpoints = rearrangement.findall('.//partner-breakpoint')
    if len(partner_breakpoints)>0:
        for partner_breakpoint in partner_breakpoints:
            chrom=partner_breakpoint.get('chromosome')
            position=int(partner_breakpoint.get('annotation-position'))
            pos1=int(rearrangement_data['Pos1'].split(":")[-1])
            pos2=int(rearrangement_data['Pos2'].split(":")[-1])
            ### Again things are weird where the annotation does not match the approximate position reported...this is the closest way I could attempt to match the two via distance
            if abs(position-pos1) < abs(position-pos2):
                dir[chrom + ":" + str(pos1)] = partner_breakpoint.get('genomic-disruption-direction')
            elif abs(position-pos2) < abs(position-pos1):
                dir[chrom + ":" + str(pos2)] = partner_breakpoint.get('genomic-disruption-direction')
            else : 
                ("PANIC nothing makes sense")
                exit(1)
    else:
        ### If we don't have breakpoint annotation, auto assign the strand
        dir[rearrangement_data['Pos1']] = rearrangement_data['pos1_strand']
        dir[rearrangement_data['Pos2']] = rearrangement_data['pos2_strand']

    data.append(rearrangement_data)

    sorted_junc = sorted(junc, key=custom_sort_key) # sort junction based on first chr, than pos

    # create a dictionary, naming bnd_# for all junctions based on sorted order
    junction = {}
    num = 0
    for item in sorted_junc:
        num +=1
        junction[item] = f'bnd_{count}'

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data)

    # Create a DataFrame for vcf file
    df_a = pd.DataFrame()

    df_a[['#CHROM', 'POS']] = df['Pos1'].str.split(':', expand=True)
    df_a['POS'] = df_a['POS'].astype(int)
    df_a['ID'] = df['Pos1'].map(junction)
    df_a['REF'] = df['Pos1'].map(ref)
    df_a['ALT'] = df.apply(lambda row: alternative(row['Pos1'], row['Pos2'],dir,ref), axis=1)
    df_a['QUAL'] = '.'
    df_a['FILTER'] = '.'
    df_a['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos2'].map(junction) + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    df_b = pd.DataFrame()

    df_b[['#CHROM', 'POS']] = df['Pos2'].str.split(':', expand=True)
    df_b['POS'] = df_b['POS'].astype(int)
    df_b['ID'] = df['Pos2'].map(junction)
    df_b['REF'] = df['Pos2'].map(ref)
    df_b['ALT'] = df.apply(lambda row: alternative(row['Pos2'], row['Pos1'],dir,ref), axis=1)
    df_b['QUAL'] = '.'
    df_b['FILTER'] = '.'
    df_b['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos1'].map(junction) + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    result = pd.concat([df_a, df_b], ignore_index=True)

    return(result)

def generate_inversion(rearrangement_data,rearrangement,fasta_file,count):
    ###To report breakpoints for an inversion, we ultimately need to report 4x breakends
    ###To do we'll make an identical set of coordinates
    ###For example original coordinates were 5'-EXON1-3':5'-EXON2-3':5'-EXON3-3'
    ###The inversion becomes 5'-EXON1-3':3'-EXON2-5':5'-EXON3-3'
    ###To report this we generate the following set of coordinates in VCF:
    ###coordinates ID mateID
    ###EXON1-3' breakendA  breakendD
    ###EXON1-3'+1 breakendB breakendC
    ###EXON3-3' breakendC breakendB
    ###EXON3-3'-1 breakendD breakendA

    ###We generate a complement
    complement_data=rearrangement_data.copy()
    ###For whatever reason coordinates are reported as "chrN:N-N+1" or "chrN", if the prior....which set of coordinates to use?
    if "-" in rearrangement_data.get('Pos1'):
        if rearrangement_data.get('targeted-gene') and rearrangement_data.get('other-gene'):
            if len(rearrangement)>0:
                description=rearrangement[0][0].get('description')
                ###Use chimeric description where if 3'-FGFR1(x18-3*)-5':FGFR1-upstream(20kB)
                if re.findall("3'|5'",description)[1]=="5'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]
                    complement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1])-1)
                elif re.findall("3'|5'",description)[1]=="3'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                    complement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos1'])[-1])+1)

                if re.findall("3'|5'",description)[2]=="5'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                    complement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1])-1)
                elif re.findall("3'|5'",description)[2]=="3'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
                    complement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos2'])[-1])+1)
            else:
                ### If we don't have annotations to work with query external API to determine strand
                gene=rearrangement.get('targeted-gene')
                url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                responseA=requests.get(url)

                if responseA.status_code!=200:
                    print("ERROR_GI1: Unable to ping %s" % url)
                    exit(1)
                rearrangement_data['pos1_strand']="+" if responseA.json()['gene']['strand']==1 else "-"
                complement_data['pos1_strand']="+" if responseA.json()['gene']['strand']==1 else "-"

                ###Assign arbitary strand if partner gene is missing
                if rearrangement.get('other_gene')=='N/A':
                    rearrangement_data['pos2_strand']="+"
                    complement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"
                else:
                    gene=rearrangement.get('other-gene')
                    url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                    responseB=requests.get(url)

                    if responseB.status_code!=200:
                        print("ERROR_GI2: Unable to ping %s" % url)
                        exit(1)
                    rearrangement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"
                    complement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"

                ###Using strand info return 5' or 3' positions
                if rearrangement_data['pos1_strand']=="+":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                    complement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos1'])[-1])+1)
                else:
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]
                    complement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1])-1)


                if rearrangement_data['pos2_strand']=="+":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                    complement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1])-1)
                else:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
                    complement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+str(int(re.findall("[0-9]+",rearrangement_data['Pos2'])[-1])+1)

    ###BECAUSE WE CANT EVEN EXPECT THEM TO PROVIDE COORDINATES WHERE POS1<POS2 CONSISTENTLY:
    if int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1]) > int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1]):
        tmpA=rearrangement_data['Pos1']
        tmpB=rearrangement_data['Pos2']
        rearrangement_data['Pos1']=tmpB
        rearrangement_data['Pos2']=tmpA
        tmpA=complement_data['Pos1']
        tmpB=complement_data['Pos2']
        complement_data['Pos1']=tmpB
        complement_data['Pos2']=tmpA


    ###Swap coordinate pairs
    tmpA=rearrangement_data['Pos2']
    tmpB=complement_data['Pos2']
    rearrangement_data['Pos2']=tmpB
    complement_data['Pos2']=tmpA

    data = [] # Create a list to store the extracted data
    junc = [] # contain junction
    dir = {} # contian direction info

    junc.append(rearrangement_data.get('Pos1'))
    junc.append(rearrangement_data.get('Pos2'))
    junc.append(complement_data.get('Pos1'))
    junc.append(complement_data.get('Pos2'))

    ref[rearrangement_data.get('Pos1')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos1'))
    ref[rearrangement_data.get('Pos2')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos2'))
    ref[complement_data.get('Pos1')] = extract_nucleotide(fasta_file, complement_data.get('Pos1'))
    ref[complement_data.get('Pos2')] = extract_nucleotide(fasta_file, complement_data.get('Pos2'))

    ####We're checking the partner breakpoint for additional context
    partner_breakpoints = rearrangement.findall('.//partner-breakpoint')
    if len(partner_breakpoints)>0:
        for partner_breakpoint in partner_breakpoints:
                chrom=partner_breakpoint.get('chromosome')
                position=int(partner_breakpoint.get('annotation-position'))
                pos1=int(rearrangement_data['Pos1'].split(":")[-1])
                pos2=int(rearrangement_data['Pos2'].split(":")[-1])
                ### Again things are weird where the annotation does not match the approximate position reported...this is the closest way I could attempt to match the two via distance
                if abs(position-pos1) < abs(position-pos2):
                    dir[chrom + ":" + str(pos1)] = partner_breakpoint.get('genomic-disruption-direction')
                    dir[complement_data["Pos1"]]=dir[chrom + ":" + str(pos1)]
                elif abs(position-pos2) < abs(position-pos1):
                    dir[chrom + ":" + str(pos2)] = partner_breakpoint.get('genomic-disruption-direction')
                    dir[complement_data["Pos2"]]=dir[chrom + ":" + str(pos2)]
                else : 
                    ("PANIC nothing makes sense")
                    exit(1)
    else:
        ### If we don't have breakpoint annotation, auto assign the strand
        dir[rearrangement_data['Pos1']]=rearrangement_data['pos1_strand']
        dir[rearrangement_data['Pos2']]=rearrangement_data['pos2_strand']
        dir[complement_data['Pos1']]=rearrangement_data['pos1_strand']
        dir[complement_data['Pos2']]=rearrangement_data['pos2_strand']

                
    data.append(rearrangement_data)
    data.append(complement_data)

    sorted_junc = sorted(junc, key=custom_sort_key) # sort junction based on first chr, than pos

    # create a dictionary, naming bnd_# for all junctions based on sorted order
    junction = {}
    num = 0
    for sec_count,item in enumerate(sorted_junc):
        num +=1
        junction[item] = f'bnd_{count}{string.ascii_uppercase[sec_count]}'

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data)

    # Create a DataFrame for vcf file
    df_a = pd.DataFrame()

    df_a[['#CHROM', 'POS']] = df['Pos1'].str.split(':', expand=True)
    df_a['POS'] = df_a['POS'].astype(int)
    df_a['ID'] = df['Pos1'].map(junction)
    df_a['REF'] = df['Pos1'].map(ref)
    df_a['ALT'] = df.apply(lambda row: alternative(row['Pos1'], row['Pos2'],dir,ref), axis=1)
    df_a['QUAL'] = '.'
    df_a['FILTER'] = '.'
    df_a['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos2'].map(junction) + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    df_b = pd.DataFrame()

    df_b[['#CHROM', 'POS']] = df['Pos2'].str.split(':', expand=True)
    df_b['POS'] = df_b['POS'].astype(int)
    df_b['ID'] = df['Pos2'].map(junction)
    df_b['REF'] = df['Pos2'].map(ref)
    df_b['ALT'] = df.apply(lambda row: alternative(row['Pos2'], row['Pos1'],dir,ref), axis=1)
    df_b['QUAL'] = '.'
    df_b['FILTER'] = '.'
    df_b['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos1'].map(junction) + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    result = pd.concat([df_a, df_b], ignore_index=True)

    return(result)

def generate_rearrangement(rearrangement_data,rearrangement,fasta_file,count):
    ###For whatever reason coordinates are reported as "chrN:N-N+1" or "chrN", if the prior....which set of coordinates to use?
    if "-" in rearrangement_data.get('Pos1'):
        if rearrangement_data.get('targeted-gene') and rearrangement_data.get('other-gene'):
            if len(rearrangement)>0:
                description=rearrangement[0][0].get('description')
                ###Use chimeric description where if 3'-FGFR1(x18-3*)-5':FGFR1-upstream(20kB)
                if re.findall("3'|5'",description)[1]=="5'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]
                elif re.findall("3'|5'",description)[1]=="3'":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]

                if "upstream" in description or "downstream" in description:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="5'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                elif re.findall("3'|5'",description)[2]=="3'":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]
            else:
                ### If we don't have annotations to work with query external API to determine strand
                gene=rearrangement.get('targeted-gene')
                url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                responseA=requests.get(url)

                if responseA.status_code!=200:
                    print("ERROR_GR1: Unable to ping %s" % url)
                    exit(1)
                rearrangement_data['pos1_strand']="+" if responseA.json()['gene']['strand']==1 else "-"

                ###Assign arbitary strand if partner gene is missing
                gene=rearrangement.get('other-gene')
                if gene.lower()=='n/a':
                    rearrangement_data['pos2_strand']="+"
                else:
                    url='https://www.genenetwork.nl/api/v1/gene/%s' % gene.lower()
                    responseB=requests.get(url)

                    if responseB.status_code!=200:
                        print("ERROR_GR2: Unable to ping %s" % url)
                        exit(1)
                    rearrangement_data['pos2_strand']="+" if responseB.json()['gene']['strand']==1 else "-"

                ###Using strand info return 5' or 3' positions
                if rearrangement_data['pos1_strand']=="+":
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[-1]
                else:
                    rearrangement_data['Pos1']=re.findall("^chr[0-9]+",rearrangement_data['Pos1'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos1'])[1]


                if rearrangement_data['pos2_strand']=="+":
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[1]
                else:
                    rearrangement_data['Pos2']=re.findall("^chr[0-9]+",rearrangement_data['Pos2'])[0]+":"+re.findall("[0-9]+",rearrangement_data['Pos2'])[-1]

    ###BECAUSE WE CANT EVEN EXPECT THEM TO PROVIDE COORDINATES WHERE POS1<POS2 CONSISTENTLY:
    if int(re.findall("[0-9]+",rearrangement_data['Pos1'])[1]) > int(re.findall("[0-9]+",rearrangement_data['Pos2'])[1]):
        tmpA=rearrangement_data['Pos1']
        tmpB=rearrangement_data['Pos2']
        rearrangement_data['Pos1']=tmpB
        rearrangement_data['Pos2']=tmpA


    data = [] # Create a list to store the extracted data
    junc = [] # contain junction
    dir = {} # contian direction info

    junc.append(rearrangement_data.get('Pos1'))
    junc.append(rearrangement_data.get('Pos2'))

    ref[rearrangement_data.get('Pos1')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos1'))
    ref[rearrangement_data.get('Pos2')] = extract_nucleotide(fasta_file, rearrangement_data.get('Pos2'))

    ####We're checking the partner breakpoint for additional context
    partner_breakpoints = rearrangement.findall('.//partner-breakpoint')
    if len(partner_breakpoints)>0:
        for partner_breakpoint in partner_breakpoints:
            chrom=partner_breakpoint.get('chromosome')
            position=int(partner_breakpoint.get('annotation-position'))
            pos1=int(rearrangement_data['Pos1'].split(":")[-1])
            pos2=int(rearrangement_data['Pos2'].split(":")[-1])
            ### Again things are weird where the annotation does not match the approximate position reported...this is the closest way I could attempt to match the two via distance
            if abs(position-pos1) < abs(position-pos2):
                dir[chrom + ":" + str(pos1)] = partner_breakpoint.get('genomic-disruption-direction')
            elif abs(position-pos2) < abs(position-pos1):
                dir[chrom + ":" + str(pos2)] = partner_breakpoint.get('genomic-disruption-direction')
            else : 
                ("PANIC nothing makes sense")
                exit(1)
    else:
        dir[rearrangement_data['Pos1']] = rearrangement_data['pos1_strand']
        dir[rearrangement_data['Pos2']] = rearrangement_data['pos2_strand']

    data.append(rearrangement_data)

    sorted_junc = sorted(junc, key=custom_sort_key) # sort junction based on first chr, than pos

    # create a dictionary, naming bnd_# for all junctions based on sorted order
    junction = {}
    num = 0
    for item in sorted_junc:
        num +=1
        junction[item] = f'bnd_{count}'

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data)

    # Create a DataFrame for vcf file
    df_a = pd.DataFrame()

    df_a[['#CHROM', 'POS']] = df['Pos1'].str.split(':', expand=True)
    df_a['POS'] = df_a['POS'].astype(int)
    df_a['ID'] = df['Pos1'].map(junction)
    df_a['REF'] = df['Pos1'].map(ref)
    df_a['ALT'] = df.apply(lambda row: alternative(row['Pos1'], row['Pos2'],dir,ref), axis=1)
    df_a['QUAL'] = '.'
    df_a['FILTER'] = '.'
    df_a['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos2'].map(junction) + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    ';STATUS=' + df['status'] + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    df_b = pd.DataFrame()

    df_b[['#CHROM', 'POS']] = df['Pos2'].str.split(':', expand=True)
    df_b['POS'] = df_b['POS'].astype(int)
    df_b['ID'] = df['Pos2'].map(junction)
    df_b['REF'] = df['Pos2'].map(ref)
    df_b['ALT'] = df.apply(lambda row: alternative(row['Pos2'], row['Pos1'],dir,ref), axis=1)
    df_b['QUAL'] = '.'
    df_b['FILTER'] = '.'
    df_b['INFO'] = 'SVTYPE=BND;MATEID=' + df['Pos1'].map(junction) + \
                    ';SRP=' + df['SRP'].astype(str) + \
                    ';AF=' + df['AF'].astype(str) + \
                    df['equivocal'].apply(lambda x: ';EQUIVOCAL' if x.lower() == 'true' else '') + \
                    df['analytical-only'].apply(lambda x: ';AO' if x.lower() == 'true' else '')

    result = pd.concat([df_a, df_b], ignore_index=True)
    return(result)

# define these two so they can be used both in main and alternative function
ref = {} # contian reference nt

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
    date_str = date.today().strftime("%Y%m%d")

    # Parse the XML file
    tree = ET.parse(input_file_name)
    root = tree.getroot()
    fasta_file = reference_file_name
    fai_file = reference_file_name2

    data = [] # Create a list to store the extracted data
    junc = [] # contain junction

    # Iterate over each "short-variant"
    variant_data=[]
    total=root.findall('.//rearrangement')
    for count,rearrangement in enumerate(root.findall('.//rearrangement')):
        print("Processing Rearrangements # %s/%s" % (str(count+1),len(total)))
        # Extract data from "short-variant"
        rearrangement_data = {
            'AF': rearrangement.get('allele-fraction'),
            'Pos1': rearrangement.get('pos1'),
            'Pos2': rearrangement.get('pos2'),
            'SRP': rearrangement.get('supporting-read-pairs'),
            'type': rearrangement.get('type'),
            'status': rearrangement.get('status'),
            'equivocal': rearrangement.get('equivocal'),
            'analytical-only': rearrangement.get('analytical-only') if rearrangement.get('analytical-only') else "false",
            'genomic-type': rearrangement.find('.//chimeric-junctions').get('genomic-type') if rearrangement.find('.//chimeric-junctions') else None, #unused ignore
            'chimeric-description': rearrangement.find('.//chimeric-junction').get('description') if rearrangement.find('.//chimeric-junctions') else None, #unused ignore
            'rearrangement-description' : rearrangement.get('description'),
            'targeted-gene' : rearrangement.get('targeted-gene'),
            'other-gene' : rearrangement.get('other-gene')
        }
        
        ###according to description and/or rearrangement type
        if "deletion" in rearrangement_data['rearrangement-description'] or "deletion" in rearrangement_data['type']:
            variant_data.append(generate_deletion(rearrangement_data,rearrangement,fasta_file,count+1))
        elif "duplication" in rearrangement_data['rearrangement-description'] or "duplication" in rearrangement_data['type']:
            variant_data.append(generate_duplication(rearrangement_data,rearrangement,fasta_file,count+1))
        elif "fusion" in rearrangement_data['rearrangement-description'] or "fusion" in rearrangement_data['type']:
            variant_data.append(generate_fusion(rearrangement_data,rearrangement,fasta_file,count+1))
        elif "inversion" in rearrangement_data['rearrangement-description'] or "inversion" in rearrangement_data['type']:
            variant_data.append(generate_inversion(rearrangement_data,rearrangement,fasta_file,count+1))
        else:
            variant_data.append(generate_rearrangement(rearrangement_data,rearrangement,fasta_file,count+1))

    result=pd.concat(variant_data)
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
    vcf_headers = create_vcf_header(date_str, chrs, chr_dic, input_file_name)

    # Write headers and data to VCF file
    with open(output_file_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        sorted_df.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()