#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2022,  icgc-argo

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

import argparse
import gzip

def main(main_vcf, end_vcf, output_vcf):
    # Create a dictionary from the end VCF file
    id_pos_dict = {}
    with gzip.open(end_vcf, 'rt') as end:
        for line in end:
            if not line.startswith('#'):  # Skip header lines
                fields = line.strip().split('\t')
                pos = fields[1]
                ids = fields[2]
                id_pos_dict[ids] = pos

    # Read the main VCF file and write the modified VCF file
    with gzip.open(main_vcf, 'rt') as main, open(output_vcf, 'w') as vcf_out:
        for line in main:
            if line.startswith('#'):  # Copy header lines directly
                vcf_out.write(line)
            else:
                fields = line.strip().split('\t')
                record_id = fields[2]

                # Modify INFO field if ID exists in the dictionary
                if record_id in id_pos_dict:
                    info_field = fields[7]
                    # Add or update the 'END_POS' key in the INFO field
                    info_field = f"END_POS={id_pos_dict[record_id]};{info_field}"
                    fields[7] = info_field

                # Optionally, reset the ID to '.'
                fields[2] = '.'

                # Write the modified record to the output VCF
                vcf_out.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify VCF files by adding END_POS field.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the main input VCF file (gzipped).")
    parser.add_argument('-i2', '--input2', type=str, required=True, help="Path to the end VCF file (gzipped).")
    parser.add_argument('-o', '--output', type=str, required=True, help="Path to the output VCF file (gzipped).")

    args = parser.parse_args()

    main(args.input, args.input2, args.output)
