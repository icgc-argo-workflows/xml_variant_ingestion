#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

 Author: Junjun Zhang <junjun.zhang@oicr.on.ca>
         Linda Xiang <linda.xiang@oicr.on.ca>
 """

import os
import sys
import json
import sys
import argparse
import subprocess
import json
import re
import hashlib
import uuid
import tarfile
from datetime import date
import copy
from glob import glob
import yaml
import io
import shutil
import csv

workflow_full_name = {
    'variant-calling': 'Variant Calling'
}

def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()


def rename_file(f, payload, seq_experiment_analysis_dict, date_str):
    experimental_strategy = payload['experiment']['experimental_strategy'].lower()

    if f.endswith('.vcf.gz'):
        file_ext = 'vcf.gz'
    elif f.endswith('.vcf.gz.tbi'):
        file_ext = 'vcf.gz.tbi'
    else:
        sys.exit('Error: unknown aligned seq extention: %s' % f)

    variant_type = ''
    if 'short_variant' in f:
        variant_type = 'snv'
    elif 'rearrangement' in f:
        variant_type = 'sv'
    else:
        sys.exit('Error: unknown variant type: %s' % f)

    new_name = "%s.%s.%s.%s.%s.%s.%s" % (
        payload['studyId'],
        seq_experiment_analysis_dict['donor_id'],
        seq_experiment_analysis_dict['sample_id'],
        experimental_strategy,
        variant_type,
        date_str,
        file_ext
    )

    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    dst = os.path.join(os.getcwd(), new_dir, new_name)
    os.symlink(os.path.abspath(f), dst)

    return dst


def get_files_info(file_to_upload,pipeline_info):
    return {
        'fileName': os.path.basename(file_to_upload),
        'fileType': 'VCF' if file_to_upload.split(".")[-2] == 'vcf' else 'TBI',
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'dataType': 'Raw Variant Calls' if file_to_upload.split(".")[-2] == 'vcf' else 'VCF Index',
        'info': {
            'data_category': 'Simple Nucleotide Variation' if file_to_upload.split(".")[-4] == 'snv' or file_to_upload.split(".")[-5] == 'snv' else 'Rearrangement Variation',
            'analysis_tools': [{key.split(":")[-1]:pipeline_info[key]} for key in pipeline_info.keys()] # to work on it later
            }
    }

def main(args):
    with open(args.seq_experiment_analysis, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_experiment_analysis_dict = row

    pipeline_info = {}
    if args.pipeline_yml:
      with open(args.pipeline_yml, 'r') as f:
        pipeline_info = yaml.safe_load(f)

    tools_dict = dict(tool.split(' ', 1) for tool in seq_experiment_analysis_dict.get('analysis_tools (tools and versions)').split(', '))
    pipeline_info['FoundationOneCDx'] = tools_dict

    for key, value in pipeline_info.items():
        for sub_key, sub_value in value.items():
            value[sub_key] = str(sub_value)
        pipeline_info[key] = value

    payload = {
        'studyId': seq_experiment_analysis_dict.get('program_id'),
        'samples': [
            {
                'submitterSampleId' : seq_experiment_analysis_dict.get('submitter_sample_id'),
                'matchedNormalSubmitterSampleId' : seq_experiment_analysis_dict.get('submitter_matched_normal_sample_id') or None,
                'sampleType' : seq_experiment_analysis_dict.get('sample_type'),
                'specimen' : {
                    'submitterSpecimenId' : seq_experiment_analysis_dict.get('submitter_specimen_id'),
                    'tumourNormalDesignation' : seq_experiment_analysis_dict.get('tumour_normal_designation'),
                    'specimenTissueSource' : seq_experiment_analysis_dict.get('specimen_tissue_source'),
                    'specimenType' : seq_experiment_analysis_dict.get('specimen_type')
                },
                'donor' : {
                    'gender' : seq_experiment_analysis_dict.get('gender'),
                    'submitterDonorId': seq_experiment_analysis_dict.get('submitter_donor_id')
                }
            }
        ],
        'analysisType': { 'name': 'variant_calling' },
        'variant_calling_strategy': 'Tumour Only', # need further confirm
        'workflow': {
            'genome_build': 'GRCh38',
            'workflow_name': 'FoundationOneCDx',
            'workflow_version': '1',
            'workflow_short_name': 'F1CDx'
        },
        'experiment' : {
            'submitter_sequencing_experiment_id': seq_experiment_analysis_dict.get('submitter_sequencing_experiment_id') or "something",
            'platform' : seq_experiment_analysis_dict.get('platform'),
            'experimental_strategy': 'Targeted-Seq',
            'target_capture_kit': seq_experiment_analysis_dict.get('target_capture_kit'),
            'primary_target_regions': seq_experiment_analysis_dict.get('primary_target_regions'),
            'capture_target_regions': seq_experiment_analysis_dict.get('capture_target_regions'),
            'coverage': ['Coding Exons','Introns','Promoters']
        },
        'variant_class' : 'Somatic',  # need further confirm
        'files': []
    }

    # get file of the payload
    date_str = date.today().strftime("%Y%m%d")
    for f in args.files_to_upload:
        renamed_file = rename_file(f, payload, seq_experiment_analysis_dict, date_str)
        payload['files'].append(get_files_info(renamed_file,pipeline_info))

    with open("%s.variant_call.payload.json" % str(uuid.uuid4()), 'w') as f:
        f.write(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool: payload-variant-call')
    parser.add_argument("-f", "--files_to_upload", dest="files_to_upload", type=str, required=True,
                        nargs="+", help="Aligned reads files to upload")
    parser.add_argument("-a", "--seq_experiment_analysis", dest="seq_experiment_analysis", required=True,
                        help="Input analysis for sequencing experiment", type=str)
    parser.add_argument("-w", "--wf_name", dest="wf_name", required=True, help="Workflow name")
    parser.add_argument("-v", "--wf_version", dest="wf_version", required=True, help="Workflow version")
    parser.add_argument("-r", "--wf_run", dest="wf_run", required=True, help="workflow run ID")
    parser.add_argument("-s", "--wf_session", dest="wf_session", required=True, help="workflow session ID")
    parser.add_argument("-b", "--genome_build", dest="genome_build", default="GRCh38", help="Genome build")
    parser.add_argument("-p", "--pipeline_yml", dest="pipeline_yml", required=False, help="Pipeline info in yaml")

    args = parser.parse_args()

    main(args)
