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

 Author:Guanqiao Feng <gfeng@oicr.on.ca>
        Junjun Zhang <junjun.zhang@oicr.on.ca>
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

    if f.endswith('xml'):
        file_ext = 'xml'
    else:
        sys.exit('Error: unknown aligned seq extention: %s' % f)

    new_name = "%s.%s.%s.%s.%s.%s.%s.%s.%s" % (
        payload['studyId'],
        seq_experiment_analysis_dict['donor_id'],
        seq_experiment_analysis_dict['sample_id'],
        experimental_strategy,
        date_str,
        payload['workflow']['workflow_short_name'],
        seq_experiment_analysis_dict['variant_class'].replace(",", "-").lower(),
        variant_type,
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


def get_files_info(file_to_upload, date_str,seq_experiment_analysis_dict,donor_id,sample_id):
    print(seq_experiment_analysis_dict)
    print(file_to_upload,donor_id,sample_id)

    file_info = {
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'dataType' : "Original submitted XML document",
        'info': {
            'data_category': 'Supplement',
            'data_subtypes': None,
            'files_in_tgz': []
        }
    }
    # Name format for index file: MONSTAR-JP.DO264836.SA626657.targeted-seq.20240926.somatic-germline.cnv.vcf.gz.tbi
    # Name format for vcf file: MONSTAR-JP.DO264836.SA626657.targeted-seq.20240926.somatic-germline.cnv.vcf.gz
    if file_to_upload.endswith('xml'):
        data_type = "Supplement"
        # Check for variations in the third-last part (e.g., .snv., .indel., etc.)
        if '.snv.' in file_to_upload:
            data_category = 'Simple Nucleotide Variation'
        elif '.indel.' in file_to_upload:
            data_category = 'Simple Nucleotide Variation'
        elif '.sv.' in file_to_upload:
            data_category = 'Structural Variation'
        elif '.cnv.' in file_to_upload:
            data_category = 'Copy Number Variation'
        else:
            raise ValueError(f"Data type not recognized for file: {file_to_upload}")
    print([
    seq_experiment_analysis_dict.get('program_id'),
    donor_id,
    sample_id,
    seq_experiment_analysis_dict['experimental_strategy'],
    date_str,
    'tgz']
    )
    new_fname = '.'.join([
    seq_experiment_analysis_dict.get('program_id'),
    donor_id,
    sample_id,
    seq_experiment_analysis_dict['experimental_strategy'],
    date_str,
    'tgz'
    ])

    file_info['fileName'] = new_fname
    file_info['fileType'] = new_fname.split('.')[-1].upper()

    with tarfile.open(file_to_upload, 'r') as tar:
      for member in tar.getmembers():
        file_info['info']['files_in_tgz'].append(member.name)

    new_dir = 'out'
    try:
      os.mkdir(new_dir)
    except FileExistsError:
      pass

    dst = os.path.join(os.getcwd(), new_dir, new_fname)
    os.symlink(os.path.abspath(file_to_upload), dst)


    return file_info

def prepare_tarball(sampleId, qc_files, tool_list):

    tgz_dir = 'tarball'
    try:
      os.mkdir(tgz_dir)
    except FileExistsError:
      pass

    files_to_tar = {}
    for tool in tool_list:
      if not tool in files_to_tar: files_to_tar[tool] = []
      for f in sorted(qc_files):
          files_to_tar[tool].append(f)

    for tool in tool_list:
      if not files_to_tar[tool]: continue
      tarfile_name = f"{tgz_dir}/{sampleId}.{tool}.tgz"
      with tarfile.open(tarfile_name, "w:gz", dereference=True) as tar:
        for f in files_to_tar[tool]:
          tar.add(f, arcname=os.path.basename(f))

def main(args):
    with open(args.seq_experiment_analysis, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_experiment_analysis_dict = row

    pipeline_info = {}
    updated_pipeline_info = {}
    if args.pipeline_yml:
      with open(args.pipeline_yml, 'r') as f:
        pipeline_info = yaml.safe_load(f)
        for tool, version in pipeline_info.items():
            # tool are like NFCORE_VARIANTCALL:VARIANTCALL:PICARD_LIFTOVERVCF_RA, NFCORE_VARIANTCALL:VARIANTCALL:XML_VCF
            new_tool_parts = tool.split(":")[-1].split("_")
            # Join the first two parts if there are multiple parts
            new_tool = "_".join(new_tool_parts[:2]) if len(new_tool_parts) > 1 else new_tool_parts[0]
            updated_pipeline_info[new_tool]=version
    if seq_experiment_analysis_dict.get('analysis_tools (tools and versions)'):
        tools_dict = dict(tool.split(' ', 1) for tool in seq_experiment_analysis_dict.get('analysis_tools (tools and versions)').split(', '))
        updated_pipeline_info['FoundationOneCDx'] = tools_dict
        for key, value in updated_pipeline_info.items():
            for sub_key, sub_value in value.items():
                value[sub_key] = str(sub_value)
            updated_pipeline_info[key] = value
    else:
        updated_pipeline_info['FoundationOneCDx'] = ""


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
        'analysisType': { 'name': 'variant_calling_supplement' },
        'variant_calling_strategy': seq_experiment_analysis_dict.get('variant_calling_strategy'),
        'workflow': {
            'genome_build': 'GRCh38',
            'workflow_name': seq_experiment_analysis_dict.get('workflow_name'),
            'workflow_version': seq_experiment_analysis_dict.get('workflow_version'),
            'workflow_short_name': seq_experiment_analysis_dict.get('workflow_short_name'),
            'pipeline_info': updated_pipeline_info,
            'run_id':args.wf_run
            
        },
        'experiment' : {
            'submitter_sequencing_experiment_id': seq_experiment_analysis_dict.get('submitter_sequencing_experiment_id'),
            'platform' : seq_experiment_analysis_dict.get('platform'),
            'experimental_strategy': seq_experiment_analysis_dict.get('experimental_strategy'),
            'target_capture_kit': seq_experiment_analysis_dict.get('target_capture_kit'),
            'primary_target_regions': seq_experiment_analysis_dict.get('primary_target_regions'),
            'capture_target_regions': seq_experiment_analysis_dict.get('capture_target_regions'),
            'coverage': seq_experiment_analysis_dict.get('coverage').split(',')
        },
        'variant_class' : seq_experiment_analysis_dict.get('variant_class'),  # update to "Somatic/Germline" seq_experiment_analysis_dict.get('coverage') after schema update
        'files': []
    }

    tmp=payload['variant_class'].split(",")
    tmp.sort()

    if len(tmp)>1:
        payload['variant_class']="+".join(tmp)
    else:
        payload['variant_class']=tmp[0]




    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    # generate date string
    date_str = date.today().strftime("%Y%m%d")

    # prepare tarball to include all QC files generated by one tool
    #sampleId=retrieveClinical(payload['samples'][0]['submitterSampleId'])
    print(args.sample_id)
    prepare_tarball(args.sample_id, args.files_to_upload, ["FoundationMedicine"])

    # get file of the payload
    date_str = date.today().strftime("%Y%m%d")

    for f in sorted(glob('tarball/*.tgz')):
      print("WHAT",args.donor_id,args.sample_id)
      file_info = get_files_info(f, date_str,seq_experiment_analysis_dict,args.donor_id,args.sample_id)
      payload['files'].append(file_info)


    # for f in args.files_to_upload:
    #     renamed_file = rename_file(f, payload, seq_experiment_analysis_dict, date_str)
    #     payload['files'].append(get_files_info(renamed_file))

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
    parser.add_argument("-c", "--sample_id", dest="sample_id", required=True, help="sample_id")
    parser.add_argument("-d", "--donor_id", dest="donor_id", required=True, help="donor_id")

    args = parser.parse_args()

    main(args)
