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
    Edmund Su
"""

from email import contentmanager
import os
import sys
import argparse
import subprocess
import json
import csv
import requests


def main():
    """
    Python implementation of tool: sanity-check

    This is auto-generated Python code, please update as needed!
    """

    parser = argparse.ArgumentParser(description='Tool: sanity-check')
    parser.add_argument('-t', '--api_token', dest='api_token', type=str,
                        help='Platform UUID token', required=True)
    parser.add_argument('-c', '--clinical_url', dest='clinical_url', type=str,
                        help='Clinical API interaction URL', required=True)
    parser.add_argument('-s', '--submission_song_url', dest='submission_song_url',type=str,
                        help='Submission SONG API interaction URL', required=True)
    parser.add_argument('-x', '--experiment_info_tsv', dest='experiment_info_tsv', type=str,required=True,
                        help='TSV containing experiment info and submitter IDs')
    parser.add_argument('-f', '--force', action='store_true',help='Skip step checking for sample + experiment duplicates')
    args = parser.parse_args()

    if args.experiment_info_tsv:
        metadata=load_tsv(args.experiment_info_tsv)

        clinical_metadata=get_clinical(
            metadata,
            args.clinical_url,
            args.api_token
        )

        final_metadata=compile_metadata(
            metadata,
            clinical_metadata
        )
        check_study_exists(
                final_metadata,
                args.submission_song_url
        )

        if not args.force:
            check_analysis_exists(
                final_metadata,
                args.submission_song_url
            )

        update_tsv(final_metadata,"updated_"+args.experiment_info_tsv)


def load_tsv(experiment_info_tsv):
    metadata_dict = {}
    with open(experiment_info_tsv, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
        if len(rows) != 1:
            sys.exit("Error found: experiment TSV expects exactly one data row, offending file: %s has %s row(s)\n" % \
                (experiment_info_tsv, len(rows)))
        metadata_dict.update(rows[0])

    # Early exit if any submitter IDs missing
    for field in ['donor','sample','specimen']:
        if metadata_dict.get("submitter_"+field+"_id")==None:
            sys.exit("submitter_%s_id was not found. Please resubmit %s with appropriate field" % (field,experiment_info_tsv))

    return metadata_dict

def get_clinical(metadata,clinical_url,api_token):

    return_metadata={}

    ### Check program exists first

    headers={"accept":"*/*","Authorization":"Bearer %s" % (api_token)}
    endpoint="%s/clinical/program/%s/donors" % (clinical_url,metadata.get('program_id'))
    response=requests.get(endpoint,headers=headers)

    if response.status_code==200:
        if len(response.text)==0:
            sys.exit("Project %s does not exist or no samples have been registered" % (metadata.get('program_id')))
    else:
        sys.exit("Unable to fetch. Status code : %s" % (str(response.status_code)))


    ###Populate ARGO IDs
    headers={"accept":"text/plain","Authorization":"Bearer %s" % (api_token)}

    for field in ['donor','sample','specimen']:
        endpoint="%s/clinical/%ss/id?programId=%s&submitterId=%s" % (clinical_url,field,metadata.get('program_id'),metadata.get("submitter_"+field+"_id"))

        response=requests.get(endpoint,headers=headers)

        if response.status_code==404:
            sys.exit("submitter_%s_id:'%s' was not found in project:'%s'. Verify sample has been registered." % (field,metadata.get("submitter_"+field+"_id"),metadata.get('program_id')))
        elif response.status_code!=200:
            sys.exit("Unable to fetch. Status code : %s" % (str(response.status_code)))
        else:
            return_metadata[field+"_id"]=response.content.decode("utf-8")

    ###Populate w/ Clinical
    endpoint="%s/clinical/program/%s/donor/%s" % (clinical_url,metadata.get('program_id'),return_metadata.get('donor_id'))
    response=requests.get(endpoint,headers=headers)
    if response.status_code!=200:
        sys.exit("Unable to fetch. Status code : %s" % (str(response.status_code)))
    else:
        return_metadata['gender']=response.json()['gender']
        return_metadata['submitter_donor_id']=response.json()['submitterId']
        return_metadata['program_id']=response.json()['programId']

        specimen_ind=[ele for ele,specimen in enumerate(response.json()['specimens']) if specimen['specimenId']==return_metadata['specimen_id']]
        if len(specimen_ind)!=1:
            sys.exit("ID Mismatch detected. Specimen_id:'%s'/'%s' was not found within Donor:'%s'/'%s' 's specimens" % (metadata['submitter_specimen_id'],return_metadata['specimen_id'],metadata['submitter_donor_id'],return_metadata['donor_id']))
        return_metadata['specimen_tissue_source']=response.json()['specimens'][specimen_ind[0]]['specimenTissueSource']
        return_metadata['tumour_normal_designation']=response.json()['specimens'][specimen_ind[0]]['tumourNormalDesignation']
        return_metadata['specimen_type']=response.json()['specimens'][specimen_ind[0]]['specimenType']
        return_metadata['submitter_specimen_id']=response.json()['specimens'][specimen_ind[0]]['submitterId']

        sample_ind=[ele for ele,sample in enumerate(response.json()['specimens'][specimen_ind[0]]['samples']) if sample['sampleId']==return_metadata['sample_id']]
        if len(sample_ind)!=1:
            sys.exit("ID Mismatch detected. Sample_id:'%s'/'%s'  was not found within Specimen:'%s'/'%s' 's samples" % (metadata['submitter_sample_id'],return_metadata['sample_id'],metadata['submitter_specimen_id'],return_metadata['specimen_id']))
        return_metadata['sample_type']=response.json()['specimens'][specimen_ind[0]]['samples'][sample_ind[0]]['sampleType']
        return_metadata['submitter_sample_id']=response.json()['specimens'][specimen_ind[0]]['samples'][sample_ind[0]]['submitterId']

    if return_metadata['tumour_normal_designation']=="Tumour":
        #WGS, WXS, RNA-Seq, Bisulfite-Seq, ChIP-Seq, Targeted-Seq
        if metadata.get("submitter_matched_normal_sample_id"):
            check_normal_sample_exists(metadata,response.json())
        else:
            if metadata.get("experimental_strategy") in ['WGS','WXS']:
                sys.exit("Null entry for `submitter_matched_normal_sample_id` detected. For tumour `experiment_strategy` type %s ,this field is required and must reference a registered normal sample." % (metadata.get("experimental_strategy")))

    return return_metadata



def check_normal_sample_exists(metadata,clinical_metadata):
    submitter_id=metadata['submitter_matched_normal_sample_id']

    return_id=None
    tumourNormalDesignation=None

    for specimen in clinical_metadata['specimens']:
        for samples in specimen['samples']:
            if samples['submitterId']==submitter_id:
                return_id=samples['submitterId']
                tumourNormalDesignation=specimen["tumourNormalDesignation"]
    if return_id==None:
        sys.exit("'submitter_matched_normal_sample_id':%s was not found in study. Please verify '%s' has been registered." % (submitter_id,submitter_id))
    if tumourNormalDesignation=="Tumour":
        sys.exit("'submitter_matched_normal_sample_id':%s detected as tumour instead of normal. Please verify correct sample." % (submitter_id))


def compile_metadata(metadata,clinical_metadata):

    ###Over-write metadata with clinical_metadata
    return_metadata=metadata.copy()
    for key in clinical_metadata.keys():
        return_metadata[key]=clinical_metadata[key]
    return return_metadata

def check_study_exists(metadata,submission_song_url):
    headers={"accept":"*/*"}
    endpoint="%s/studies/%s" % (submission_song_url,metadata['program_id'])
    response=requests.get(endpoint,headers=headers)

    if response.status_code==404:
        sys.exit("Program %s does not exist in SONG. Please verify program code is correct. Otherwise contact DCC-admin for help to troubleshoot." % (metadata.get('program_id')))
    elif response.status_code!=200:
        sys.exit("Unable to fetch. Status code : %s" % (str(response.status_code)))
    else:
        return True

def check_analysis_exists(metadata,submission_song_url):
    headers={"accept":"*/*"}
    endpoint="%s/studies/%s/analysis/search/id?sampleId=%s" % (submission_song_url,metadata.get('program_id'),metadata.get('sample_id'))
    error_message=[]
    analysis_exists=False
    response=requests.get(endpoint,headers=headers)

    if response.status_code!=200:
        sys.exit("Unable to fetch. Status code : %s" % (str(response.status_code)))
    elif response.status_code==200:
        ### If first submission of sample, query will return empty list
        if len(response.json())==0:
            return True
        else:
        ### Handles instances where sample has existing analyses
            for analysis in response.json():
                ### If analysis is suppressed; we ignore
                if \
                analysis["analysisState"]=="PUBLISHED" and \
                analysis["experiment"]["experimental_strategy"]==metadata.get('experimental_strategy') and \
                analysis['analysisType']['name']=="sequencing_experiment"\
                :
                    analysis_exists=True
                    error_message.append(
                        "Sample '%s'/'%s' has an existing published analysis '%s' for experiment_strategy '%s.'" \
                        % \
                        (metadata.get('submitter_sample_id'),metadata.get('sample_id'),analysis['analysisId'],metadata.get('experimental_strategy'))
                    )
            if analysis_exists:
                if len(error_message)>1:
                    sys.exit("Too many conflict analyses detected. Displaying subset:\n"+"\n".join(error_message[:5]))
                else:
                    sys.exit(error_message[0])
            else:
                return True
    else:
        return False

def update_tsv(metadata,file):

    # Remove ARGO IDs, payload-gen doesn't like them
    # metadata.pop("donor_id")
    # metadata.pop("specimen_id")
    # metadata.pop("sample_id")
    with open(file, 'w', newline='') as csvfile:
        fieldnames = metadata.keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames,delimiter='\t')
        writer.writeheader()
        writer.writerow(metadata)

if __name__ == "__main__":
    main()
