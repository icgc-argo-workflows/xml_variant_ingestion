# import pandas as pd
# import yaml
import csv

# with open('updated_sample.tsv', 'r') as f:
#         seq_experiment_analysis_df = pd.read_csv(f, sep='\t')
#         seq_experiment_analysis_dict = seq_experiment_analysis_df.iloc[0].to_dict()

# print(seq_experiment_analysis_dict)

# with open('collated_versions.yml', 'r') as f:
#         pipeline_info = yaml.safe_load(f)

# for key, value in pipeline_info.items():
#     for sub_key, sub_value in value.items():
#         value[sub_key] = str(sub_value)
#     pipeline_info[key] = value

# print(pipeline_info)

# tools_dict = dict(tool.split(' ', 1) for tool in seq_experiment_analysis_dict.get('analysis_tools (tools and versions)').split(', '))
# pipeline_info['FoundationOneCDx'] = tools_dict

# print(pipeline_info)

with open('updated_sample.tsv', 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_experiment_analysis_dict = row

print(seq_experiment_analysis_dict)
