# %%
import pandas as pd
import numpy as np
import os

os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_mapping/erysipelothrix_rhusiopathiae')
# %%
biosamples_text_path='erysipelothrix_biosamples.txt'

output_dict={}

with open(biosamples_text_path) as biosamples_summary:
    for line in biosamples_summary.readlines():
        line=line.strip()
        if line.strip().split(':')[0].isdigit():
            identifier=line.strip().split(':')[0]
            output_dict[identifier]={}
        else: 
            if line.startswith('Identifiers:'):
                ids=line.strip('Identifiers:').split('; ')
                for id_entry in ids:
                    id_identity=id_entry.strip().split(': ')[0]
                    id_value=id_entry.strip().split(': ')[1]
                    output_dict[identifier][id_identity]=id_value
            else:
                dict_entry=line.strip().strip('/').split('=')
                if len(dict_entry)==2:
                    output_dict[identifier][dict_entry[0]]=dict_entry[1]
# %%
def exclude_due_to_check_flags(data_entry):
    flags_that_exclude_raw_data=['MAG','MAG coverage depth']
    for f in flags_that_exclude_raw_data:
        if f in data_entry:
            return True
    return False

def exclude_due_to_data_in_fields(data_entry):
    if 'sample type' in data_entry:
        if data_entry['sample type'] == '"metagenomic assembly"':
            return True
    if 'investigation type' in data_entry: 
        investigation_type=data_entry['investigation type']
        if investigation_type == '"metagenome"' or investigation_type == '"metagenome-assembled genome"':
            return True
    return False

def exclude_due_to_non_sra(data_entry):
    if 'SRA' not in data_entry:
        return True
    elif data_entry['SRA'].startswith('DRS'):
        return True
    else:
        return False

def check_host_or_source(data):
    na_values_host=['"unknown"','"not available"','"missing"','"Not Applicable"','"Vaccine"','']
    na_values_isolation=['"not available: to be reported later"', '"missing"', '']
    host_data=''
    isolate_data=''
    if 'host' in data:
        host_data=data['host']
    if 'isolation source' in data:
        isolate_data=data['isolation source']
    # compare
    if host_data in na_values_host and isolate_data in na_values_isolation:
        to_return='NA'
    elif host_data in na_values_host:
        to_return= isolate_data
    elif isolate_data in na_values_host:
        to_return= host_data
    else:
        to_return= host_data
    to_return=to_return.strip('"')
    return to_return

def check_geographic(data):
    geographic_info=''
    if 'geographic location' in data:
        geographic_info=data['geographic location']
    else: 
        geographic_info='NA'
    geographic_info=geographic_info.split(':')[0].strip('"')
    if geographic_info in ['Not Applicable','missing','not determined']:
        geographic_info='NA'
    return geographic_info

def extract_sample_name(data):
    if 'Sample name' in data:
        to_export=data['Sample name']
    elif 'sample name' in data:
        to_export=data['sample name']
    elif 'strain' in data: 
        to_export=data['strain']
    return to_export.strip('"').replace(' ','_')

sra_info=[]
host_info=[]
geographic_info=[]
sample_names=[]
biosample_info=[]
for entry in output_dict:
    data=output_dict[entry]
    # removal of non-raw data
    if not exclude_due_to_check_flags(data) and not exclude_due_to_data_in_fields(data) and not exclude_due_to_non_sra(data):
        host_info.append(check_host_or_source(data))
        geographic_info.append(check_geographic(data))
        sra_info.append(data['SRA'])
        sample_names.append(extract_sample_name(data))
        biosample_info.append(data['BioSample'])
# output for downloading, manipulation

parsed_metadata=pd.DataFrame(data={'sample_name':sample_names,
                   'SRA':sra_info,
                   'host':host_info,
                   'geographic_loc':geographic_info,
                   'biosample_info':biosample_info})

parsed_metadata.to_csv('comparison_samples_to_download.csv',index=False)

# %%
biosamples_to_check=pd.read_csv('rhusio_filtered_genbank_samplenames_accessions.tsv',header=None, names=['GCA','biosample_acc'],sep='\t')

# %%
output_list=[]
for index,row in biosamples_to_check.iterrows():
    gca_accession=row['GCA']
    biosamples_accession=row['biosample_acc']
    metadata_row=parsed_metadata.query(f'biosample_info == "{biosamples_accession}"').reset_index()
    if len(metadata_row)==1:
        host=metadata_row["host"].loc[0].replace(' ','_')
        geo=metadata_row["geographic_loc"].loc[0].replace(' ','_')
        output_list.append(f'{gca_accession}_Host_{host}_Geo_{geo}')
    else:
        output_list.append(f'{gca_accession}_Host_unknown_Geo_unknown')

pd.DataFrame({'GCA':biosamples_to_check['GCA'],'update_name':output_list}).to_csv('name_conversion_genbank.csv',sep=',')
# %%

