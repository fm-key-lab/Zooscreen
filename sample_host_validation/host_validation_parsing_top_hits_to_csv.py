# %% 
import numpy as np
import pandas as pd
import os
import gzip
import glob
import csv
# %% 
os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/sample_host_validation')

# %%
# supporting metadata matching
ids_to_zoological_id_path='/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/sample_host_validation/sample_to_zooid.tsv'
ids_to_zoological_id=pd.read_csv(ids_to_zoological_id_path,sep='\t',header=None).fillna('')
ids_to_zoological_id_conversion_dict={x:y for x,y in zip(ids_to_zoological_id[0],ids_to_zoological_id[1])}

unique_zooids=np.unique([x for x in ids_to_zoological_id_conversion_dict.values()])


tartu_ids_to_azp_path='/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/sample_host_validation/tartu_to_azp.tsv'
tartu_ids_to_azp=pd.read_csv(tartu_ids_to_azp_path,sep='\t',header=None)
tartu_ids_to_azp_conversion_dict={x:y for x,y in zip(tartu_ids_to_azp[0],tartu_ids_to_azp[1])}
all_to_fill_double_extract_sample={}
for x,y in zip(tartu_ids_to_azp[0],tartu_ids_to_azp[1]):
    if x in all_to_fill_double_extract_sample:
        all_to_fill_double_extract_sample[x].append(y.replace('_','-'))
    else:
        all_to_fill_double_extract_sample[x]=[y.replace('_','-')]
# %%
def parse_gzipped_tsv_kingdom_metazoa(file_path):
    data = []
    collecting = False
    with gzip.open(file_path, 'rt') as f:
        reader = csv.reader(f,  delimiter='\t')
        for row in reader:
            if row[5].strip() == 'K':
                if collecting:
                    break
                if row[7].strip() == 'Metazoa':
                    collecting = True
            if collecting:
                data.append([x.strip() for x in row])

    return pd.DataFrame(data, columns=['proportion_reads', 'total_reads_at_or_below_node', 'total_reads_assigned_node', 'minimizers', 'distinct_minimizers', 'rank', 'taxid', 'name'])

def parse_gzipped_tsv_total_reads_assinged(file_path):
    results_dict={}
    with gzip.open(file_path, 'rt') as f:
        reader = csv.reader(f,  delimiter='\t')
        for row in reader:
            if row[7].strip() == 'unclassified':   
                results_dict['unclassified']=row[1]
            if row[7].strip() == 'root':
                results_dict['classified']=row[1]            
                break
    return results_dict

# %% 
contaminated_eukaryotic_genuses=['Cyprinus']
name=[]
zoological_assignment=[]
top_1=[]
top_2=[]
top_3=[]
top_1_per=[]
top_2_per=[]
top_3_per=[]
classified=[]
unclassified=[]

for index,file_path in enumerate(glob.glob('k2_nt/*.gz')):
    this_file_basename=file_path.split('/')[1].split('_pe_k2_nt')[0]
    if 'AZP' in this_file_basename and '-E' not in this_file_basename:
        azp_id='AZP-'+this_file_basename.split('AZP-')[1].split('_')[0]
        zoological_assign=ids_to_zoological_id_conversion_dict[azp_id]
        this_file_parsed=parse_gzipped_tsv_kingdom_metazoa(file_path)
        classified_unclassified=parse_gzipped_tsv_total_reads_assinged(file_path)
        classified.append(int(classified_unclassified['classified']))
        unclassified.append(int(classified_unclassified['unclassified']))
        this_file_parsed["proportion_reads"] = pd.to_numeric(this_file_parsed["proportion_reads"])
        subset_df=this_file_parsed.query('rank == "G"').sort_values(by='proportion_reads',ascending=False)
        filtered_df = subset_df[~subset_df['name'].isin(contaminated_eukaryotic_genuses)]
        name.append(azp_id)
        zoological_assignment.append(zoological_assign)
        top3=np.array(filtered_df[['name','proportion_reads']].sort_values('proportion_reads',ascending=False)[:3]).flatten().tolist()
        top_1.append(top3[0])
        top_1_per.append(top3[1])
        if len(top3)>2:
            top_2.append(top3[2])
            top_2_per.append(top3[3])
        else:
            top_2.append('NA')
            top_2_per.append(np.NaN)
        if len(top3)>4:
            top_3.append(top3[4])
            top_3_per.append(top3[5])
        else:
            top_3.append('NA')
            top_3_per.append(np.NaN)
    elif 'Pool' not in this_file_basename:
        azp_id_list=all_to_fill_double_extract_sample[this_file_basename]
        for azp_id in azp_id_list:
            zoological_assign=ids_to_zoological_id_conversion_dict[this_file_basename]
            this_file_parsed=parse_gzipped_tsv_kingdom_metazoa(file_path)
            this_file_parsed["proportion_reads"] = pd.to_numeric(this_file_parsed["proportion_reads"])
            classified_unclassified=parse_gzipped_tsv_total_reads_assinged(file_path)
            classified.append(int(classified_unclassified['classified']))
            unclassified.append(int(classified_unclassified['unclassified']))
            subset_df=this_file_parsed.query('rank == "G"').sort_values(by='proportion_reads',ascending=False)
            filtered_df = subset_df[~subset_df['name'].isin(contaminated_eukaryotic_genuses)]
            name.append(azp_id)
            zoological_assignment.append(zoological_assign)
            top3=np.array(filtered_df[['name','proportion_reads']].sort_values('proportion_reads',ascending=False)[:3]).flatten().tolist()
            top_1.append(top3[0])
            top_1_per.append(top3[1])
            if len(top3)>2:
                top_2.append(top3[2])
                top_2_per.append(top3[3])
            else:
                top_2.append('NA')
                top_2_per.append(np.NaN)
            if len(top3)>4:
                top_3.append(top3[4])
                top_3_per.append(top3[5])
            else:
                top_3.append('NA')
                top_3_per.append(np.NaN)



dict_for_df={'Name':name,
             'ZooID':zoological_assignment,
             'kraken2_top1_euk_genus':top_1,
             'kraken2_top1_euk_genus_perc':top_1_per,
             'kraken2_top2_euk_genus':top_2,
             'kraken2_top2_euk_genus_perc':top_2_per,
             'kraken2_top3_euk_genus':top_3,
             'kraken2_top3_euk_genus_perc':top_3_per,
             'classified':classified,
             'unclassified':unclassified}

test=pd.DataFrame(data=dict_for_df)


# %%

zooid_to_genus={x:y for x,y in [('Aurochs', 'Bos'),
('Beaver', 'Meles'),
('Cattle', 'Bos'),
('Cattle subadult (ca. 6 years)','Bos'),
('Dog', 'Canis'),
('Dog?', 'Canis'),
('Equid', 'Equus'),
('Fox', 'Vulpes'),
('Fox adult','Vulpes'),
('Gazella', 'Gazella'),
('Goat', 'Capra'),
('Spermophilus musicus','Spermophilus'),
('Indeterminate', 'Unknown'),
('Mammals', 'Mammals'),
('Human', 'Homo'),
('Pig', 'Sus'),
('Pig-size', 'Sus'),
('Pig/wild pig', 'Sus'),
('Pig (domestic or wild)','Sus'),
('Red deer', 'Cervus'),
('Sheep', 'Ovis'),
('Sheep/goat', 'Ovis'),
('Sheep-size', 'Ovis'),
('Sheep goat or red deer','Ovis'),
('Small ruminant', 'Ruminant'),
('Wild pig', 'Sus'),
('cow-size', 'Bos'),
('small ruminant', 'Ruminant'),
('Horse','Equus'),
('Bear','Ursus'),
('','Unknown'),
('Serpentes','Vipera')]}

multiple_zooid_to_genus={
    'Sheep/goat': ['Ovis','Capra'],
    'Cattle or elk': ['Bos','Cervus'],
    'Cattle/ sheep/ goat':['Bos','Ovis','Capra'],
    'Sheep/ goat':['Ovis','Capra'],
    'small ruminant': ['Ovis','Capra','Bos'],
    'Small ruminant': ['Ovis','Capra','Bos'],
    'Mammals':['Ovis','Capra','Bos','Homo','Sus','Cervus','Meles','Canis','Vulpes','Equus'],
    '':['Ovis','Capra','Bos','Sus','Cervus','Meles','Canis','Vulpes','Equus','Chionomys'],
    'Indeterminate':['Ovis','Capra','Bos','Sus','Cervus','Meles','Canis','Vulpes','Equus','Chionomys']
}
multiple_zooid_to_output={
    'Cattle or elk': 'Bos/Cervus',
    'Sheep/ goat': 'Ovis/Capra',
    'Cattle/ sheep/ goat':'Bos/Ovis/Capra',
    'Sheep/goat': 'Ovis/Capra',
    'small ruminant': 'Small ruminant',
    'Small ruminant': 'Small ruminant',
    'Mammals':'Mammal',
    '':'Unknown',
    'Indeterminate':'Unknown'
}

# %%
# check for confirmation if >0.1% reads assign to same genus
# alert for swap if >1% reads assign to different genus
swap_status=[]
swap_cutoff=1
confirmation_cutoff=0.1
final_id=[]
for index,row in test.iterrows():
    zooid=row['ZooID']
    top_genus=row['kraken2_top1_euk_genus']
    top_genus_perc=row['kraken2_top1_euk_genus_perc']
    if zooid in multiple_zooid_to_genus:
        # check for possibilities of match to consistent genus
        if top_genus in multiple_zooid_to_genus[zooid]:
            final_id.append(top_genus)
            swap_status.append('multiple - agree')
        elif top_genus_perc >= swap_cutoff:
            final_id.append(top_genus)
            swap_status.append('multiple - swap')
        else: 
            final_id.append(multiple_zooid_to_output[zooid])
            swap_status.append('multiple -no swap')
    else:
        if zooid_to_genus[zooid]==top_genus and top_genus_perc >= confirmation_cutoff:
            print('confirmed')
            swap_status.append('agree')
            final_id.append(top_genus)
        elif zooid_to_genus[zooid]==top_genus and top_genus_perc < confirmation_cutoff:
            print('confirmed, but low data')
            swap_status.append('agree - low data')
            final_id.append(top_genus)
        elif zooid_to_genus[zooid]!=top_genus and top_genus_perc >= swap_cutoff:
            print('swapped',f'ZooID Genus: {zooid}, Kraken genus: {top_genus}, {top_genus_perc}% reads assigned')
            swap_status.append('disagree - swap supported')
            final_id.append(top_genus)
        elif zooid_to_genus[zooid]!=top_genus and top_genus_perc < confirmation_cutoff:
            print('inconsistent, below both cutoffs',f'ZooID Genus: {zooid}, Kraken genus: {top_genus}, {top_genus_perc}% reads assigned')
            swap_status.append('disagree - low data')
            final_id.append(zooid_to_genus[zooid])
        elif zooid_to_genus[zooid]!=top_genus and top_genus_perc < swap_cutoff:
            print('inconsistent, below swap cutoff',f'ZooID Genus: {zooid}, Kraken genus: {top_genus}, {top_genus_perc}% reads assigned')
            swap_status.append('disagree - swap not supported')
            final_id.append(zooid_to_genus[zooid])

test['comparison_result']=swap_status
test['final_id']=final_id
comparison_dataset=test[['Name', 'ZooID', 'kraken2_top1_euk_genus',
       'kraken2_top1_euk_genus_perc','final_id','comparison_result','classified','unclassified']]
# %%
to_clean_for_output=pd.DataFrame(data=test)

clean_name=[]
for x in to_clean_for_output['Name']:
    clean_azp=x.split('AZP-')[1].split('_')[0]
    clean_name.append(f'AZP-{clean_azp}')
to_clean_for_output['Name']=clean_name

import re

def extract_numeric(name):
    match = re.search(r'\d+', name)
    return int(match.group()) if match else float('inf')

# Create a new column with the numeric part
to_clean_for_output['NumericPart'] = to_clean_for_output['Name'].apply(extract_numeric)

# Sort the DataFrame by the new column
df_sorted = to_clean_for_output.sort_values(by='NumericPart')
df_sorted = df_sorted.drop(columns=['NumericPart'])

df_for_output=df_sorted[['Name','final_id', 'ZooID', 'kraken2_top1_euk_genus',
       'kraken2_top1_euk_genus_perc', 'kraken2_top2_euk_genus',
       'kraken2_top2_euk_genus_perc', 'kraken2_top3_euk_genus',
       'kraken2_top3_euk_genus_perc', 'comparison_result','classified','unclassified'
       ]]
# %%
df_for_output.to_csv('kraken2_host_top_3_eukaryotic_genus_with_comparisson.csv',index=False)
