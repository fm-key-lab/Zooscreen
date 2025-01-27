# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import os

# %% 
os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_screening/zooscreen_final')

# %% 
# get sample order data
#get name conversions
sample_name_conversions=np.loadtxt('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_screening/sample_names_conversion_final.tsv',dtype=str)
dict_convert_sample_names={x[2]:x[3] for x in sample_name_conversions}

zooscreen_heatmap_data = pd.read_csv('zooscreen_inhouse/heatmap_overview_Wevid.tsv', sep='\t', header=0,index_col=0)
tartu_heatmap_data = pd.read_csv('200423_tartu_data/heatmap_overview_Wevid.tsv', sep='\t', header=0,index_col=0)
first_tartu_data = None
# %%
# combine data
for_export=pd.concat([tartu_heatmap_data, zooscreen_heatmap_data], axis=1, join='outer').fillna(1)

# rename sample names
columns=[]
for c in for_export.columns:
    c_stripped=c.replace('.unmapped.rma6','')
    if c_stripped in dict_convert_sample_names:
        columns.append(dict_convert_sample_names[c_stripped])
    elif 'UDI' in c_stripped:
        pass
    else:
        columns.append(f'{c_stripped.split(".")[3]}-{c_stripped.split(".")[4].split("_")[0]}')


columns=[dict_convert_sample_names[x.replace('.unmapped.rma6','')] for x in for_export.columns]
for_export.columns=columns
# set sample order
#for_export=for_export[[x for x in sample_order_data if x in columns]]
# reorder species
#for_export=for_export.reindex(species_order_data)
# standardize naming to human readable
for_export=for_export.rename(index=lambda x: x.replace('_', ' '))

# %%
# convert to boolean and export
for_export = for_export==4
for_export=for_export.sort_index()

for_export.to_csv('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_screening/zooscreen_final/zooscreen_heatmap_merged.tsv',sep='\t')

