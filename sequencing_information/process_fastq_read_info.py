# %%
# load modules
import numpy as np
import pandas as pd
import os

# %% paths
os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/sequencing_informatino')

path_raw_fastq_counts='filtred_reads.tsv'

path_fastqs_to_sample='seq_id_to_azp_id.csv'
# %% 
# import dataframes
raw_fastq_counts=pd.read_csv(path_raw_fastq_counts,sep='\t',header=0)
fastqs_to_sample=pd.read_csv(path_fastqs_to_sample,sep=',')
# %% 
# collapse files into samples
uncollapsed_reads_passed={}
total_reads={}
for lab_id,seq_id in zip(fastqs_to_sample['Lab ID'],fastqs_to_sample['Seq ID']):
    entries_with_labid=raw_fastq_counts.loc[raw_fastq_counts['Sample'].str.contains(seq_id)]
    uncollapsed_reads_passed[lab_id]=np.sum(entries_with_labid['Passed Filter'])
    total_reads[lab_id]=np.sum(entries_with_labid['Passed Filter'])+np.sum(entries_with_labid['Low Quality'])+np.sum(entries_with_labid['Too Many N'])+np.sum(entries_with_labid['Too Short'])+np.sum(entries_with_labid['Too Long'])


reads_df_for_output=pd.DataFrame({'Total Reads':total_reads,'Passed Reads':uncollapsed_reads_passed},dtype=int)
reads_df_for_output.to_csv('raw_to_qc_reads.tsv',sep='\t')
# %%
