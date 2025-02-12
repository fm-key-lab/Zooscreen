# %% 
import numpy as np
import pandas as pd
import os
import gzip
import glob
import csv
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats as sts
from statsmodels.stats.multitest import multipletests
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

zooscreen_pathogen_samples=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_screening/zooscreen_final/zooscreen_heatmap_merged.tsv', sep='\t',index_col=0)

# metadata on sample origins:
zooscreen_sample_sites_metadata=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/sequencing_informatino/AZP supplementary tables - Table S1_ Site info.tsv', sep='\t',skiprows=[0,36,37,38,39,40,41,42,43,44,45,46,47],header=0)
zooscreen_samples_metadata=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/sequencing_informatino/AZP supplementary tables - Table S2_ Sample info.tsv', sep='\t',skiprows=[0],header=0)

# fill counties for site data
countries=[]
for x in zooscreen_sample_sites_metadata['Country']:
    if pd.isna(x):
        countries.append(previous)
    else:
        previous=x
        countries.append(previous)
zooscreen_sample_sites_metadata['Country']=countries

# dictionary mapping site to country:
countries_not_considered_for_pathology=np.unique(zooscreen_sample_sites_metadata[(zooscreen_sample_sites_metadata['Investigated for Paleopathology']!='Yes')]['Country'].to_list())
site_to_country_mapping={x:y for x,y in zip(zooscreen_sample_sites_metadata['Site name'].to_list(),zooscreen_sample_sites_metadata['Country'].to_list())}
site_to_country_mapping['Augsburg']='Germany (Augsburg)'
labid_to_site_mapping={x:y for x,y in zip(zooscreen_samples_metadata['Lab ID'].to_list(),zooscreen_samples_metadata['Archaeological site'].to_list())}

# fill in if assessed site-by-site pathology
pathology_assessed=[]
for x in zooscreen_samples_metadata['Archaeological site']:
    if site_to_country_mapping[x] not in countries_not_considered_for_pathology:
        pathology_assessed.append(True)
    else:
        pathology_assessed.append(False)

zooscreen_samples_metadata['pathology_assessed']=pathology_assessed
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

def parse_gzipped_tsv_bacteria(file_path):
    data = []
    collecting = False
    with gzip.open(file_path, 'rt') as f:
        reader = csv.reader(f,  delimiter='\t')
        for row in reader:
            if row[5].strip() == 'D':
                if collecting:
                    break
                if row[7].strip() == 'Bacteria':
                    collecting = True
            if collecting:
                data.append([x.strip() for x in row])
    return pd.DataFrame(data, columns=['proportion_reads', 'total_reads_at_or_below_node', 'total_reads_assigned_node', 'minimizers', 'distinct_minimizers', 'rank', 'taxid', 'name'])

def mean_confidence_interval(data, confidence=0.95): 
    # adopted from https://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data#15034143
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), sts.sem(a)
    h = se * sts.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h
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

parsed_bacteria={}

def helper_extract_genus(file_path,azp_id):
    this_file_parsed=parse_gzipped_tsv_kingdom_metazoa(file_path)
    this_file_parsed["proportion_reads"] = pd.to_numeric(this_file_parsed["proportion_reads"])
    classified_unclassified=parse_gzipped_tsv_total_reads_assinged(file_path)
    classified.append(int(classified_unclassified['classified']))
    unclassified.append(int(classified_unclassified['unclassified']))
    subset_df=this_file_parsed.query('rank == "G"').sort_values(by='proportion_reads',ascending=False)
    filtered_df = subset_df[~subset_df['name'].isin(contaminated_eukaryotic_genuses)]
    name.append(azp_id)
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
    
def helper_extract_genus_bacteria(file_path):
    this_file_parsed_bacteria=parse_gzipped_tsv_bacteria(file_path)
    this_file_parsed_bacteria["proportion_reads"] = pd.to_numeric(this_file_parsed_bacteria["proportion_reads"])

for index,file_path in enumerate(glob.glob('k2_nt/*.gz')):
    this_file_basename=file_path.split('/')[1].split('_pe_k2_nt')[0]
    if 'AZP' in this_file_basename and '-E' not in this_file_basename:
        azp_id='AZP-'+this_file_basename.split('AZP-')[1].split('_')[0]
        if azp_id != 'AZP-289': # false mapping of library indexes, skip
            helper_extract_genus(file_path,azp_id)
            zoological_assign=ids_to_zoological_id_conversion_dict[azp_id]
            zoological_assignment.append(zoological_assign)
            parsed_bacteria[azp_id]=parse_gzipped_tsv_bacteria(file_path)
    elif 'UDI_54' in this_file_basename:
        print('UDI_54 processing')
        zoological_assignment.append('Small ruminant')
        helper_extract_genus(file_path,'AZP-289')
        parsed_bacteria['AZP-289']=parse_gzipped_tsv_bacteria(file_path)
    elif 'Pool' not in this_file_basename:
        azp_id_list=all_to_fill_double_extract_sample[this_file_basename]
        for azp_id in azp_id_list:
            zoological_assign=ids_to_zoological_id_conversion_dict[this_file_basename]
            zoological_assignment.append(zoological_assign)
            helper_extract_genus(file_path,azp_id)
            parsed_bacteria[azp_id]=parse_gzipped_tsv_bacteria(file_path)

        
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

eukaryotic_classifications_top3=pd.DataFrame(data=dict_for_df)


# %%
# mean and median of top genus:
print('mean top euk perc.',np.mean(eukaryotic_classifications_top3['kraken2_top1_euk_genus_perc']))
print('median top euk perc.',np.median(eukaryotic_classifications_top3['kraken2_top1_euk_genus_perc']))
# %% Correlation of top eukaryote DNA % vs pathogen DNA %
plotting_dict={}

zooscreen_pathogen_samples=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_screening/zooscreen_final/zooscreen_heatmap_merged.tsv', sep='\t',index_col=0)
for x in zooscreen_pathogen_samples.columns:
    plotting_dict[x]={'top_genus_percent':0,'top_proportion_reads_bacteria':0}
    name=x.replace('-00','-').replace('-0','-')
    this_file_parsed_bacteria=parsed_bacteria[name]
    pathogens=zooscreen_pathogen_samples[x]
    print(name,eukaryotic_classifications_top3.loc[np.where(eukaryotic_classifications_top3.Name==name)[0]][['kraken2_top1_euk_genus_perc','classified']])
    plotting_dict[x]['top_genus_percent']=eukaryotic_classifications_top3.loc[np.where(eukaryotic_classifications_top3.Name==name)[0]]['kraken2_top1_euk_genus_perc'].values[0]
    for pathogen in pathogens[pathogens].index:
        index_this_pathogen_in_kraken2=np.where(this_file_parsed_bacteria['name']==pathogen)[0]
        if len(index_this_pathogen_in_kraken2)>0:
            this_bacterial_hit=this_file_parsed_bacteria.loc[np.where(this_file_parsed_bacteria['name']==pathogen)[0]].to_dict()
            proportion_this_bacterial_hit=[int(bact_reads) for bact_reads in this_bacterial_hit['total_reads_at_or_below_node'].values()][0]/(eukaryotic_classifications_top3.loc[np.where(eukaryotic_classifications_top3.Name==name)[0]]['classified'].values[0]+eukaryotic_classifications_top3.loc[np.where(eukaryotic_classifications_top3.Name==name)[0]]['unclassified'].values[0])
            plotting_dict[x]['top_proportion_reads_bacteria']=max(plotting_dict[x]['top_proportion_reads_bacteria'],proportion_this_bacterial_hit)
            print(this_file_parsed_bacteria.loc[np.where(this_file_parsed_bacteria['name']==pathogen)[0]].to_dict())
    this_file_parsed_bacteria.loc[np.where(this_file_parsed_bacteria['name']==pathogen)[0]].to_dict()




all_data_for_plotting=pd.DataFrame.from_dict(plotting_dict).transpose()

plotting_data_from_dict=all_data_for_plotting[ (all_data_for_plotting.top_genus_percent>0.1)]

country_to_plotting_data={}
for sample_id,plotting_values in plotting_data_from_dict.iterrows():
    this_sample_country=site_to_country_mapping[labid_to_site_mapping[sample_id]]
    if this_sample_country not in country_to_plotting_data:
        country_to_plotting_data[this_sample_country]={'x':[],'y':[]}
    country_to_plotting_data[this_sample_country]['x'].append(np.log10(plotting_values['top_genus_percent']))
    country_to_plotting_data[this_sample_country]['y'].append(np.log10(plotting_values['top_proportion_reads_bacteria']))

# plotting
fig, ax = plt.subplots(facecolor='white',figsize=(6,5))
colors=['#332288','#117733','#88CCEE','#CC6677']
for country,color in zip(country_to_plotting_data,colors):
    plt.scatter(country_to_plotting_data[country]['x'],country_to_plotting_data[country]['y'],label=country,c=color)
plt.ylabel('Percentage reads assigned to top pathogen (log$_{10}$)')
plt.xlabel("Percentage reads assigned to top eukaryotic genus (log$_{10}$)")
plt.legend()
plt.savefig('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/production_figures/percentage_euk_vs_pathogen.svg',bbox_inches='tight')

X,y=plotting_data_from_dict['top_genus_percent'] ,plotting_data_from_dict['top_proportion_reads_bacteria'] 

X2 = sm.add_constant(X)
est = sm.OLS(y, X2)
est2 = est.fit()
print(est2.summary())

# %% pathology for specific animal samples?
# ONLY consider samples which were assessed for pathologies
collapsing_dict={
    'Cattle':'Bos'
}
patholgoy={'ruminants':{'pathology':0,'no_pathology':0},
           'non_ruminants':{'pathology':0,'no_pathology':0}}
totals={'pathology':0,'no_pathology':0}

pathology_assessed_samples_metadata=zooscreen_samples_metadata[(zooscreen_samples_metadata['pathology_assessed']) & ~(zooscreen_samples_metadata['is_tooth']) ]
for index,zooscreen_samples_metadata_row in pathology_assessed_samples_metadata.iterrows():
    taxonomy=zooscreen_samples_metadata_row['Taxonomic classification']
    # correct cattle --> bos
    if taxonomy in collapsing_dict:
        taxonomy=collapsing_dict[taxonomy]
    if taxonomy not in patholgoy:
        patholgoy[taxonomy]={'pathology':0,'no_pathology':0}
    if zooscreen_samples_metadata_row['Palaeopathology'] != 'No' :
        patholgoy[taxonomy]['pathology']+=1
        totals['pathology']+=1
    else:
        patholgoy[taxonomy]['no_pathology']+=1
        totals['no_pathology']+=1
    # collapse ruminants
    if taxonomy in ['Ovis','Ovis/Capra','Bos','Small ruminant','Cervus','Capra','Gazelle','Bos/Cervus','Gazella']:
        if zooscreen_samples_metadata_row['Palaeopathology'] != 'No' :
            patholgoy['ruminants']['pathology']+=1
        else:
            patholgoy['ruminants']['no_pathology']+=1
    else:
        print(taxonomy)
        if zooscreen_samples_metadata_row['Palaeopathology'] != 'No' :
            patholgoy['non_ruminants']['pathology']+=1
        else:
            patholgoy['non_ruminants']['no_pathology']+=1

p_values_pathology_species_comparisons=[]
for site in patholgoy:
    chisq_contingency_dict={site:patholgoy[site],'other_sites':{'pathology':0,'no_pathology':0}}
    chisq_contingency_dict['other_sites']['pathology']=totals['pathology']-patholgoy[site]['pathology']
    chisq_contingency_dict['other_sites']['no_pathology']=totals['no_pathology']-patholgoy[site]['no_pathology']
    chi_sq_array=np.array(pd.DataFrame.from_dict(chisq_contingency_dict).transpose())
    chi_sq_results=sts.chi2_contingency(chi_sq_array)
    p_values_pathology_species_comparisons.append(chi_sq_results.pvalue)
    print(site)
    print('Chi-Sq contigency results:',chi_sq_results)
    print('Observation:',chisq_contingency_dict)

successful_comps_fdr=multipletests(p_values_pathology_species_comparisons,alpha=0.05,method='fdr_bh')[0]

comparisons_conducted=np.array([x for x in patholgoy.keys()])
comparisons_conducted[successful_comps_fdr]

# %% Any sites with unexpectedly high/low preservation? (cutoff >1%)
cutoff_well_preserved=1

preservation_comparison_dataset={}
totals={'success':0,'fail':0}
for index,eukaryotic_classifications_top3_row in eukaryotic_classifications_top3.iterrows():
    name=eukaryotic_classifications_top3_row['Name']
    while len(name) < 7:
        name=name.replace('-','-0')
    if name in labid_to_site_mapping:
        site=labid_to_site_mapping[name]
        if site not in preservation_comparison_dataset:
            preservation_comparison_dataset[site]={'success':0,'fail':0}
        if eukaryotic_classifications_top3_row['kraken2_top1_euk_genus_perc']>cutoff_well_preserved:
            preservation_comparison_dataset[site]['success']+=1
            totals['success']+=1
        else:
            preservation_comparison_dataset[site]['fail']+=1
            totals['fail']+=1

p_values_site_comparisons_preservation=[]
for site in preservation_comparison_dataset:
    chisq_contingency_dict={site:preservation_comparison_dataset[site],'other_sites':{'success':0,'fail':0}}
    chisq_contingency_dict['other_sites']['success']=totals['success']-preservation_comparison_dataset[site]['success']
    chisq_contingency_dict['other_sites']['fail']=totals['fail']-preservation_comparison_dataset[site]['fail']
    chi_sq_array=np.array(pd.DataFrame.from_dict(chisq_contingency_dict).transpose())
    chi_sq_results=sts.chi2_contingency(chi_sq_array)
    p_values_site_comparisons_preservation.append(chi_sq_results.pvalue)
    if chi_sq_results.pvalue<(0.05/27):
        print(site)
        print('Chi-Sq contigency results:',chi_sq_results)
        print('Observation:',chisq_contingency_dict)

multipletests(p_values_site_comparisons_preservation,alpha=0.05,method='fdr_bh')
# %% site vs pathogen recovery?
pathogen_recovery={}
totals={'success':0,'fail':0}
for index,zooscreen_samples_metadata_row in zooscreen_samples_metadata.iterrows():
    site=zooscreen_samples_metadata_row['Archaeological site']
    if site not in pathogen_recovery:
        pathogen_recovery[site]={'success':0,'fail':0}
    if zooscreen_samples_metadata_row['pathogen_recovered']:
        pathogen_recovery[site]['success']+=1
        totals['success']+=1
    else:
        pathogen_recovery[site]['fail']+=1
        totals['fail']+=1

p_values_site_comparisons_recovery=[]
for site in pathogen_recovery:
    chisq_contingency_dict={site:pathogen_recovery[site],'other_sites':{'success':0,'fail':0}}
    chisq_contingency_dict['other_sites']['success']=totals['success']-pathogen_recovery[site]['success']
    chisq_contingency_dict['other_sites']['fail']=totals['fail']-pathogen_recovery[site]['fail']
    chi_sq_array=np.array(pd.DataFrame.from_dict(chisq_contingency_dict).transpose())
    chi_sq_results=sts.chi2_contingency(chi_sq_array)
    p_values_site_comparisons_recovery.append(chi_sq_results.pvalue)
    if chi_sq_results.pvalue<(0.05/27):
        print(site)
        print('Chi-Sq contigency results:',chi_sq_results)
        print('Observation:',chisq_contingency_dict)

multipletests(p_values_site_comparisons_recovery,alpha=0.05,method='fdr_bh')

# %% skeletal element vs pathogen recovery?
pathogeny_recovery_skeletal_elements={'tooth':{'success':0,'fail':0},'nontooth':{'success':0,'fail':0}}
for index,zooscreen_samples_metadata_row in zooscreen_samples_metadata.iterrows():
    element=zooscreen_samples_metadata_row['Element']
    if 'Tooth' in element:
        element='tooth'
    else:
        print(element)
        element='nontooth'
    if zooscreen_samples_metadata_row['pathogen_recovered']:
        pathogeny_recovery_skeletal_elements[element]['success']+=1
    else:
        pathogeny_recovery_skeletal_elements[element]['fail']+=1

chi_sq_array=np.array(pd.DataFrame.from_dict(pathogeny_recovery_skeletal_elements).transpose())
chi_sq_results=sts.chi2_contingency(chi_sq_array)
print('Chi-Sq contigency results:',chi_sq_results)
print('Observation:',pathogeny_recovery_skeletal_elements)


# %% Pathologies vs not for pathogen recovery?
# ONLY consider sites which were assessed for pathologies
zooscreen_samples_metadata['pathogen_recovered']=np.isin(zooscreen_samples_metadata['Plotting ID'],zooscreen_pathogen_samples.columns)
tooth_bool=[]
for x in zooscreen_samples_metadata['Element']:
    if 'Tooth' in x:
        tooth_bool.append(True)
    else:
        tooth_bool.append(False)
zooscreen_samples_metadata['is_tooth']=tooth_bool

pathology_assessed_samples_metadata=zooscreen_samples_metadata[(zooscreen_samples_metadata['pathology_assessed']) & ~(zooscreen_samples_metadata['is_tooth'])]


pathology_and_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']!='No') & (pathology_assessed_samples_metadata['pathogen_recovered']))
pathology_no_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']!='No') & ~(pathology_assessed_samples_metadata['pathogen_recovered']))
no_pathology_and_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']=='No') & (pathology_assessed_samples_metadata['pathogen_recovered']))
no_pathology_no_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']=='No') & ~(pathology_assessed_samples_metadata['pathogen_recovered']))

contingency_chisq_pathology_to_recovery=np.array([[pathology_and_recovery,pathology_no_recovery],[no_pathology_and_recovery,no_pathology_no_recovery]])
results_bone_only_pathology_to_hit=sts.chi2_contingency(contingency_chisq_pathology_to_recovery)
print('Considering ONLY bones with vs bones without paleopathological lesions',results_bone_only_pathology_to_hit)


# with teeth included

pathology_assessed_samples_metadata=zooscreen_samples_metadata[(zooscreen_samples_metadata['pathology_assessed'])]

pathology_and_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']!='No') & (pathology_assessed_samples_metadata['pathogen_recovered']))
pathology_no_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']!='No') & ~(pathology_assessed_samples_metadata['pathogen_recovered']))
no_pathology_and_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']=='No') & (pathology_assessed_samples_metadata['pathogen_recovered']))
no_pathology_no_recovery=np.sum((pathology_assessed_samples_metadata['Palaeopathology']=='No') & ~(pathology_assessed_samples_metadata['pathogen_recovered']))

contingency_chisq_pathology_to_recovery=np.array([[pathology_and_recovery,pathology_no_recovery],[no_pathology_and_recovery,no_pathology_no_recovery]])
results_bone_and_tooth_pathology_to_hit=sts.chi2_contingency(contingency_chisq_pathology_to_recovery)

print('considering also teeth',results_bone_and_tooth_pathology_to_hit)


# %% # non pathology bones vs pathology bones
pathogeny_recovery_skeletal_elements_nontooth_pathology={'nontooth_pathology':{'success':0,'fail':0},'nontooth_no_pathology':{'success':0,'fail':0}}
for index,zooscreen_samples_metadata_row in zooscreen_samples_metadata.iterrows():
    element=zooscreen_samples_metadata_row['Element']
    if 'Tooth' not in element:
        if zooscreen_samples_metadata_row['Palaeopathology'] != 'No' and zooscreen_samples_metadata_row['pathogen_recovered']:
            pathogeny_recovery_skeletal_elements_nontooth_pathology['nontooth_pathology']['success']+=1
        elif zooscreen_samples_metadata_row['Palaeopathology'] != 'No' and not zooscreen_samples_metadata_row['pathogen_recovered']:
            pathogeny_recovery_skeletal_elements_nontooth_pathology['nontooth_pathology']['fail']+=1
        elif zooscreen_samples_metadata_row['Palaeopathology'] == 'No' and zooscreen_samples_metadata_row['pathogen_recovered']:
            pathogeny_recovery_skeletal_elements_nontooth_pathology['nontooth_no_pathology']['success']+=1
        elif zooscreen_samples_metadata_row['Palaeopathology'] == 'No' and not zooscreen_samples_metadata_row['pathogen_recovered']:
            pathogeny_recovery_skeletal_elements_nontooth_pathology['nontooth_no_pathology']['fail']+=1

chi_sq_array=np.array(pd.DataFrame.from_dict(pathogeny_recovery_skeletal_elements_nontooth_pathology).transpose())
chi_sq_results=sts.chi2_contingency(chi_sq_array)
print('Chi-Sq contigency results:',chi_sq_results)
print('Observation:',pathogeny_recovery_skeletal_elements_nontooth_pathology)



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
for index,row in eukaryotic_classifications_top3.iterrows():
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

eukaryotic_classifications_top3['comparison_result']=swap_status
eukaryotic_classifications_top3['final_id']=final_id
comparison_dataset=eukaryotic_classifications_top3[['Name', 'ZooID', 'kraken2_top1_euk_genus',
       'kraken2_top1_euk_genus_perc','final_id','comparison_result','classified','unclassified']]


# getting summary stats:
number_above_tenthpercent_cutoff=len(eukaryotic_classifications_top3[(eukaryotic_classifications_top3['kraken2_top1_euk_genus_perc']>0.1)])
percentage_above_tenthpercent_cutoff=number_above_tenthpercent_cutoff/len(eukaryotic_classifications_top3)
print('number above tenth percentage cutoff for consideration:',number_above_tenthpercent_cutoff)
print('percentage above tenth percentage cutoff for consideration:',percentage_above_tenthpercent_cutoff)

agree_or_multiple=[x for x in eukaryotic_classifications_top3[(eukaryotic_classifications_top3['kraken2_top1_euk_genus_perc']>0.1)]['comparison_result'].to_list() if 'disagree' not in x]
agree_only=[x for x in agree_or_multiple if x != 'multiple -no swap']
print('number above 0.1 threshold, which agree:',len(agree_only))
print('percentage above 0.1 threshold of which agree:',len(agree_only)/number_above_tenthpercent_cutoff)

# 1% cutoff
number_above_onepercent_cutoff=len(eukaryotic_classifications_top3[(eukaryotic_classifications_top3['kraken2_top1_euk_genus_perc']>1)])
percentage_above_onepercent_cutoff=number_above_onepercent_cutoff/len(eukaryotic_classifications_top3)
print('number above one percentage cutoff for consideration:',number_above_onepercent_cutoff)
print('percentage above one percentage cutoff for consideration:',percentage_above_onepercent_cutoff)

update_number=len([x for x in eukaryotic_classifications_top3[(eukaryotic_classifications_top3['kraken2_top1_euk_genus_perc']>1)]['comparison_result'].to_list() if 'swap supported' in x])
percentage_updated=update_number/number_above_onepercent_cutoff
print('number above one percentage cutoff with host swap:',update_number)
print('percentage above one percentage cutoff with host swap:',percentage_updated)

# %%
to_clean_for_output=pd.DataFrame(data=eukaryotic_classifications_top3)

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
