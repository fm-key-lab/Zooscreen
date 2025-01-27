##########################################################################################################################################
####################### Zooscreen Projection Erysipelothrix_rhusiopathiae ############################################################################################
##########################################################################################################################################

# %% 
## import libraries not gotten from analysispy_module
from datetime import date
from importlib import reload
import sys
import os
import re
import glob
from tokenize import group
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
import statistics
import math
from statsmodels.stats.multitest import multipletests
import datetime
import time
import upsetplot
import seaborn as sns

## declare paths for finding analysispy_modules and reference genome folder

os.makedirs('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_mapping/erysipelothrix_rhusiopathiae/analysis_py/genbank_modern_comps',exist_ok=True)
os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_mapping/erysipelothrix_rhusiopathiae/analysis_py/genbank_modern_comps')



sys.path.append("/Users/ad_loris/Documents/key_lab")
import analysispy_modules_devmain as apy
import analysispy_modules_main as apy_main
import analysispy_modules_testdev as apyt

first_run = False
coverage_based_cutoffs = True

ref_genome_folder='/Users/ad_loris/Nextcloud/keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845'
    
# %% Filtering parameters
# =============================================================================

# Much of this will need to vary according to your reference genome,
# coverage, and particular samples
# consider to add SD filter in addition to coverage. Isolates of other species that harbor plasmid present in pangenome at very high cov can otherwise meet the coverage threshold even though FP

# for finding fixed mutations between samples
filter_parameter_sample_across_sites = {\
                                        'min_average_coverage_to_include_sample': 4,
                                        'min_ancient_covg_for_projection':0.05,  # low, since aDNA
                                        'min_basecalls_to_include_sample\ minimum_percent_genome_with_cov8plus': 0.1, # remove samples that have too many undefined base (ie. N). added this filter.
}

filter_parameter_site_per_sample = {\
                                    'min_maf_for_call' : 0.90, #on individual samples, calls
                                    'min_cov_per_strand_for_call' : 2,  # off since paired end and single end data
                                    'min_qual_for_call' : 30,  #on individual samples, calls
                                    'min_cov_on_pos' : 4,
                                    'min_indel_gl_diff': 20,
                                    'min_indel_af': 0.8,
                                    'max_prop_0_covg_ancient': 0.9,
                                    'max_percentile_cov_ancient': 0.95}

filter_parameter_site_across_samples = {\
                                        'max_fraction_ambigious_samples' : 0.25, #across samples per position
                                        'min_median_coverage_position' : 2, #across samples per position
                                        'correlation_threshold_recombination' : 0.98,
                                        'distance_for_recombination_checking': 500}

indel_filtering_params={'min_covg_on_pos': 3,
                        'major_allele_freq_indel_call': 0.85,
                        'min_gl_diff': 5}

filtering_dicts=[filter_parameter_sample_across_sites,filter_parameter_site_per_sample,filter_parameter_site_across_samples,indel_filtering_params]

#FQ_cutoff=60; #min -FQ that samples supporting both alleles must have
#min_average_coverage_to_include_sample = 8;
#min_maf_for_call = .85; #on individual samples, calls
#min_cov_per_strand_for_call = 2;  #on individual samples, calls # off since paired end and single end data
#min_qual_for_call = 30;  #on individual samples, calls
#max_fraction_ambigious_samples = .25; #across samples per position
#min_median_coverage_position = 3; #across samples per position

## how far upstream of the nearest gene to annotate something a promoter
## mutation (not used if no annotation)
promotersize=250;

# %% 
# Define variables to store subject-specific results 
# =============================================================================
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'

para_evo_cand = np.array([],dtype=object)
mut_genes_df = {}
snp_freq_shifts_anno_lod = []
regression_res_lod = [] # stores all the results of the regression, and turned to pandas DF and saved
annotation_mutation_allParaSignal = {} # store all annotation_mutation for all detected candidates
mol_clock_data_dc = {} # dc to store for each subject the calculated data for mol clock estimation
allP_gene_paraevo = {} # store all para evo candidates that fullfill bonferroni

######################################################
### SETUP DONE ###### SETUP DONE ###### SETUP DONE ###
######################################################

######################################################
## START PROCESSING DATA #### START PROCESSING DATA ##
######################################################
# %% 
# Read in genome information
# =============================================================================
[chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder);
refgenome=ref_genome_folder.split('/')[-1]
# # %% Load data from candidate_mutation_
# =============================================================================
#py_file='candidate_mutation_table.pickle.gz'
#[quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats] = apy.read_candidate_mutation_table_pickle_gzip(py_file)
    # # %% Breakpoint: candidate_mutation_table problems
    # # currenty only called in function when counts lacks three dimensions
    # if failFlag:
    #     continue

py_file_indels=f'candidate_mutation_table.pickle.gz'
[quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats,indel_p,indel_depth,indel_support,indel_size] = apy.read_candidate_mutation_table_pickle_gzip(py_file_indels)

sampleName_corrections={
        'EASI241':'AZP-12' # ADJUSTED TO 284 due to index mislabelling
    } 

for idx,n in enumerate(sampleNames):
    if n in sampleName_corrections:
        sampleNames[idx] = sampleName_corrections[n]

p_chr_indices=[0]+[np.max(np.where(p < x + 1))+1 for x in chrStarts[1:]]

# %% remove redundant samples
unique_sample_ends={}
unique_samples_to_use=[]
for s in sampleNames:
    if s.startswith('G'): 
        tail = s.split('_')[1]
        if tail not in unique_sample_ends:
            unique_sample_ends[tail]=[s]
        else: 
            unique_sample_ends[tail].append(s)
    else:
        unique_samples_to_use.append(s)

for s in unique_sample_ends:
    unique_samples_to_use.append(unique_sample_ends[s][0])

# outgroup too distant for use in tree
pattern = re.compile('009906265.1') 
outgroup_name = np.array(list(filter(pattern.search, list(sampleNames))))
outgroup_bool = np.in1d(sampleNames , outgroup_name)

unique_samples=np.isin(sampleNames,unique_samples_to_use)

# remove suppressed records from NCBI:
pattern = re.compile('(?:GCA_029073425|GCF_029073305|GCA_029102595.1|GCA_029102585.1|GCA_029102545.1|GCA_029102525.1|GCA_029102535.1|GCA_029074585.1|GCA_029074515.1|GCA_029074505.1|GCA_029074475.1|GCA_029074465.1|GCA_029074445.1|GCA_029074425.1|GCA_029074405.1|GCA_029074385.1|GCA_029074345.1|GCA_029074355.1|GCA_029074315.1|GCA_029074285.1|GCA_029074305.1|GCA_029074245.1|GCA_029074265.1|GCA_029074215.1|GCA_029074205.1|GCA_029074185.1|GCA_029074155.1|GCA_029074145.1|GCA_029074095.1|GCA_029074105.1|GCA_029074085.1|GCA_029074065.1|GCA_029074035.1|GCA_029074005.1|GCA_029073995.1|GCA_029073985.1|GCA_029073945.1|GCA_029073955.1|GCA_029073915.1|GCA_029073895.1|GCA_029073885.1|GCA_029073865.1|GCA_029073845.1|GCA_029073815.1|GCA_029073805.1|GCA_029073785.1|GCA_029073765.1|GCA_029073745.1|GCA_029073695.1|GCA_029073715.1|GCA_029073685.1|GCA_029073665.1|GCA_029073645.1|GCA_029073625.1|GCA_029073595.1|GCA_029073585.1|GCA_029073565.1|GCA_029073545.1|GCA_029073525.1|GCA_029073505.1|GCA_029073475.1|GCA_029073455.1|GCA_029073445.1|GCA_029073425.1|GCA_029073405.1|GCA_029073365.1|GCA_029073345.1|GCA_029073355.1|GCA_029073325.1|GCA_029073305.1|GCA_029073265.1|GCA_029073275.1|GCA_029073175.1|GCA_029073165.1|GCA_029073145.1|GCA_029073125.1|GCA_029073105.1|GCA_029073065.1|GCA_029073055.1|GCA_029073045.1|GCA_029073025.1|GCA_029073005.1|GCA_029073235.1|GCA_029073215.1|GCA_029073205.1|GCA_029074545.1|GCA_029074555.1|GCA_029072975.1|GCA_029072965.1|GCA_029072945.1|GCA_029072925.1|GCA_029072905.1|GCA_029072885.1|GCA_029072845.1|GCA_029072765.1|GCA_029072745.1|GCA_029072705.1|GCA_029072715.1|GCA_029072675.1|GCA_029072665.1|GCA_029072645.1|GCA_029072615.1|GCA_029072605.1|GCA_029072555.1|GCA_029072565.1|GCA_029072545.1|GCA_029072525.1|GCA_029072505.1|GCA_029072485.1|GCA_029072465.1|GCA_029072445.1|GCA_029072425.1|GCA_029072405.1|GCA_029072385.1|GCA_029072365.1|GCA_029072345.1|GCA_029072325.1|GCA_029072295.1|GCA_029072285.1|GCA_029072265.1|GCA_029072235.1|GCA_029072225.1|GCA_029072205.1|GCA_029072135.1|GCA_029072155.1|GCA_029072145.1|GCA_029072125.1|GCA_029072095.1|GCA_029072085.1|GCA_029072045.1|GCA_029072035.1|GCA_029072025.1|GCA_029072005.1|GCA_029075055.1|GCA_029075045.1|GCA_029075015.1|GCA_029075005.1|GCA_029074955.1|GCA_029074965.1|GCA_029074945.1|GCA_029074915.1|GCA_029074905.1|GCA_029074885.1|GCA_029074835.1|GCA_029074865.1|GCA_029074825.1|GCA_029074805.1|GCA_029074785.1|GCA_029074765.1|GCA_029074745.1|GCA_029074715.1|GCA_029074705.1|GCA_029074685.1|GCA_029074625.1|GCA_029074635.1|GCA_029074615.1|GCA_029074605.1|GCA_029077745.1|GCA_029077715.1|GCA_029077705.1|GCA_029077685.1|GCA_029077665.1|GCA_029077635.1|GCA_029077625.1|GCA_029077605.1|GCA_029077585.1|GCA_029077565.1|GCA_029077545.1|GCA_029077525.1|GCA_029077505.1|GCA_029077485.1|GCA_029077455.1|GCA_029077445.1|GCA_029077425.1|GCA_029077405.1|GCA_029077385.1|GCA_029072785.1|GCA_029072855.1|GCA_029072805.1|GCA_029072815.1)') 
suppressed_name = np.array(list(filter(pattern.search, list(sampleNames))))
suppressed_bool = np.in1d(sampleNames , suppressed_name)

unique_samples_to_retain=(~outgroup_bool & unique_samples & ~suppressed_bool)

sampleNames=sampleNames[unique_samples_to_retain]
counts = counts[unique_samples_to_retain , : , : ] # keep only level (samples) that fullfil filter!
quals = quals[ : , unique_samples_to_retain ]
indel_counter = indel_counter[unique_samples_to_retain,:,:]

# %% 
# Mean cov per samples based on all positions in counts (aka p)
# =============================================================================

a = np.column_stack( (sampleNames,np.mean( np.sum(counts,axis=1),axis=1)) )
mean_cov_p_per_sample = np.array(['sample','mean_p_cov'])

mean_cov_p_per_sample = np.vstack( (mean_cov_p_per_sample,a) )    
np.savetxt('./mean_cov_per_sample.csv', mean_cov_p_per_sample, delimiter=',',fmt='%s')
#    # # analysis/plot /Users/u_key/Documents/mit/stapAD/readme/staphAD_figures_tables.r l.1040+

# %% 
# convert indels from pileups
# =============================================================================
# % indel counter reduction to indel count
# indel_counter: The first statistic is the number of reads (at this position and in this sample) that 
# support an indel. The second statistics is the number of reads (at this position and in this sample) that support a deletion.
# for now we need only the count for indels overall 
indel_counter = indel_counter[:,0,:].transpose()
# indel counter >> 50% indel row:1; row2 >> deletion

# %% 
# convert outputs to more permanent datatypes
# =============================================================================

sampleNames_all=np.asarray(sampleNames,dtype=object)
quals_all=-quals;
counts_all=counts;
coverage_all = counts_all.sum(axis=1).transpose() # axis=1 == rows; transpose needed > rows: pos and col: samples
indel_depth_all=indel_depth
indel_support_all=indel_support
indel_size_all=indel_size

# %%
pattern = re.compile('(?:ldo|EASI)') # string tag in sample Name to identify outgroup in staphAD data
ancient_name_all = np.array(list(filter(pattern.search, list(sampleNames))))
ancient_bool_all = np.in1d(sampleNames , ancient_name_all)

goodsamples_modern =  ((np.all([coverage_all.mean(axis=0) >= filter_parameter_sample_across_sites['min_average_coverage_to_include_sample']],axis=0) & ~ancient_bool_all))
goodsamples_ancient =  ((np.all([coverage_all.mean(axis=0) >= filter_parameter_sample_across_sites['min_ancient_covg_for_projection']],axis=0) & ancient_bool_all))

goodsamples = ((goodsamples_modern | goodsamples_ancient))

#Breakpoint: Too few samples passed filter, checking that at least 2 samples pass QC
if np.sum(goodsamples) < 2:
    print("Too few samples fullfill filter criteria! >> skip: " + refgenome)

sampleNames = sampleNames_all[goodsamples]
counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
quals = quals_all[ : , goodsamples ]
coverage = coverage_all[ : ,goodsamples]
indels = indel_counter[:,goodsamples]

#indel_depth=indel_depth_all[:,goodsamples,:]
#indel_support=indel_support_all[:,goodsamples]
#indel_size=indel_size_all[:,goodsamples]
#indel_total_depth=np.nansum(indel_depth,axis=2)
num_samples = len(sampleNames)

coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()
   


# %% 
# Extract refnt and define out/in-group bools
# =============================================================================

## Note ancnti/outs_nti defined below after filtered calls has been generated!
## get reference allele for all p; NOTE: analysis.m stored ATCG as double-digit-numeric
# use ref nt for mutation calls. important if multiple outgroups called 
refnt = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p,chrStarts));
refnti = apy.nts2idx(refnt)
refnti_m = np.tile(refnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele
# print(np.unique(refnt)) # sanity check

# When no outgroup defined: refnt ~= ancnt:
#ancnt = refnt   

## Estimate outgroup (ancestral allele) from ALL samples added as outgroup to SM pipeline (ancnti* == major allele)
##  NOTE: NOTE: NOTE:

## TODO: update how this is done
##  NOTE: 
# CAN ALSO JUST GRAB OUTGROUP BY OUTPUTTING THE SAMPLES CSV from case
pattern = re.compile('GCF_000373785.1') # string tag in sample Name to identify outgroup in staphAD data
outgroup_name = np.array(list(filter(pattern.search, list(sampleNames))))
outgroup_bool = np.in1d(sampleNames , outgroup_name)
outgroup_idx=np.nonzero(outgroup_bool)[0]

# ingroup array (bool, idx) used later
ingroup_bool = np.invert(outgroup_bool)
ingroup_idx = np.nonzero(ingroup_bool)[0]

# ancient array and bool
pattern = re.compile('(?:ldo|EASI)') # string tag in sample Name to identify outgroup in staphAD data
ancient_name = np.array(list(filter(pattern.search, list(sampleNames))))
ancient_bool = np.in1d(sampleNames , ancient_name)
ancient_idx=np.nonzero(ancient_bool)[0]
modern_idx=np.nonzero(~ancient_bool)[0]

# ingroup modern
ingroup_modern_bool = ( ~ancient_bool & ingroup_bool )
ingroup_modern_idx = np.nonzero(ingroup_modern_bool)
# %%
contig_positions = apy.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
[maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 
# NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
# NOTE: maf==0 -> no data;minorAF==0 -> no minor allele/ or no major allele; NT number corresponds to index in NTs [ATCG] or if maf==0 > NA == 4  
    
#  Make some basic structures for finding mutations
mutantAF = np.zeros(maNT.shape)
mutantAF[maNT != refnti_m] = maf[ maNT != refnti_m]; 
## generate mutantAF --> all positions that are not the reference, fetch major allele frequency, if reference at this pos, put zero

# mutantAF[ (minorNT != refnti_m) & (minorNT != 4) ] = mutantAF[  (minorNT != refnti_m) & (minorNT != 4)] + minorAF[  (minorNT != refnti_m) & (minorNT != 4) ] #this construction allows for positions with two different non-ancestral-mutations (when analysing data of more than one colony...plate sweeps)   

# Define mutations we do not trust in each and across samples.
# goodpos are indices of p that we trust
## Filter per mutation
calls = maNT
failed_quals = (quals[:,~ancient_bool] < filter_parameter_site_per_sample['min_qual_for_call'])
failed_maf=(maf[:,~ancient_bool] < filter_parameter_site_per_sample['min_maf_for_call'])
failed_forward=(coverage_forward_strand[:,~ancient_bool] < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
failed_reverse=(coverage_reverse_strand[:,~ancient_bool] < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
failed_cov=(coverage_reverse_strand[:,~ancient_bool] + coverage_forward_strand[:,~ancient_bool] < filter_parameter_site_per_sample['min_cov_on_pos'])
failed_indels=(indels[:,~ancient_bool] > (0.5*coverage[:,~ancient_bool]) )

calls_for_modern=calls[:,~ancient_bool]
calls_for_modern[( failed_quals | failed_maf | failed_forward | failed_reverse | failed_cov | failed_indels  )] = 4
calls[:,~ancient_bool]=calls_for_modern

failed_maf_all=(maf < filter_parameter_site_per_sample['min_maf_for_call'])
calls[failed_maf_all]=4
# NOTE: reassigning calls to 4 changes also maNT!

## Filter per site across samples
# Ignore here outgroup samples!
called_100_modern_ingroup=np.where(np.sum(calls[:,((ingroup_modern_bool | outgroup_bool ))]<4,axis=1) >= (np.sum((ingroup_modern_bool | outgroup_bool ))))[0]

# ID any sites to set to N in ancient samples
filter_parameter_site_per_projection = {\
                                    'min_maf_for_call' : 0.90, 
                                    'min_qual_for_call' : 30,  
                                    'min_cov_on_pos' : 1}

projection_failed_quals = (quals[:,ancient_bool] < filter_parameter_site_per_projection['min_qual_for_call'])
projection_failed_maf=(maf[:,ancient_bool] < filter_parameter_site_per_projection['min_maf_for_call'])
projection_failed_cov=(coverage_reverse_strand[:,ancient_bool] + coverage_forward_strand[:,ancient_bool] < filter_parameter_site_per_projection['min_cov_on_pos'])

calls_for_projection=calls[:,ancient_bool]
calls_for_projection[( projection_failed_quals | projection_failed_maf | projection_failed_cov )] = 4
calls[:,ancient_bool]=calls_for_projection


# %% QC Plotting
#################################################
## Covered positions histogram ancient samples ##
#################################################
plt.hist(np.sum(coverage[:,ancient_bool]>0,axis=1))
plt.title('Positions with any coverage across ancient sample')
plt.xlabel('Number of ancient sample with coverage')
plt.ylabel('Number of positions')
plt.savefig('histogram_positions_covered_ancient.png',bbox_inches='tight')

#################################################
## Called positions histogram ancient samples ##
#################################################
plt.hist(np.sum((calls[:,ancient_bool]<4),axis=1))
plt.title('Positions with call across ancient samples')
plt.xlabel('Number of ancient sample with call')
plt.ylabel('Number of positions')
plt.savefig('histogram_positions_called_ancient.png',bbox_inches='tight')

plt.hist(np.sum((calls[:,ancient_bool][called_100_modern_ingroup]<4),axis=1),bins=85)
plt.title('100 deletion Modern-ingroup positions with call across ancient sample')
plt.xlabel('Number of ancient samples with call')
plt.ylabel('Number of positions')
plt.savefig('histogram_positions_called_100del_ancient.png',bbox_inches='tight')



#######################################################
## Positions covered overlap with modern ref genomes ##
#######################################################
# calculate proporiton of positions overlapping with each modern, for each ancient
summary_table={}
max_table={}
for a_index in ancient_idx:
    ancient_name=sampleNames[a_index]
    summary_table[ancient_name]={}
    current_max=-1
    for m_index in modern_idx:
        modern_name=sampleNames[m_index]
        overlapping_positions=np.sum((calls[:,m_index] == calls[:,a_index]) & ( calls[:,a_index] <  4 ))
        non_overlapping_positions=np.sum((calls[:,m_index] != calls[:,a_index]) & ( calls[:,a_index]<  4 ) &  ( calls[:,m_index] < 4 ))
        diff_overlap_non_overlap=overlapping_positions-non_overlapping_positions
        summary_table[ancient_name][modern_name]=(diff_overlap_non_overlap)/np.sum((calls[:,m_index] < 4 ) & ( calls[:,a_index] < 4))
        if (diff_overlap_non_overlap)/np.sum((calls[:,m_index] < 4 ) & ( calls[:,a_index] < 4)) > current_max:
            current_max = (diff_overlap_non_overlap)/np.sum((calls[:,m_index] < 4 ) & ( calls[:,a_index] < 4))
            max_table[ancient_name] = modern_name

summary_table_df=pd.DataFrame(summary_table)
plt.subplots(figsize=(20,10))
sns.heatmap(summary_table_df)
plt.title('(Positions agree - Positions disagree)/Positions called in both')
plt.savefig('heatmap_positions_agreement.png',bbox_inches='tight')

summary_table={}
max_table={}
for a_index in ancient_idx:
    ancient_name=sampleNames[a_index]
    summary_table[ancient_name]={}
    current_max=-1
    for m_index in modern_idx:
        modern_name=sampleNames[m_index]
        overlapping_positions=np.sum((calls[:,m_index] == calls[:,a_index]) & ( calls[:,a_index] <  4 ))
        non_overlapping_positions=np.sum((calls[:,m_index] != calls[:,a_index]) & ( calls[:,a_index]<  4 ) &  ( calls[:,m_index] < 4 ))
        diff_overlap_non_overlap=overlapping_positions-non_overlapping_positions
        summary_table[ancient_name][modern_name]=diff_overlap_non_overlap
        if (diff_overlap_non_overlap)/np.sum((calls[:,m_index] < 4 ) & ( calls[:,a_index] < 4)) > current_max:
            current_max = (diff_overlap_non_overlap)/np.sum((calls[:,m_index] < 4 ) & ( calls[:,a_index] < 4))
            max_table[ancient_name] = modern_name


summary_table_df=pd.DataFrame(summary_table)
plt.subplots(figsize=(20,10))
sns.heatmap(summary_table_df)
plt.title('Positions compared')

# %%
[mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning
# %% 
hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!

pattern = re.compile('(?:ldo050)') # string tag in sample Name to identify outgroup in staphAD data
ldo050_index = np.array(list(filter(pattern.search, list(sampleNames))))
ldo050_bool = np.in1d(sampleNames , ldo050_index)

candpos_modern_ldo = np.where( np.sum(hasmutation[:,(~ancient_bool | ldo050_bool)], axis=1)>0  )[0]

goodpos=np.intersect1d(candpos_modern_ldo,called_100_modern_ingroup)

def get_read_info(qualimap_path):
    with open(qualimap_path) as f:
        for l in f.readlines():
            if l.strip().startswith('number of reads'):
                result=int(l.strip().split('number of reads = ')[1].replace(',',''))
                break
    return result
def get_depth_info(qualimap_path):
    with open(qualimap_path) as f:
        for l in f.readlines():
            if l.strip().startswith('mean coverageData '):
                coverage_result=l.strip().split('mean coverageData = ')[1].replace('X','')
                break
    return coverage_result

def get_breadth_info(qualimap_path):
    with open(qualimap_path) as f:
        for l in f.readlines():
            if l.strip().startswith('There is a '):
                coverage_result=l.strip().split('There is a ')[1].split('%')[0]
                break
    return coverage_result


reads_dict={}
for uncorrected_sample_name in sampleNames[ancient_bool]:
    if f'../../qualimap/{uncorrected_sample_name}.trimmed_stats' in glob.glob('../../qualimap/*'):
        reads_dict[uncorrected_sample_name]= {'reads_mapped':get_read_info(f'../../qualimap/{uncorrected_sample_name}.trimmed_stats/genome_results.txt')}
        reads_dict[uncorrected_sample_name]['depth']=get_depth_info(f'../../qualimap/{uncorrected_sample_name}.trimmed_stats/genome_results.txt')
        reads_dict[uncorrected_sample_name]['breadth']=get_breadth_info(f'../../qualimap/{uncorrected_sample_name}.trimmed_stats/genome_results.txt')

for sample in reads_dict:
    # calculate total number of variant sites
    sample_index=np.where(sampleNames==sample)
    reads_dict[sample]['sites_called']=np.sum(calls[:,sample_index][goodpos]<4)
    reads_dict[sample]['proportion_informative_sites_called']=np.sum(calls[:,sample_index][goodpos]<4)/len(goodpos)
    reads_dict[sample]['expected_breadth']=1-math.e**(float(reads_dict[sample]['depth'])*-1)
    
pd.DataFrame(reads_dict).transpose().to_csv('rhusio_ancient_stats.csv')


# %%
calls_outgroup = calls[:,outgroup_bool] 
ancnti = apy.major_allele(calls_outgroup) 

#TODO: check w/ felix: calls already had filters applied, so following NOTE statement wrong? # NOTE: the filter criteria (cov,qual etc.) are not applied before major allele call
ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

# %%
goodpos2useTree = goodpos #(1:1000); %trim for easier tree view; TDL called quality_positions


# intersection of present on all canonical brucellosis, positions with mutation in ancient samples (after MAF cutoff), and >8 outgroup ones have coverage on pos (of 11))

# get data and filter for goodpos
calls_for_treei = calls; 
calls_for_treei = calls_for_treei[ goodpos2useTree, : ]
calls_for_treei_no_outgroup = calls_for_treei[ : , ingroup_bool ]

# build sampleNames  (that passed filter, see above from sampleNames_all --> sampleNames) w/ metainfo
treesampleNamesLong_no_outgroup = sampleNames[ingroup_bool] # remove outgroup samples
treesampleNamesLong = sampleNames.copy() # include all samples 

num_pos_projected_ancient=np.sum(calls[:,ancient_bool][goodpos]<4,axis=0)
for index_in_ancient_index,index in enumerate(ancient_idx):
    treesampleNamesLong[index]=treesampleNamesLong[index]


# translate index to nucleotide
calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
calls_for_tree_no_outgroup= apy.idx2nts(calls_for_treei_no_outgroup)
# add reference nucleotide for all positions
## TODO: figure out why this is erroring (why add chrStarts to the call, isnt this what refnts is???)

## outgroup already in samples, so unnecessary (CO92)
refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
refgenome_nts_for_tree=refnt[goodpos2useTree]


calls_for_tree_anc_ref_no_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree_anc[:, None]),calls_for_tree_anc_no_outgroup),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both
calls_for_tree_anc_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree_anc[:, None]),calls_for_tree_anc),axis=1)

calls_for_tree_ref_no_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree_no_outgroup),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both
calls_for_tree_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

treesampleNamesLong_ref_no_outgroup = np.append(['Sref'],treesampleNamesLong_no_outgroup) # add name for outgroup
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)


apyt.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'100_deletion_modern_goodpos_projected_ancient')


apy_main.generate_tree(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,treesampleNamesLong_ref_outgroup,refgenome,filetag=f'100_deletion_modern_ldo_goodpos_projected_ancient',buildTree='ML',additional_raxml_parameters=f'--threads 4 --bs-trees 100 --bs-metric TBE --all --redo')


# %% generate tree just for EASI tree -- allowing for SNVs called in that specific sample
hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] 

for s_indx in np.where(ancient_bool)[0]:
    ancient_sample_to_subset_to=sampleNames[s_indx]
    if ancient_sample_to_subset_to == 'EASI241':
        called_100_ancient_sample=np.where(calls[:,s_indx]<4)[0]

        # getting samples to include in THIS tree
        subsample_to_use= ~ancient_bool
        subsample_to_use[s_indx]=True

        candpos_modern_100del_called_ancient=np.intersect1d(called_100_ancient_sample,candpos_modern_100del)

        goodpos2useTree=np.intersect1d(candpos,candpos_modern_100del_called_ancient)

        calls_outgroup = calls[:,outgroup_bool] 
        ancnti = apy.major_allele(calls_outgroup) 
        ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

        # get data and filter for goodpos
        calls_for_treei = calls[:,subsample_to_use]; 
        calls_for_treei = calls_for_treei[ goodpos2useTree, : ]

        # build sampleNames  (that passed filter, see above from sampleNames_all --> sampleNames) w/ metainfo
        treesampleNamesLong = sampleNames.copy()[subsample_to_use] # include all samples 

        # translate index to nucleotide
        calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
        
        # add reference nucleotide for all positions
        ## outgroup already in samples, so unnecessary (CO92)
        refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
        refgenome_nts_for_tree=refnt[goodpos2useTree]

        calls_for_tree_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

        treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

        apyt.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'100_single_sample_{ancient_sample_to_subset_to}_goodpos')
        apy_main.generate_tree(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,treesampleNamesLong_ref_outgroup,refgenome,filetag=f'100_single_sample_{ancient_sample_to_subset_to}_goodpos',buildTree='ML',additional_raxml_parameters=f'--threads 2 --bs-trees 100 --all --redo')
