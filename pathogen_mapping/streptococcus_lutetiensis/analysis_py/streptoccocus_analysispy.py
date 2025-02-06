##########################################################################################################################################
####################### Zooscreen Projection Streptococcus lutetiensis ###################################################################
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

os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/ak_ancient_zoonosis_project/pathogen_mapping/streptococcus_lutetiensis/analysis_py')

ref_genome_folder='/Users/ad_loris/Nextcloud/keylab/reference_genomes/Slutetiensis_GCF900475675'

sys.path.append("/Users/ad_loris/Documents/key_lab")
import analysispy_modules_devmain as apy

first_run = False
coverage_based_cutoffs = True


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
refgenome='Slutetiensis_GCF900475675'
# # %% Load data from candidate_mutation_
# =============================================================================
#py_file='candidate_mutation_table.pickle.gz'
#[quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats] = apy.read_candidate_mutation_table_pickle_gzip(py_file)
    # # %% Breakpoint: candidate_mutation_table problems
    # # currenty only called in function when counts lacks three dimensions
    # if failFlag:
    #     continue

py_file_indels='candidate_mutation_table.pickle.gz'
[quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats,indel_p,indel_depth,indel_support,indel_size] = apy.read_candidate_mutation_table_pickle_gzip(py_file_indels)

sampleName_corrections={
        'Streptococcus_lutetiensis_BR_UFV02':'Streptococcus_infantarius_BR_UFV02',
        'Streptococcus_lutetiensis_BR_UFV05':'Streptococcus_infantarius_BR_UFV05',
        'Streptococcus_lutetiensis_CMB014':'Streptococcus_equinus_CMB014',
        'Streptococcus_lutetiensis_CMB062':'Streptococcus_equinus_CMB062',
        'Pool-3-010-AZP-195_S14':'AZP-195',
        'Pool-3-002-AZP-187_S6':'AZP-187',
        'Pool-3-038-AZP-223_S42':'AZP-223',
        'Pool-3-041-AZP-226_S45':'AZP-226',
        'Pool-4-007-AZP-244_S8':'AZP-244',
        'Pool-4-046-AZP-283_S47':'AZP-283' 
    } 

for idx,n in enumerate(sampleNames):
    if n in sampleName_corrections:
        sampleNames[idx] = sampleName_corrections[n]

p_chr_indices=[0]+[np.max(np.where(p < x + 1))+1 for x in chrStarts[1:]]
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
pattern = re.compile('AZP') # string tag in sample Name to identify outgroup in staphAD data
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
pattern = re.compile('salivarius') # string tag in sample Name to identify outgroup in staphAD data
outgroup_name = np.array(list(filter(pattern.search, list(sampleNames))))
outgroup_bool = np.in1d(sampleNames , outgroup_name)
outgroup_idx=np.nonzero(outgroup_bool)[0]

# ingroup array (bool, idx) used later
ingroup_bool = np.invert(outgroup_bool)
ingroup_idx = np.nonzero(ingroup_bool)[0]

# ancient array and bool
pattern = re.compile('AZP') # string tag in sample Name to identify outgroup in staphAD data
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

# %%
[mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning
# %% 
hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
candpos_ancient = np.where( np.sum(hasmutation[:,ancient_idx], axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!
candpos_ancient_over_half = np.where( np.sum(hasmutation[:,ancient_idx], axis=1) > len(ancient_idx)*0.5 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

candpos_modern = np.where( np.sum(hasmutation[:,(  ~ancient_bool)], axis=1)>0  )[0]
#candpos_canonical = np.where( np.sum(hasmutation[:,canonical_bool], axis=1)>0  )[0]
#candpos_modern_excluding_falso = np.where( np.sum(hasmutation[:,(canonical_bool|ochrobactrum_bool)], axis=1)>0  )[0]


candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] 
candpos_modern_100del=np.intersect1d(candpos_modern,called_100_modern_ingroup)
goodpos=np.intersect1d(candpos_modern,candpos_modern_100del)

# %%
calls_outgroup = calls[:,outgroup_bool] 
ancnti = apy.major_allele(calls_outgroup) 

#TODO: check w/ felix: calls already had filters applied, so following NOTE statement wrong? # NOTE: the filter criteria (cov,qual etc.) are not applied before major allele call
ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

# %%
goodpos2useTree = goodpos #(1:1000); %trim for easier tree view; TDL called quality_positions

# get data and filter for goodpos
calls_for_treei = calls; 
calls_for_treei = calls_for_treei[ goodpos2useTree, : ]

# build sampleNames  (that passed filter, see above from sampleNames_all --> sampleNames) w/ metainfo
treesampleNamesLong = sampleNames.copy() # include all samples 


# sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )

# translate index to nucleotide
calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
# add reference nucleotide for all positions
## TODO: figure out why this is erroring (why add chrStarts to the call, isnt this what refnts is???)

## outgroup already in samples, so unnecessary (CO92)
refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
refgenome_nts_for_tree=refnt[goodpos2useTree]

calls_for_tree_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)


apy.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'100_deletion_modern_goodpos_projected_ancient')
# building tree on cluster using script: generate_s_lutetiensis_trees

# %% collect information for table (reads mapped, num variant sites, % of total variant sites)
# read information
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
for uncorrected_sample_name in sampleName_corrections:
    if f'../qualimap/{uncorrected_sample_name}.trimmed_stats' in glob.glob('../qualimap/*'):
        reads_dict[sampleName_corrections[uncorrected_sample_name]]= {'reads_mapped':get_read_info(f'../qualimap/{uncorrected_sample_name}.trimmed_stats/genome_results.txt')}
        reads_dict[sampleName_corrections[uncorrected_sample_name]]['depth']=get_depth_info(f'../qualimap/{uncorrected_sample_name}.trimmed_stats/genome_results.txt')
        reads_dict[sampleName_corrections[uncorrected_sample_name]]['breadth']=get_breadth_info(f'../qualimap/{uncorrected_sample_name}.trimmed_stats/genome_results.txt')

for sample in reads_dict:
    # calculate total number of variant sites
    sample_index=np.where(sampleNames==sample)
    reads_dict[sample]['sites_called']=np.sum(calls[:,sample_index][goodpos]<4)
    reads_dict[sample]['proportion_sites_called']=np.sum(calls[:,sample_index][goodpos]<4)/len(goodpos)
    reads_dict[sample]['expected_breadth']=1-math.e**(float(reads_dict[sample]['depth'])*-1)

pd.DataFrame(reads_dict).transpose().to_csv('streptococcus_lutetiensis_ancient_stats.csv')

# %% 
# Individual sample projection

hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
candpos_ancient = np.where( np.sum(hasmutation[:,ancient_idx], axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!
candpos_ancient_over_half = np.where( np.sum(hasmutation[:,ancient_idx], axis=1) > len(ancient_idx)*0.5 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

candpos_modern = np.where( np.sum(hasmutation[:,((~outgroup_bool & ~ancient_bool))], axis=1)>0  )[0]
#candpos_canonical = np.where( np.sum(hasmutation[:,canonical_bool], axis=1)>0  )[0]
#candpos_modern_excluding_falso = np.where( np.sum(hasmutation[:,(canonical_bool|ochrobactrum_bool)], axis=1)>0  )[0]


candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] 
candpos_modern_100del=np.intersect1d(candpos_modern,called_100_modern_ingroup)

for s_indx in np.where(ancient_bool)[0]:
    ancient_sample_to_subset_to=sampleNames[s_indx]
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


    apy.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'100_single_sample_{ancient_sample_to_subset_to}_goodpos')
