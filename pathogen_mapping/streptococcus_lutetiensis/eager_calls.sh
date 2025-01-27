#! /bin/bash

# these commands were used to start the eager processing of the aDNA sequencing data and modern comparison sequencing data from Forde et al 2016
##

tsv_loc=/u/iclight/eager/other_input_sheets

run_location=/ptmp/iclight/eager/mapping/zooscreen/Slutetiensis
mkdir -p ${run_location} && cd ${run_location}
grep -E 'Pool-3-038-AZP-223_S42|Pool-3-041-AZP-226_S45|Pool-4-046-AZP-283_S47|Pool-4-007-AZP-244_S8' /u/iclight/eager/screening_data/zooscreen.tsv > ${tsv_loc}/zooscreen_Slutetiensis.tsv
cat ${tsv_loc}/header.tsv ${tsv_loc}/zooscreen_Slutetiensis.tsv > ${run_location}/tsv.tmp.final.tsv
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './tsv.tmp.final.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' \
    --run_trim_bam \
    --run_bam_filtering \
    --bam_mapping_quality_threshold '0' \
    --bam_unmapped_type 'discard' \
    --bamutils_clip_double_stranded_half_udg_left '2' \
    --bamutils_clip_double_stranded_half_udg_right '2' \
    --bamutils_clip_double_stranded_none_udg_left '3' \
    --bamutils_clip_double_stranded_none_udg_right '3' 



## 

tsv_loc=/u/iclight/eager/other_input_sheets

run_location=/ptmp/iclight/eager/mapping/zooscreen/modern_Slutetiensis
mkdir -p ${run_location} && cd ${run_location}
cat ${tsv_loc}/header.tsv ${tsv_loc}/strept_modern_seqs.tsv > ${run_location}/tsv.tmp.final.tsv
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './tsv.tmp.final.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Slutetiensis_GCF900475675/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --mapper 'bwamem' \
    --skip_adapterremoval \
    --skip_fastqc \
    --skip_damage_calculation 
