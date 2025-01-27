#! /bin/bash

# these commands were used to start the eager processing of the aDNA sequencing data and modern comparison sequencing

## REFERENCE GENOME (CLADE 2)
tsv_loc=/u/iclight/eager/other_input_sheets

run_location=/ptmp/iclight/eager/mapping/zooscreen/easi241_for_erysipelothrix
mkdir -p ${run_location} && cd ${run_location}
grep -E 'EASI241' /u/iclight/eager/screening_data/200423_tartu_data.tsv > ${tsv_loc}/easi241_for_erysipelothrix.tsv
cat ${tsv_loc}/header.tsv ${tsv_loc}/easi241_for_erysipelothrix.tsv > ${run_location}/tsv.tmp.final.tsv
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './tsv.tmp.final.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta.dict' \
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


# other ancient genomes
tsv_loc=/u/iclight/eager/other_input_sheets

run_location=/ptmp/iclight/eager/mapping/zooscreen/other_ancient_seqs_for_erysipelothrix
mkdir -p ${run_location} && cd ${run_location}
cat ${tsv_loc}/header.tsv ${tsv_loc}/zooscreen_rhusiopathiae_ancient_comps.tsv > ${run_location}/tsv.tmp.final.tsv
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './tsv.tmp.final.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta.dict' \
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


######## all refseq genomes
tsv_loc=/u/iclight/eager/other_input_sheets

run_location=/ptmp/iclight/eager/mapping/zooscreen/modern_erysipelothrix_genbank
mkdir -p ${run_location} && cd ${run_location}
cat ${tsv_loc}/header.tsv ${tsv_loc}/zooscreen_rhusiopathiae_modern.tsv > ${run_location}/tsv.tmp.final.tsv
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './tsv.tmp.final.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Erysipelothrix_rhusiopathiae_GCF_900637845/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --mapper 'bwamem' \
    --skip_adapterremoval \
    --skip_fastqc \
    --skip_damage_calculation 
