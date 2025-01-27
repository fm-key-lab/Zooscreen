#! /bin/bash
########################################################################
## 200423_tartu_data screening data for zooscreen ##
########################################################################

## to declare before running
tsv_loc=/u/iclight/eager/screening_data
tsv_file_main=tartu_screening_names_corrected.tsv
sample_file_main=210617_utartu_unique_samples.csv
num_per_loop=5
## common for all runs
name=$(echo ${tsv_file_main} | sed 's/\.tsv//g')
total_file_len=$(wc -l <${tsv_loc}/${sample_file_main})
num_loops_needed=$(( $(( ${total_file_len} + ${num_per_loop} - 1)) / ${num_per_loop}))

## running the loop
for i in $(seq 1 ${num_loops_needed}) 
    do
    ## pre_nextflow things
    ## nextflow variables
    name_of_run=step_${i}_${name}
    ## create tsv
    start_i=$(($((${i}-1))*${num_per_loop}))
    end_i=$((${i}*${num_per_loop}))

    ## change to new folder in /ptmp and start the run
    mkdir -p /ptmp/iclight/eager/screening/${name}/${i} && cd /ptmp/iclight/eager/screening/${name}/${i}
    
    run_loc=/ptmp/iclight/eager/screening/${name}/${i}
    cat ${tsv_loc}/intermediate_files/header.tsv > ${run_loc}/tsv.tmp.final.tsv
    awk -v start_i="$start_i" -v end_i="$end_i" 'NR > start_i && NR <= end_i {print}' ${tsv_loc}/${sample_file_main} > to_grep.tsv 
    cat to_grep.tsv | while read line
        do
        grep "$line" ${tsv_loc}/${tsv_file_main} >> ${run_loc}/tsv.tmp.final.tsv
    done

    ## call nextflow
    NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 -bg -resume \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './tsv.tmp.final.tsv' \
    --fasta '/u/iclight/eager/reference_genome/GRCh38.p13_genomic.fna' \
    --bwa_index '/u/iclight/eager/reference_genome/' \
    --fasta_index '/u/iclight/eager/reference_genome/GRCh38.p13_genomic.fna.fai' \
    --seq_dict '/u/iclight/eager/reference_genome/GRCh38.p13_genomic.fna.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --run_bam_filtering \
    --bam_unmapped_type 'fastq' \
    --run_metagenomic_screening \
    --metagenomic_tool 'malt' \
    --database '/nexus/posix0/MPIIB-keylab/reference_genomes/malt_clustered_database_step4' \
    --malt_min_support_mode 'reads' \
    --run_maltextract \
    --maltextract_taxon_list '/u/iclight/zooscreen/zooscreen_final_production_run.txt' \
    --maltextract_ncbifiles '/u/iclight/bin/HOPS/Resources' \
    --maltextract_destackingoff

    sleep 15m
done

