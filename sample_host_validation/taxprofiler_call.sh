#! /bin/bash 
# nf-core/taxprofiler run

# unified samplesheet for in-house + second tartu set
nextflow run nf-core/taxprofiler -resume -r 1.1.0 -c /u/iclight/eager/mpcdf_raven.config \
    -w '/ptmp/iclight/taxprofiler/zooscreen/work' \
    -profile singularity,mpcdf,raven,resources_normal,taxprofiler \
    --input /u/iclight/taxprofiler/zooscreen/samplesheet_unified_azp.tsv --databases /u/iclight/taxprofiler/zooscreen/databases.csv \
    --outdir /ptmp/iclight/taxprofiler/zooscreen \
    --perform_shortread_qc \
    --run_kraken2 \
    --kraken2_save_minimizers

# first set tartu data
nextflow run nf-core/taxprofiler -resume -r 1.1.0 -c /u/iclight/eager/mpcdf_raven.config \
    -w '/ptmp/iclight/taxprofiler/zooscreen/work' \
    -profile singularity,mpcdf,raven,resources_normal,taxprofiler \
    --input /u/iclight/taxprofiler/zooscreen/samplesheet_tartu_first_batch.csv --databases /u/iclight/taxprofiler/zooscreen/databases.csv \
    --outdir /ptmp/iclight/taxprofiler/zooscreen_tartu_first_batch \
    --perform_shortread_qc \
    --perform_runmerging \
    --run_kraken2 \
    --kraken2_save_minimizers

