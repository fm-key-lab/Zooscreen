#! /bin/bash

# generate for each fastq file, the total number of reads in the file
cd /ptmp/iclight/200423_tartu_data || exit
ls -1 *.fastq.gz | while read file
    do
    total_lines=$(zcat $file | wc -l)
    total_reads=$(expr $total_lines / 4)
    echo $file,$total_reads >> ~/zooscreen/fastq_read_counts.csv
done

cd /ptmp/iclight/210617_utartu_easiagenomics_batch1 || exit
ls -1 *.fastq.gz | while read file
    do
    total_lines=$(zcat $file | wc -l)
    total_reads=$(expr $total_lines / 4)
    echo $file,$total_reads >> ~/zooscreen/fastq_read_counts.csv
done

cd /ptmp/iclight/zooscreen_data || exit
ls -1 *.fastq.gz | while read file
    do
    total_lines=$(zcat $file | wc -l)
    total_reads=$(expr $total_lines / 4)
    echo $file,$total_reads >> ~/zooscreen/fastq_read_counts.csv
done
