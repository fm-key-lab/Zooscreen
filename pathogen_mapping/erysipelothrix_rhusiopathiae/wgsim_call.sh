#! /bin/bash
# download genomes from ncbi
datasets download genome taxon 1648

conda activate basic_genomics

num_reads=1000000
ls -1 ncbi_dataset/data | while read accession
do
    fasta_file_dir=ncbi_dataset/data/${accession}
    fastq_r1=ncbi_dataset/data/${accession}/${accession}_R1.fq
    fastq_r2=ncbi_dataset/data/${accession}/${accession}_R2.fq
    wgsim -N ${num_reads} -r 0 -e 0 -R 0 -1 100 -2 100 ${fasta_file_dir}/*.fna ${fastq_r1} ${fastq_r2}
done

# saving into appropriate format for eager
ls -1 ncbi_dataset/data/ | \
awk -F' ' '{print $1, $1, "1", "4", "PE", "NA","double","none","/ptmp/iclight/zooscreen_rhusiopathiae_seqs/ncbi_dataset/data/"$1"/"$1"_R1.fq","/ptmp/iclight/zooscreen_rhusiopathiae_seqs/ncbi_dataset/data/"$1"/"$1"_R2.fq","NA"}' | \
sed 's/ /\t/g' > ~/eager/other_input_sheets/zooscreen_rhusiopathiae_genbank_seqs.tsv
