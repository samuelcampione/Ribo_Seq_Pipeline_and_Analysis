#!/bin/bash

# check if the input argument is provided
if [ $# -eq 0 ]; then
    echo "No file path provided. Exiting."
    exit 1
fi

input_dir=$(dirname "$1") # get the directory name from input file path
cd $input_dir

while IFS= read -r sra; do # loop through each line of input file with 1 SRA per line
    echo "Ribo-seq Analysis of ${sra}"

    
    ########################
    # Step 1 
    ########################
    echo "1. Aquire raw data"

    mkdir -p $sra
    cd $sra || exit
    
    if [ ! -f "${sra}.fastq" ]; then
        fasterq-dump --threads 8 --progress ${sra}
    else
        echo "${sra}.fastq already exists. Skipping download."
    fi

    echo "Running quality check on ${sra}.fastq"
    fastqc -q -t 8 ${sra}.fastq
    
    
    ########################
    # Step 2
    ########################
    echo "2. Trim adapters, remove duplicates, trim unique barcodes"

    # Trim adapter (given in publication)
    cutadapt --report minimal \
             -a CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             --trimmed-only -q 20 \
             --cores 8 \
             -o trimmed_adapter_${sra}.fastq \
             ${sra}.fastq
    
    # Deduplicate
    seqkit rmdup -s trimmed_adapter_${sra}.fastq -o deduplicated_${sra}.fastq

    # Remove unique barcodes
    cutadapt --report minimal --cores 8 -q 20 --minimum-length 20 -u -4 -o trimmed_${sra}.fastq deduplicated_${sra}.fastq

    echo "Running quality check report on trimmed_${sra}.fastq"
    fastqc -q -t 8 trimmed_${sra}.fastq


    ########################
    # Step 3
    ########################
    echo "3. Alignment of ribosomal RNA (output non-rRNA reads)"

    echo "Aligning against custom_rRNA_MN880429_index"
    bowtie2 -x /Users/scampione/data/custom_rRNA_MN880429.1_index/custom_rRNA_index \
            -p 8 \
            -U trimmed_${sra}.fastq \
            -S custom_rrna_alignment_${sra}.sam \
            --un custom_non_rrna_reads_${sra}.fastq
        
    echo "Aligning against custom_rRNA_AH005487_index"
    bowtie2 -x /Users/scampione/data/custom_rRNA_AH005487.2_index/custom_rRNA_index_2 \
            -p 8 \
            -U custom_non_rrna_reads_${sra}.fastq \
            -S custom_rrna_alignment_2_${sra}.sam \
            --un custom_non_rrna_reads_2_${sra}.fastq
    
    echo "Aligning against custom_rRNA_8BN3_2_index"
    bowtie2 -x /Users/scampione/data/custom_rRNA_8BN3_2_index/custom_rRNA_index_3 \
            -p 8 \
            -U custom_non_rrna_reads_2_${sra}.fastq \
            -S custom_rrna_alignment_3_${sra}.sam \
            --un custom_non_rrna_reads_3_${sra}.fastq
        
    echo "Aligning against custom_rRNA_5S_index"
    bowtie2 -x /Users/scampione/data/custom_rRNA_5S_index/custom_5S_rRNA_index \
            -p 8 \
            -U custom_non_rrna_reads_3_${sra}.fastq \
            -S custom_rrna_alignment_4_${sra}.sam \
            --un custom_non_rrna_reads_4_${sra}.fastq
    
    echo "Aligning against custom_rRNA_18S_rRNA_8CIV_c_index"
    bowtie2 -x /Users/scampione/data/custom_rRNA_18S_rRNA_8CIV_c_index/custom_18S_rRNA_index \
            -p 8 \
            -U custom_non_rrna_reads_4_${sra}.fastq \
            -S custom_rrna_alignment_5_${sra}.sam \
            --un custom_non_rrna_reads_5_${sra}.fastq

    echo "Aligning against custom_rRNA_5S_rRNA_7V08_3_index"
    bowtie2 -x /Users/scampione/data/custom_rRNA_5S_rRNA_7V08_3_index/custom_5S_rRNA_7V08_3_index \
            -p 8 \
            -U custom_non_rrna_reads_5_${sra}.fastq \
            -S custom_rrna_alignment_6_${sra}.sam \
            --un custom_non_rrna_reads_6_${sra}.fastq
        
    echo "Aligning against custom_contaminant_index"
    bowtie2 -x /Users/scampione/data/custom_contaminant_index/custom_contaminant_index \
            -p 8 \
            -U custom_non_rrna_reads_6_${sra}.fastq \
            -S custom_rrna_alignment_7_${sra}.sam \
            --un final_custom_non_rrna_reads_${sra}.fastq

    echo "Running quality check report on final_custom_non_rrna_reads_${sra}.fastq"
    fastqc -q -t 8 final_custom_non_rrna_reads_${sra}.fastq


    ########################
    # Step 4
    ########################
    echo "4. Alignment and qunatification of non rRNA reads"
    echo "Aligning against s_cerevisiae_R64 using STAR"

    STAR --runThreadN 8 \
         --genomeDir /Users/scampione/data/s_cerevisiae_R64_release_95 \
         --readFilesIn final_custom_non_rrna_reads_${sra}.fastq \
         --outFileNamePrefix ${sra}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts


    echo "End of analysis of ${sra}"
    
    cd ..

done < "$1"