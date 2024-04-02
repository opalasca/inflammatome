# Title

A snakemake pipeline for processing public RNA-Seq datasets, adapted for use on computerome (https://www.computerome.dk)

## Required input

An SRA accession txt file with the accessions of the samples of interest (download from SRA)
The name of the file should contain the dataset identifier. 

## Steps

Fasterq - download fastq files
trimmomatic - quality and adapter trimming
salmon index - build transcriptome index 
salmon quantify - estimate transcript abundance
counts - summarise counts at gene level using the tximport R package 



