A snakemake pipeline for processing public RNA-Seq datasets, adapted for use on computerome (https://www.computerome.dk)
----------------------------
----------------------------

**Required input**

- An SRA accession txt file with the accessions of the samples of interest (download from SRA). The name of the file should contain the dataset identifier (e.g. PRJNA733490.txt)

**Pipeline steps**

- Fasterq - download fastq files
- trimmomatic - quality and adapter trimming
- salmon index - build transcriptome index 
- salmon quantify - estimate transcript abundance
- create tx2gene file (TO DO - include in snakemake file)   
- get counts - summarise counts at gene level using the tximport R package (TO DO - test it is working as part of snakemake)



