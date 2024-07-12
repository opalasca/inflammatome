#!/usr/bin/env python3
import gzip

fasta_file = "../data/ref/gencode.v45.transcripts.fa.gz"
# Output file for Tx2gene mapping
output_file = "../data/ref/tx2gene/Tx2gene.txt"

# Open the input gzipped FASTA file
with gzip.open(fasta_file, "rt") as fasta:
    # Open the output file for writing
    with open(output_file, "w") as out_file:
        # Write column headers
        out_file.write("tx_id\tgene_id\n")
        
        # Initialize variables to store transcript ID and gene ID
        transcript_id = None
        gene_id = None
        
        # Iterate through each line in the FASTA file
        for line in fasta:
            # Check if the line is a header (starts with '>')
            if line.startswith('>'):
                # If we have previously encountered transcript ID and gene ID, write them to the output file
                if transcript_id and gene_id:
                    out_file.write(f"{transcript_id}\t{gene_id}\n")
                # Extract transcript ID and gene ID from the header
                fields = line.strip().lstrip('>').split('|')
                transcript_id = fields[0]
                gene_id = fields[1]
            else:
                # If the line is not a header, continue to the next line
                continue
        
        # Write the last transcript ID and gene ID to the output file
        if transcript_id and gene_id:
            out_file.write(f"{transcript_id}\t{gene_id}\n")

