#!/bin/bash

#Genome download:  wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#This genome seems to already be ready to be used.

#TODO: Add link for the genome.annotation

#The GTF and fasta file have different names for chr
#sed -i 's/^chr//' genome/gencode.v43.annotation.gtf 
#sed -i 's/^M\t/MT\t/' gencode.v43.annotation.gtf #For GTF

# Note: Requires tidyverse to be installed (R)

#  Alignment
#Generate index 
STAR --runMode genomeGenerate --runThreadN 32 --genomeDir genome_index --sjdbGTFfile genome/gencode.v43.annotation.gtf --sjdbOverhang 49 --genomeFastaFiles genome/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Create alignment folder 
mkdir alignment


while IFS= read -r file; do
    # Run the command with the current file
    STAR --runThreadN 32 --genomeDir genome_index --quantMode TranscriptomeSAM --readFilesIn "RawData/$file" --outFileNamePrefix "alignment/$file";
done < file_list.txt

# Get statistics from these files

# Output file for the summary
output_file="star_alignment_summary.txt"

# Print the header
echo "File\tNumber of input reads\tUniquely mapped reads %\tNumber of reads mapped to multiple loci\t% of reads mapped to too many loci\t% of reads unmapped: too short\t% of reads unmapped: too many mismatches" > $output_file

# Loop through each .Log.final.out file in the alignment directory
for file in alignment/*Log.final.out; do
    # Extract the required statistics using grep and awk
    num_input_reads=$(grep "Number of input reads" "$file" | awk '{print $6}')
    uniquely_mapped_reads_percent=$(grep "% Uniquely mapped reads" "$file" | awk '{print $6}')
    reads_mapped_to_multiple_loci=$(grep "% of reads mapped to multiple loci" "$file" | awk '{print $9}')
    percent_mapped_to_too_many_loci=$(grep "% of reads mapped to too many loci" "$file" | awk '{print $10}')
    percent_unmapped_too_short=$(grep "% of reads unmapped: too short" "$file" | awk '{print $8}')
    percent_unmapped_too_many_mismatches=$(grep "% of reads unmapped: too many mismatches" "$file" | awk '{print $9}')

    # Get the base name of the file (remove the directory and .Log.final.out extension)
    base_name=$(basename "$file" Log.final.out)

    # Print the extracted values to the summary file
    echo "$base_name\t$num_input_reads\t$uniquely_mapped_reads_percent\t$reads_mapped_to_multiple_loci\t$percent_mapped_to_too_many_loci\t$percent_unmapped_too_short\t$percent_unmapped_too_many_mismatches" >> $output_file
done

echo "Summary of STAR alignments saved to $output_file."

# Create ref folder for RSEM
mkdir ref

# Make reference for RSEM
rsem-prepare-reference genome/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa ref/GRCh38 --gtf  genome/gencode.v43.annotation.gtf #--trusted-sources ENSEMBL

# Quantify the transcript and gene
mkdir quant

while IFS= read -r file; do
    # Run the command with the current file
	sample_name="${file%%.*}"
	
	mkdir "quant/${sample_name}"
	
	rsem-calculate-expression -p 32 --seed 32 --alignments "alignment/${file}Aligned.toTranscriptome.out.bam" ref/GRCh38 "quant/${sample_name}/${sample_name}"
done < file_list.txt

# Compile the expression into a single file 
Rscript src/merge_transcript_expression.R # Requires tidyverse to be installed

## SUPPA

# Generate events
# Make dir
mkdir SUPPA

# Compute psi per isoform
python3 /home/flavia/utils/SUPPA/suppa.py psiPerIsoform -g genome/gencode.v43.annotation.gtf -e quant/tpm_isoform.csv -o SUPPA/psiPerIsoform

# Generate events
python3 /home/flavia/utils/SUPPA/suppa.py generateEvents -i genome/gencode.v43.annotation.gtf -o SUPPA/events -f ioe -e SE SS MX RI FL

# Compute psi per event
for splicing_event in SE MX AF AL A5 A3 RI; do
	splicing_event_file=SUPPA/events_${splicing_event}_strict.ioe
	output_file=SUPPA/${splicing_event}
	python3 /home/flavia/utils/SUPPA/suppa.py psiPerEvent --ioe-file $splicing_event_file --expression-file  quant/tpm_isoform.csv -o $output_file --save_tpm_events 
done

## Process data in R
# Generate a "biotype" .csv, as the analysis will focused on pcod and lincRNA genes

# Extract gene_id and gene_type from the GTF file
grep -w "gene" "genome/gencode.v43.annotation.gtf" | \
awk -F '\t' '{print $9}' | \
awk -F ';' '{ 
    for(i=1; i<=NF; i++) {
        if ($i ~ /gene_id/) { split($i, a, "\""); gene_id = a[2] }
        if ($i ~ /gene_type/) { split($i, b, "\""); gene_type = b[2] }
    }
    print gene_id "," gene_type
}' > "genome/biotype.csv"

echo "biotype.csv file created successfully."

Rscript src/compare_entropy_per_gene.R