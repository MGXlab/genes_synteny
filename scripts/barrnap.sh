#!/bin/bash

#Loop over each .fna file (genome) 
for file in genomes/*.fna; do

  #Extract the base name of the file (remove the .fna extension)
  base_name=$(basename "$file" .fna)
  
  #echo "$base_name"

  #Get 16S genes from barrnap
  barrnap genomes/${base_name}.fna > 16S_genes/${base_name}.gff

  #Get Fasta sequences
  bedtools getfasta -s -fo 16S_genes/${base_name}.fasta -fi genomes/${base_name}.fna -bed 16S_genes/${base_name}.gff

  #Copy all sequences to file, that will be processed manually
  cp 16S_genes/${base_name}.fasta 16S_genes/${base_name}_16S.fasta

done
