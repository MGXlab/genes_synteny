#!/bin/bash

#Activate environment
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
conda activate barrnap_v0.9

#Loop over each .fna file (genome) 
for file in genomes/*.fna; do

  #Extract the base name of the file (remove the .fna extension)
  base_name=$(basename "$file" .fna)

  #echo "$base_name"
  
  #Get 16S genes from barrnap
  barrnap genomes/${base_name}.fna > 16S_genes/${base_name}.gff

  #Get Fasta sequences
  /home/groups/VEO/tools/bedtools2/v2.31.0/bin/bedtools getfasta -s -fo 16S_genes/${base_name}.fasta -fi genomes/${base_name}.fna -bed 16S_genes/${base_name}.gff

  #Copy all sequences to file, that will be processed manually
  cp 16S_genes/${base_name}.fasta 16S_genes/${base_name}_16S.fasta

done
