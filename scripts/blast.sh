#!/bin/bash

#Loop over each protein file (proteins of interest) 
for file in proteins/*_proteins_of_interest.fasta; do

  #Extract the base name of the file (remove the extension)
  base_name=$(basename "$file" _proteins_of_interest.fasta)

  echo "$base_name"

  #Create BLAST database of Prodigal's annotated proteins
  /home/groups/VEO/tools/blast/v2.14.0/bin/makeblastdb -in proteins/${base_name}.fasta -dbtype prot -out proteins/${base_name}

  #BLAST proteins of interest against Prodigal's annotations
  /home/groups/VEO/tools/blast/v2.14.0/bin/blastp -query proteins/${base_name}_proteins_of_interest.fasta -db proteins/${base_name} -outfmt 6 -out proteins/${base_name}_proteins_of_interest_against_prodigal.blastout
done
