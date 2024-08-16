#!/bin/bash

#Loop over each .faa file (prodigal's gene/protein predictions) 
for file in ../genomes/*.faa; do

  #Extract the base name of the file (remove the .faa extension)
  base_name=$(basename "$file" .faa)
  
  echo "$base_name"

  #Create BLAST database of Prodigal's annotated proteins
  /home/groups/VEO/tools/blast/v2.14.0/bin/makeblastdb -in ../genomes/${base_name}.faa -dbtype prot -out ../genomes/${base_name}

  #BLAST Alvaro's proteins against Prodigal's annotations
  #/home/groups/VEO/tools/blast/v2.14.0/bin/tblastn -query ../genomes/${base_name}.fna -db ../genomes/${base_name} -out ../genomes/${base_name}.out
  /home/groups/VEO/tools/blast/v2.14.0/bin/blastp -query ../proteins/${base_name}.fasta -db ../genomes/${base_name} -outfmt 6 -out ../genomes/${base_name}.blastout
done

