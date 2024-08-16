#!/bin/bash

#Create proteins file
mkdir proteins

# Loop over each .fna file 
for file in genomes/*.fna; do

  # Extract the base name of the file (remove the .fna extension)
  base_name=$(basename "$file" .fna)
  
  echo "$base_name"

  # Run Prodigal with the specified options
  /home/groups/VEO/tools/prodigal/v2.6.3/prodigal -i "$file" -o proteins/"${base_name}.gff" -f gff -a proteins/"${base_name}.faa"
done

