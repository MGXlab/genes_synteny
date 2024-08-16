#!/bin/bash

# Loop over each file 
for file in proteins/*.gff; do

  # Extract the base name of the file (remove the .fna extension)
  base_name=$(basename "$file" .gff)
  
  #Run script
  python3 scripts/get_coordinates.py proteins/"$base_name".blastoutbest proteins/"$base_name".gff
done

