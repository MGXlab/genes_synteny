import sys

# Ensure the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python process_files.py <blastoutbest_file> <prodigal_gff_file>")
    sys.exit(1)

# Get the file paths from the command line arguments
blastoutbest_file = sys.argv[1]
prodigal_gff_file = sys.argv[2]

# Extract the part before '.blastoutbest' from the blastoutbest file path
bacteria_name = blastoutbest_file.split('.blastoutbest')[0]
bacteria_name = bacteria_name.split('genomes/')[1]

# Step 1: Read the blastoutbest file and store mappings
blastoutbest_data = {}
with open(blastoutbest_file, 'r') as blast_file:
    for line in blast_file:
        parts = line.strip().split()
        protein_name = parts[0]
        protein_id = parts[1]
        
        #Store protein ID (prodigal ID) as key and name (Alvaro's name) as value
        blastoutbest_data[protein_id] = protein_name

#print(blastoutbest_data)

# Step 2: Read the prodigal GFF file and store lines indexed by ID
with open(prodigal_gff_file, 'r') as gff_file:

    for line in gff_file:
        #Ignore header
        if line.startswith('#'):
            continue
        
        #Separate columns and get 9th column
        parts = line.strip().split('\t')
        start = parts[3]
        end = parts[4] 
        strand = parts[6]
        cds = parts[2]
        bin_id = parts[0]

        attributes = parts[8]
        attr_parts = attributes.split(';')

        #Get prodigal's ID
        protein_id = attr_parts[0]
        protein_id = protein_id.replace('ID=', '')  # Remove 'ID=' prefix

        #Check if prodigal/protein ID is one of the keys in the blastoutbest dictionary 
        if protein_id in blastoutbest_data:

            protein_name = blastoutbest_data[protein_id]
            feat_id = protein_name + '_' + bacteria_name
            width = (int(end) - int(start)) + 1

            #Print line following gggenomes's GFF
            #seq_id(bacteria) bin_id(contig) start end strand feat_id(moeA_ZONMW-30) name(moaE) 
            print(bacteria_name, bin_id, start, end, strand, feat_id, protein_name, sep=',')

