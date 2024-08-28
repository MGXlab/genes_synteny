# Visualizing Synteny and Homology of Genes

This repository contains a tutorial and scripts to create a gene synteny and homology figure for a set of species. Additionally, it contains a tutorial to create a phylogenetic tree.   

Figures 1 b and c of the paper of Escobar Doncel and collaborators can be reproduced following this tutorial. The input files are FASTA genomes in folder ```genomes``` and proteins of interest in folder ```proteins```. They contain data for 8 bacterial species.   

<p align="center">
  <img src="figures/synteny.png" alt="Alt Text" width="850"/>
</p>

# Required software

If you have access to the draco high-performance cluster and are part of the VEO group, you only need to install gggenomes and Rstudio locally in your computer for creating the synteny and homology figure. Otherwise, install the software mentioned below according to the developer's recommendations and adapt the command lines to your needs. Inkscape is a good option to do final touches to your figure if needed and should also be installed locally. Other software are already installed in draco and can be used as mentioned in this tutorial.   

For synteny and homology figure:  

- [gggenomes](https://github.com/thackl/gggenomes)
- [Rstudio](https://posit.co/download/rstudio-desktop/)
- [Prodigal](https://github.com/hyattpd/Prodigal) 
- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [python3](https://www.python.org/downloads/) 
- [Jupyter notebook](https://jupyter.org/install)
  
To create a taxonomy tree:

- [Barrnap](https://github.com/tseemann/barrnap)
- [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/installation_without_root.html)
- [iqtree](http://www.iqtree.org/)
  
For final touches:

- [Inkscape](https://wiki.inkscape.org/wiki/Installing_Inkscape)

# Synteny and homology figure

To start, log in to draco, if you are part of the VEO group, and allocate a node to work on. If you are not part of the VEO group, adapt the commands to your needs.   

```
ssh <fsu_id>@login2.draco.uni-jena.de
salloc --partition=standard 
```

Clone this repository and move to the repository folder:

```
git clone https://github.com/MGXlab/genes_synteny.git
cd genes_synteny
```

As input, you should have FASTA files with genomes (as the files in ```genomes``` in this repository). If you already have proteins in FASTA and coordinates in BED or GFF formats, you can skip the next steps.   

If you do not have protein sequences and coordinates yet, predict genes with Prodigal with script ```scripts/prodigal.sh```, as indicated below. As output, you will obtain files named ```proteins/<SPECIES_ID>.fasta``` and ```proteins/<SPECIES_ID>.gff```.    

```
#Run Prodigal
bash scripts/prodigal.sh
```

If you already have proteins of interest (as the files in ```proteins/<SPECIES_ID>_proteins_of_interest.fasta```) and need their coordinates, BLAST them against Prodigal's protein predictions with script ```scripts/blast.sh```. You will obtain as output files named ```proteins/<SPECIES_ID>_proteins_of_interest_against_prodigal.blastout```.      

```
#Run BLAST
bash scripts/blast.sh
```

Afterwads, select the proteins of interest and name them ```proteins/<SPECIES_ID>.blastoutbest```. As an example, see below the BLAST hits of protein moeA from file ```proteins/UW101_proteins_of_interest_against_prodigal.blastout``` for species UW101:   

```
moeA    1_4772  100.000 390     0       0       1       390     1       390 0.0      787
moeA    1_5127  20.667  150     83      3       179     292     5       154 0.019    34.3
moeA    1_4765  21.642  134     87      3       169     291     148     274 0.050    32.7
moeA    1_4726  27.941  136     79      5       210     326     165     300 1.1      28.5
moeA    1_1744  54.545  22      10      0       155     176     104     125 2.6      26.9
moeA    1_4824  33.333  57      35      2       219     273     17      72  2.8      26.9
moeA    1_3329  40.000  50      23      2       113     157     1226    12733.2      27.3
moeA    1_5001  44.118  34      17      1       132     165     87      118 4.7      26.6
```

Selecting the best hit by hand, you have:

```
moeA    1_4772  100.000 390     0       0       1       390     1       390 0.0      787
```

If you do this for all proteins of interest, you will end up with files ```proteins/<SPECIES_ID>.blastoutbest``` (which are given as example in this repository).  

Next, you should adapt the format of ```.blastoutbest``` files to be compatible to gggenomes using script ```scripts/get_coordinates.py``` (the output will be a file called ```objects/alv_genes_<SPECIES_ID>.csv```). The usage of this script follows below with an example file given in this repository. You can run this script for all files within a folder using script ```scripts/get_coordinates.sh``` as below.  

```
#Optional to learn the usage of the script: run get_coordinates.py just for one file
python3 scripts/get_coordinates.py proteins/UW101.blastoutbest proteins/UW101.gff > objects/alv_genes_UW101.csv
#Run get_coordinates.py for all files
bash scripts/get_coordinates.sh > objects/alv_genes1.csv
```

After this point, you produced file ```objects/alv_genes1.csv```, which is all you need to work locally. Log out of draco and work locally. Clone the repository locally and move to the ```scripts``` folder. Next, open jupyter notebook ```get_objects.ipynb```.   

```
git clone https://github.com/MGXlab/genes_synteny.git
cd genes_synteny/scripts
jupyter notebook
#After the browser opens, select file "get_objects.ipynb"
```

Run the cells of the notebook one by one to produce the objects required by gggenomes (they are also given as example in this repository, in folder ```objects```):   

- objects/alv_seqs3.csv: contains lengths of genomes
- objects/alv_ava2.csv: contains homologies between genomes
- objects/alv_prot_ava2.csv: contains homologies between ortholog proteins
- objects/alv_genes3.csv: contains the same genes of ```objects/alv_genes1.csv```, but with shortened coordinates and spacers (if genes are too far apart) for better visualization  

Object ```alv_operons.csv``` is optional for gggenomes and indicates operon coordinates. It was written by hand for the files given in this repository.    

Next, move back to the main folder and list the files in folder ```objects``` to make sure they were successfully created.

```
cd ..
ls -lh objects
```

Now you can visualize the synteny and homology with gggenomes in RStudio using script ```scripts/synteny.r```. The output will be figures ```figures/synteny.pdf``` and ```figures/synteny_tree.jpg```. To create the figures, open the R script in RStudio and run the commands.  

# Taxonomy tree

To create a taxonomy tree of your species of interest, you start with extracting 16S rRNA genes from the genomes of the bacteria using Barrnap. This can be done with script ```scripts/barrnap.sh```, which also uses BEDtools getfasta to extract the FASTA sequences (the output of Barrnap is GFF files with coordinates):

```
#Run Barrnap for all files
bash scripts/barrnap.sh
#Change files by hand to choose one copy of each species' 16S genes
vim 16S_genes/*_16S.fasta
```

It may be that more than one 16S rRNA gene is predicted, so you should choose one copy for the subsequent steps. Concatenate all 16S files for the individual species using the bash command below:

```
cat 16S_genes/*fasta > 16S_genes/all_species_16S.fasta
```

The next step is to align the 16S rRNA genes with MAFFT using the command below. The command for iqtree is also given below.

```
/home/groups/VEO/tools/mafft/v7.505/bin/mafft all_16S.fasta > all_16S.alg
/home/groups/VEO/tools/iqtree/1.6.12/bin/iqtree -s 16S_gene/all_16S.alg
```

iqtree produces a taxonomy tree in Newick format (16S_genes/species.treefile), which you could visualize in RStudio (script ```scripts/synteny.r```) or alternatively in the online tool [iToL](https://itol.embl.de/).   

Any final touches to the figures can be done using Inkscape.
