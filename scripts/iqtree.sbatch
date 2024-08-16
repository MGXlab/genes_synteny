#!/bin/bash
#SBATCH --job-name=tree_16S
#SBATCH --output=tree_16S_%j.out
#SBATCH --error=tree_16S_%j.err
#SBATCH --cpus-per-task=40
#SBATCH --hint=nomultithread
#SBATCH --partition=standard
#SBATCH --mem=240G

#Using all single ortholog genes aligned with GTDBtk
#/home/groups/VEO/tools/iqtree/1.6.12/bin/iqtree -s gtdbtk_align_output/align/gtdbtk.bac120.msa.fasta.gz -nt 40

#Using only the 16S gene predicted with Barrnap
/home/groups/VEO/tools/iqtree/1.6.12/bin/iqtree -s 16S_genes/all_16S.alg -nt 40
