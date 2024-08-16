#This script gets as input objects compatible to gggenomes, a NEWICK formatted taxonomy tree and outputs a synteny and homology figure and a taxonomy tree figure.
#The objects can be created by jupyter notebook get_objects.ipynb

library(tibble)
library(dplyr)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)     # plot phylogenetic trees

# Open input files ########################################################################

#Open file with object containing genome's lengths, bin/species IDs and sequence/contig IDs
alv_seqs <- read.csv("objects/alv_seqs3.csv")

#Convert seq_id to strings and object to tibble object type
alv_seqs$seq_id <- as.character(alv_seqs$seq_id)
alv_seqs <- as_tibble(alv_seqs)

#Open file with object containing gene names and coordinates
alv_genes <- read.csv("objects/alv_genes3.csv") 
#Convert bin_id to strings and object to tibble object type
alv_genes$seq_id <- as.character(alv_genes$seq_id)
alv_genes <- as_tibble(alv_genes)

#Open file with object containing homologies between genomes
alv_ava <- read.csv("objects/alv_ava2.csv")
#Convert columns to strings and object to tibble object type
alv_ava$seq_id <- as.character(alv_ava$seq_id)
alv_ava$seq_id2 <- as.character(alv_ava$seq_id2)
alv_ava <- as_tibble(alv_ava)

#Open file with object containing homologies between ortholog proteins
alv_prot_ava <- read.csv("objects/alv_prot_ava2.csv") 
#Convert object to tibble object type
alv_prot_ava <- as_tibble(alv_prot_ava)

#Open file with object containing operon coordinates
alv_operons <- read.csv("objects/alv_operons.csv") 

# Plot synteny and homologies #############################################################

plot <- gggenomes(
  genes = alv_genes,
  feats = alv_operons,                                                     # operons
  seqs = alv_seqs, 
  links = alv_ava
) |>
  add_sublinks(alv_prot_ava) |>                                            # Add grey synteny
  sync() +
  geom_feat(size = 5) +                                                    # operons
  geom_feat_note(aes(label = ""), nudge_y = -.1) +                         # operons
  geom_seq() +
  geom_link(data=links(2)) +                                               # Add grey synteny
  #  add_clusters(links=cluster) +
  geom_bin_label() +                                                       # Genes label on the right side
  geom_gene(aes(fill=name)) +                                              # removes middle black line on colored genes
  #  geom_gene_tag(aes(label=name), nudge_y=0.1, check_overlap = TRUE) +   # Adds gene names on top of colored blocks
  scale_fill_brewer("Genes", palette = "Paired", na.value = "cornsilk3") +
  theme(text = element_text(size = 12))                                    # Increase axis text size

#Save figure
ggsave("figures/synteny.pdf", plot, width = 8, height = 5, units = "in", dpi = 300)

# Plot taxonomy tree ######################################################################

#Open file containing the taxonomy tree in Newick format
alv_tree <- read.tree("files/species.treefile")
t <- ggplot(alv_tree, aes(x, y)) + geom_tree() + theme_tree() +   hexpand(.35)+ geom_tiplab(align=T, size=4, linesize=.5) 

#Save tree figure
ggsave("figures/taxonomy_tree.pdf", t, width = 20, height = 40, units = "in", dpi = 300)

#Concatenate synteny and homology and tree plots
all <- t + plot + plot_layout(widths = c(15,10))

#Save combined figure
ggsave("figures/synteny_tree.jpg", all, width = 20, height = 40, units = "in", dpi = 300)
