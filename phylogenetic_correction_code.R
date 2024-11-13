library(ape)
library(caper)
library(geiger)
library(phylolm)

bootstrap_trees <- read.nexus("birdtree_output.nex")
consensus_tree <- consensus(bootstrap_trees)
plot(consensus_tree)

s_tree <- drop.tip(consensus_tree, setdiff(consensus_tree$tip.label, s_plot$species_ID))
s_tree <- root(s_tree, outgroup = "phyllergates_cucullatus")

comp_data <- comparative.data(s_tree, s_plot, species_ID)
