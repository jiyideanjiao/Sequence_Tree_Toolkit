rm(list = ls(all = TRUE))

library(ape)
library(tidyverse)
library(dplyr)
library(dispRity)
library(geiger)
library(nlme)
library(mvMORPH)
library(phytools)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide the species list CSV file as an argument.")
}

species_list_file <- args[1]

# Read the tree from the input.tre file
anttree <- read.tree("input.tre")
invisible(Ntip(anttree))

# Create a data frame from the tree tip labels
anttreespp <- as.data.frame(anttree$tip.label)
colnames(anttreespp) <- c("species")

# Read the species list CSV file
data <- read.csv(species_list_file, header = FALSE, stringsAsFactors = FALSE)
colnames(data) <- c("species")

# Ensure the data is in the correct format
data <- subset(data, species != "")
valid.list <- unique(data$species)
traits <- data.frame(valid.list)
rownames(traits) <- traits$valid.list

# Print any species that are in the CSV but not in the tree
missing_species <- setdiff(data$species, anttreespp$species)
if (length(missing_species) > 0) {
  cat("Species in CSV but not in tree:", missing_species, "\n")
}

# Clean the data and tree
cleaned_data <- dispRity::clean.data(traits, anttree)
ctraits <- cleaned_data$data
ctree <- cleaned_data$tree

# Plot the cleaned tree (optional, can be removed if not needed)
# plot(ctree)
# plot(ctree, type = "fan", cex = 0.5)
# dev.off()

# Write the cleaned tree to the standard output without additional content
cat(write.tree(ctree))
