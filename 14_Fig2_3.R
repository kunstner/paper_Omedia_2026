# Credentials -------------------------------------------------------------

#
# Author: Axel Künstner
# Project: Otitis media microbiome
# Data: Microbiome

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(patchwork)
library(ggpubr)

# Data --------------------------------------------------------------------

seedID <- 138
source('_helper.R')

ps_bact <- readRDS('data/bacteria_phyloseq_processed.rds')
ps_fung <- readRDS('data/fungi_phyloseq_processed.rds')

# Taxa plots --------------------------------------------------------------

# Bacteria
# Test with top 2 phyla
phylum_plot_top2_bact <- create_taxa_plot(ps_bact, rank = "Phylum", group_var = "soledum",
                                          top_n = 2, other_label = "Other")
phylum_plot_top2_bact

# all phyla
phylum_plot_all1 <- create_taxa_plot(ps_bact, rank = "Phylum", group_var = "soledum",
                                     top_n = 20, other_label = "Other Phyla")

# Create genus and species plots
genus_plot1 <- create_taxa_plot(ps_bact, rank = "Genus", group_var = "soledum",
                                top_n = 20, other_label = "Other Genera")
species_plot1 <- create_taxa_plot(ps_bact, rank = "Species", group_var = "soledum",
                                  top_n = 20, other_label = "Other Species")

phylum_plot_all1 +
    genus_plot1 + species_plot1 +
    plot_annotation(
        title = "Bacteria: Taxonomic Composition",
        subtitle = "Comparing microbial composition before and after treatment",
        tag_levels = 'A'
    ) 
ggsave("plots/Fig2.pdf", width = 15, height = 9)

# Fungi
phylum_plot_all2 <- create_taxa_plot(ps_fung, rank = "Phylum", group_var = "soledum")
genus_plot2 <- create_taxa_plot(ps_fung, rank = "Genus", group_var = "soledum",
                                top_n = 20, other_label = "Other Genera")
species_plot2 <- create_taxa_plot(ps_fung, rank = "Species", group_var = "soledum",
                                  top_n = 20, other_label = "Other Species")

# Display the plots

phylum_plot_all2 +
    genus_plot2 + species_plot2 +
    plot_annotation(
        title = "Fungi: Taxonomic Composition",
        subtitle = "Comparing microbial composition before and after treatment",
        tag_levels = 'A'
    ) 
ggsave("plots/Fig3.pdf", width = 15, height = 9)
