# Credentials -------------------------------------------------------------

#
# Network Analysis Script
# Author: Axel Künstner
# Project: Otitis media microbiome
# Data: Microbiome

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(SpiecEasi)
library(igraph)
library(ggnetwork)
library(patchwork)

# Data --------------------------------------------------------------------

source('_helper.R') # Ensure the helper functions are available
seedID <- 138

ps_bact <- readRDS('data/bacteria_phyloseq_processed.rds')
ps_fung <- readRDS('data/fungi_phyloseq_processed.rds')

# Network analysis functions ----------------------------------------------

# Function to create filtered SpiecEasi networks from phyloseq objects
create_network <- function(ps, method = "spiec.easi", thresh = 0.05, 
                           min_prevalence = 0.2, min_abundance = 0, 
                           max_taxa = 50) {
    
    message("Starting with ", ntaxa(ps), " taxa")
    
    prev <- rowSums(otu_table(ps) > 0) / nsamples(ps)
    ps_rel <- microbiome::transform(ps, "compositional")
    abund <- rowSums(otu_table(ps_rel)) / nsamples(ps)
    
    taxa_filter <- prev >= min_prevalence & abund >= min_abundance
    message("After prevalence and abundance filtering: ", sum(taxa_filter), " taxa")
    
    ps_filtered <- prune_taxa(taxa_filter, ps)
    
    if (ntaxa(ps_filtered) > max_taxa) {
        top_taxa <- names(sort(abund[taxa_filter], decreasing = TRUE)[1:max_taxa])
        ps_filtered <- prune_taxa(top_taxa, ps_filtered)
        message("After limiting to top ", max_taxa, " taxa by abundance")
    }
    
    message("Final taxa count for network analysis: ", ntaxa(ps_filtered))
    
    if (phyloseq::nsamples(ps_filtered) < 5 || phyloseq::ntaxa(ps_filtered) < 5) {
        stop("Not enough samples or taxa for network inference. Need at least 5 of each.")
    }
    
    ps_rel <- microbiome::transform(ps_filtered, "compositional")
    
    otu_table <- as.matrix(phyloseq::otu_table(ps_rel))
    if (phyloseq::taxa_are_rows(ps_rel)) {
        otu_table <- t(otu_table)
    }
    
    var_taxa <- apply(otu_table, 2, var)
    otu_table <- otu_table[, var_taxa > 0, drop = FALSE]
    
    if (ncol(otu_table) < 5) {
        stop("Not enough taxa with variance for network inference. Need at least 5.")
    }
    
    spiec_params <- list(
        method = method,
        lambda.min.ratio = 1e-2,
        nlambda = 20,
        pulsar.params = list(
            thresh = 0.05,
            subsample.ratio = 0.8,
            rep.num = 10,
            ncores = 1
        )
    )
    
    tryCatch({
        spiec_net <- SpiecEasi::spiec.easi(
            data = otu_table, 
            method = spiec_params$method,
            lambda.min.ratio = spiec_params$lambda.min.ratio,
            nlambda = spiec_params$nlambda,
            pulsar.params = spiec_params$pulsar.params,
            verbose = TRUE
        )
        
        adj_matrix <- SpiecEasi::getRefit(spiec_net)
        adj_matrix[abs(adj_matrix) < thresh] <- 0
        
        species_names <- colnames(otu_table)
        rownames(adj_matrix) <- species_names
        colnames(adj_matrix) <- species_names
        
        g <- igraph::graph_from_adjacency_matrix(
            adj_matrix, 
            mode = "undirected", 
            weighted = TRUE,
            diag = FALSE
        )
        
        tax_table <- as.data.frame(phyloseq::tax_table(ps))
        tax_table_filtered <- tax_table[match(species_names, rownames(tax_table)), , drop = FALSE]
        
        for (col in colnames(tax_table_filtered)) {
            if (!is.null(tax_table_filtered[[col]])) {
                igraph::V(g)$name <- species_names
                igraph::vertex_attr(g, col) <- tax_table_filtered[[col]]
            }
        }
        
        igraph::V(g)$degree <- igraph::degree(g)
        igraph::V(g)$betweenness <- tryCatch(
            igraph::betweenness(g, normalized = TRUE),
            error = function(e) {
                message("Cannot calculate betweenness: ", e$message)
                rep(0, igraph::vcount(g))
            }
        )
        igraph::V(g)$eigen_centrality <- tryCatch(
            igraph::eigen_centrality(g)$vector,
            error = function(e) {
                message("Cannot calculate eigenvector centrality: ", e$message)
                rep(0, igraph::vcount(g))
            }
        )
        
        # Use absolute weights for community detection
        g_abs <- g
        igraph::E(g_abs)$weight <- abs(igraph::E(g_abs)$weight)
        comm <- tryCatch(
            igraph::cluster_fast_greedy(g_abs),
            error = function(e) {
                message("Cannot calculate communities: ", e$message)
                igraph::make_clusters(g, membership = rep(1, igraph::vcount(g)))
            }
        )
        
        igraph::V(g)$module <- igraph::membership(comm)
        
        return(list(
            graph = g,
            adjacency = adj_matrix,
            spiec = spiec_net,
            communities = comm,
            filtered_taxa = species_names
        ))
        
    }, error = function(e) {
        message("Error in SpiecEasi network creation: ", e$message)
        
        cor_matrix <- cor(otu_table, method = "spearman")
        cor_matrix[is.na(cor_matrix)] <- 0
        cor_matrix[abs(cor_matrix) < thresh] <- 0
        diag(cor_matrix) <- 0
        
        g <- igraph::graph_from_adjacency_matrix(
            cor_matrix,
            mode = "undirected",
            weighted = TRUE,
            diag = FALSE
        )
        
        species_names <- colnames(otu_table)
        igraph::V(g)$name <- species_names
        igraph::V(g)$degree <- igraph::degree(g)
        
        tax_table <- as.data.frame(phyloseq::tax_table(ps))
        tax_table_filtered <- tax_table[match(species_names, rownames(tax_table)), , drop = FALSE]
        for (col in colnames(tax_table_filtered)) {
            if (!is.null(tax_table_filtered[[col]])) {
                igraph::vertex_attr(g, col) <- tax_table_filtered[[col]]
            }
        }
        
        # Use absolute weights for community detection
        g_abs <- g
        igraph::E(g_abs)$weight <- abs(igraph::E(g_abs)$weight)
        comm <- tryCatch(
            igraph::cluster_fast_greedy(g_abs),
            error = function(e) {
                message("Cannot calculate communities: ", e$message)
                igraph::make_clusters(g, membership = rep(1, igraph::vcount(g)))
            }
        )
        
        igraph::V(g)$module <- igraph::membership(comm)
        message("Falling back to correlation-based network")
        
        return(list(
            graph = g,
            adjacency = cor_matrix,
            spiec = NULL,
            communities = comm,
            filtered_taxa = species_names
        ))
    })
}

# Split phyloseq objects by treatment status ------------------------------

# Bacteria
ps_bact_before <- phyloseq::subset_samples(ps_bact, soledum == "before")
ps_bact_after <- phyloseq::subset_samples(ps_bact, soledum == "after")

# Remove taxa with zero counts after filtering
ps_bact_before <- phyloseq::prune_taxa(taxa_sums(ps_bact_before) > 0, ps_bact_before)
ps_bact_after <- phyloseq::prune_taxa(taxa_sums(ps_bact_after) > 0, ps_bact_after)

# Fungi
ps_fung_before <- phyloseq::subset_samples(ps_fung, soledum == "before")
ps_fung_after <- phyloseq::subset_samples(ps_fung, soledum == "after")

# Remove taxa with zero counts after filtering
ps_fung_before <- phyloseq::prune_taxa(taxa_sums(ps_fung_before) > 0, ps_fung_before)
ps_fung_after <- phyloseq::prune_taxa(taxa_sums(ps_fung_after) > 0, ps_fung_after)

# Create networks ---------------------------------------------------------

# For bacterial communities
set.seed(seedID)  # For reproducibility

message("Creating bacterial networks...")
bact_before_net <- create_network(ps_bact_before, 
                                  thresh = 0.05, 
                                  min_prevalence = 0.1,  # Present in 10% of samples
                                  min_abundance = 0.001, # At least 0.1% abundance
                                  max_taxa = 100)         # Top 40 taxa by abundance

bact_after_net <- create_network(ps_bact_after, 
                                 thresh = 0.05, 
                                 min_prevalence = 0.1,  # Present in 10% of samples
                                 min_abundance = 0.001, # At least 0.1% abundance
                                 max_taxa = 100)         # Top 40 taxa by abundance

# For fungal communities
message("Creating fungal networks...")
fung_before_net <- create_network(ps_fung_before,
                                  thresh = 0.05,
                                  min_prevalence = 0.1, 
                                  min_abundance = 0.001,
                                  max_taxa = 100)
fung_after_net <- create_network(ps_fung_after,
                                 thresh = 0.05,
                                 min_prevalence = 0.1, 
                                 min_abundance = 0.001,
                                 max_taxa = 100)

# Save network objects for further analysis ------------------------------

saveRDS(list(
    bacteria_before = bact_before_net,
    bacteria_after = bact_after_net,
    fungi_before = fung_before_net,
    fungi_after = fung_after_net
), "data/network_analysis_results.rds")

