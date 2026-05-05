# Credentials -------------------------------------------------------------

#
# Author: Axel Künstner
# Project: Otitis media microbiome
# Data: Microbiome

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(igraph)
library(ggnetwork)
library(ggrepel)
library(patchwork)

# Data --------------------------------------------------------------------

source('_helper.R')
seedID <- 138

results <- readRDS("data/network_analysis_results.rds")

bact_before_net <- results$bacteria_before
bact_after_net  <- results$bacteria_after
fung_before_net <- results$fungi_before
fung_after_net  <- results$fungi_after

# Fig 4: Genus-level interaction networks ---------------------------------

create_genus_interactions <- function(network_obj, title = "",
                                      label_size = 4,
                                      label_top_n = 15,
                                      label_by = "degree") {
    g <- network_obj$graph
    
    taxa_info <- tibble(
        taxon  = igraph::V(g)$name,
        genus  = igraph::vertex_attr(g, "genus"),
        phylum = igraph::vertex_attr(g, "phylum")
    )
    
    edges_df <- igraph::as_data_frame(g, what = "edges") |>
        dplyr::left_join(taxa_info |> dplyr::select(taxon, genus),
                         by = c("from" = "taxon")) |>
        dplyr::rename(from_genus = genus) |>
        dplyr::left_join(taxa_info |> dplyr::select(taxon, genus),
                         by = c("to" = "taxon")) |>
        dplyr::rename(to_genus = genus)
    
    genus_interactions <- edges_df |>
        dplyr::filter(!is.na(from_genus) & !is.na(to_genus)) |>
        dplyr::group_by(from_genus, to_genus) |>
        dplyr::summarise(
            weight     = sum(abs(weight)),
            count      = dplyr::n(),
            avg_weight = mean(abs(weight)),
            .groups    = "drop"
        ) |>
        dplyr::filter(count >= 2 | avg_weight >= 0.1) |>
        dplyr::arrange(dplyr::desc(weight)) |>
        dplyr::filter(from_genus != to_genus)
    
    genus_nodes  <- unique(c(genus_interactions$from_genus, genus_interactions$to_genus))
    genus_graph  <- igraph::graph_from_data_frame(
        d        = genus_interactions,
        vertices = genus_nodes,
        directed = FALSE
    ) |>
        igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE,
                         edge.attr.comb = "sum")
    
    igraph::V(genus_graph)$degree      <- igraph::degree(genus_graph)
    igraph::V(genus_graph)$betweenness <- igraph::betweenness(genus_graph, normalized = TRUE)
    igraph::V(genus_graph)$genus_name  <- igraph::V(genus_graph)$name
    
    genus_phylum <- taxa_info |>
        dplyr::filter(!is.na(genus)) |>
        dplyr::select(genus, phylum) |>
        dplyr::distinct()
    
    for (i in seq_len(nrow(genus_phylum))) {
        genus  <- genus_phylum$genus[i]
        phylum <- genus_phylum$phylum[i]
        if (genus %in% igraph::V(genus_graph)$name) {
            igraph::V(genus_graph)[genus]$phylum <- phylum
        }
    }
    
    set.seed(42)
    layout       <- igraph::layout_with_fr(genus_graph)
    genus_net_df <- ggnetwork::ggnetwork(genus_graph, layout = layout)
    nodes_df     <- genus_net_df |> dplyr::distinct(name, .keep_all = TRUE)
    
    message("Number of unique nodes found: ", nrow(nodes_df))
    
    if (label_by == "degree") {
        label_threshold <- sort(nodes_df$degree, decreasing = TRUE)[min(label_top_n, nrow(nodes_df))]
        label_df <- nodes_df |> dplyr::filter(degree >= label_threshold)
    } else if (label_by == "betweenness") {
        label_threshold <- sort(nodes_df$betweenness, decreasing = TRUE)[min(label_top_n, nrow(nodes_df))]
        label_df <- nodes_df |> dplyr::filter(betweenness >= label_threshold)
    } else {
        label_df <- nodes_df
    }
    
    message("Number of labels to show: ", nrow(label_df))
    
    p <- ggplot(genus_net_df, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(alpha = weight, linewidth = weight),
                   color = "gray40", curvature = 0.1) +
        geom_nodes(aes(color = phylum, size = degree)) +
        geom_text(data = label_df,
                  aes(label = genus_name),
                  size = label_size, fontface = "italic",
                  color = "black", vjust = -0.8, hjust = 0.5) +
        scale_size_continuous(range = c(0.5, 3)) +
        scale_alpha(range = c(0.2, 0.8)) +
        scale_color_brewer(palette = "Set1") +
        scale_x_continuous(expand = expansion(mult = 0.15)) +
        scale_y_continuous(expand = expansion(mult = 0.15)) +
        theme_void() +
        theme(
            legend.position  = "right",
            legend.text      = element_text(size = 10),
            legend.title     = element_text(size = 11, face = "bold"),
            plot.title       = element_text(size = 12),
            plot.margin      = margin(10, 10, 10, 10)
        ) +
        labs(
            title    = title,
            subtitle = paste0("Genera: ", igraph::vcount(genus_graph),
                              " | Connections: ", igraph::ecount(genus_graph)),
            color     = "Phylum",
            linewidth = "Interaction\nStrength"
        ) +
        guides(alpha = "none")
    
    return(p)
}

p_bact_before_genus <- create_genus_interactions(bact_before_net,
                                                 "Bacterial Genus Interactions (before treatment)",
                                                 label_size = 5, label_top_n = 100, label_by = "degree")

p_bact_after_genus <- create_genus_interactions(bact_after_net,
                                                "Bacterial Genus Interactions (after treatment)",
                                                label_size = 5, label_top_n = 100, label_by = "degree")

p_fung_before_genus <- create_genus_interactions(fung_before_net,
                                                 "Fungal Genus Interactions (before treatment)",
                                                 label_size = 5, label_top_n = 100, label_by = "degree")

p_fung_after_genus <- create_genus_interactions(fung_after_net,
                                                "Fungal Genus Interactions (after treatment)",
                                                label_size = 5, label_top_n = 100, label_by = "degree")

(p_bact_before_genus + p_bact_after_genus) /
    (p_fung_before_genus + p_fung_after_genus) +
    plot_annotation(tag_levels = 'A')
ggsave("plots/Fig4.pdf", width = 16, height = 14)

# Fig 5: Cross-domain interactions ----------------------------------------

analyze_cross_domain_interactions <- function(ps_bact, ps_fung,
                                              cor_method = "spearman",
                                              min_cor = 0.3,
                                              max_p = 0.1,
                                              max_taxa = 100) {
    common_samples   <- intersect(sample_names(ps_bact), sample_names(ps_fung))
    ps_bact_filtered <- phyloseq::prune_samples(common_samples, ps_bact)
    ps_fung_filtered <- phyloseq::prune_samples(common_samples, ps_fung)
    
    bact_prev  <- phyloseq::taxa_sums(ps_bact_filtered) / phyloseq::nsamples(ps_bact_filtered)
    bact_taxa  <- names(sort(bact_prev, decreasing = TRUE)[1:min(max_taxa, length(bact_prev))])
    ps_bact_filtered <- phyloseq::prune_taxa(bact_taxa, ps_bact_filtered)
    
    fung_prev  <- phyloseq::taxa_sums(ps_fung_filtered) / phyloseq::nsamples(ps_fung_filtered)
    fung_taxa  <- names(sort(fung_prev, decreasing = TRUE)[1:min(max_taxa, length(fung_prev))])
    ps_fung_filtered <- phyloseq::prune_taxa(fung_taxa, ps_fung_filtered)
    
    bact_abund <- otu_table(ps_bact_filtered)
    if (taxa_are_rows(ps_bact_filtered)) bact_abund <- t(bact_abund)
    fung_abund <- otu_table(ps_fung_filtered)
    if (taxa_are_rows(ps_fung_filtered)) fung_abund <- t(fung_abund)
    
    bact_df        <- as.data.frame(bact_abund)
    fung_df        <- as.data.frame(fung_abund)
    bact_taxa_info <- as.data.frame(tax_table(ps_bact_filtered))
    fung_taxa_info <- as.data.frame(tax_table(ps_fung_filtered))
    
    message("Calculating correlations for ", ncol(bact_df), " bacteria and ",
            ncol(fung_df), " fungi...")
    
    cors_matrix <- cor(bact_df, fung_df, method = cor_method,
                       use = "pairwise.complete.obs")
    sig_cors    <- which(abs(cors_matrix) >= min_cor, arr.ind = TRUE)
    
    if (nrow(sig_cors) == 0) {
        message("No significant correlations found.")
        return(tibble::tibble(
            bact_taxon = character(), fung_taxon = character(),
            correlation = numeric(), p_value = numeric(),
            p_adjusted = numeric(), bact_phylum = character(),
            bact_genus = character(), fung_phylum = character(),
            fung_genus = character()
        ))
    }
    
    message("Testing ", nrow(sig_cors), " potentially significant correlations...")
    
    cors_df <- purrr::map_dfr(seq_len(nrow(sig_cors)), function(i) {
        row_idx    <- sig_cors[i, 1]
        col_idx    <- sig_cors[i, 2]
        bact_name  <- rownames(cors_matrix)[row_idx]
        fung_name  <- colnames(cors_matrix)[col_idx]
        test_result <- cor.test(bact_df[[bact_name]], fung_df[[fung_name]],
                                method = cor_method, exact = FALSE)
        tibble::tibble(
            bact_taxon  = bact_name,
            fung_taxon  = fung_name,
            correlation = cors_matrix[row_idx, col_idx],
            p_value     = test_result$p.value
        )
    }) |>
        dplyr::mutate(p_adjusted = p.adjust(p_value, method = "BH")) |>
        dplyr::filter(abs(correlation) >= min_cor & p_adjusted <= max_p) |>
        dplyr::arrange(dplyr::desc(abs(correlation)))
    
    if (nrow(cors_df) == 0) {
        message("No significant correlations after FDR correction.")
        return(cors_df)
    }
    
    cors_df |>
        dplyr::left_join(
            bact_taxa_info |>
                dplyr::mutate(taxon = rownames(bact_taxa_info)) |>
                dplyr::select(taxon, bact_phylum = phylum,
                              bact_genus = genus, bact_species = species),
            by = c("bact_taxon" = "taxon")
        ) |>
        dplyr::left_join(
            fung_taxa_info |>
                dplyr::mutate(taxon = rownames(fung_taxa_info)) |>
                dplyr::select(taxon, fung_phylum = phylum,
                              fung_genus = genus, fung_species = species),
            by = c("fung_taxon" = "taxon")
        )
}

plot_cross_domain_interactions <- function(interactions_df, title = "",
                                           labelsize = 2) {
    plot_data <- interactions_df |>
        dplyr::mutate(
            bact_label       = ifelse(!is.na(bact_genus), bact_genus, bact_taxon),
            fung_label       = ifelse(!is.na(fung_genus), fung_genus, fung_taxon),
            interaction_type = ifelse(correlation > 0, "Positive", "Negative")
        ) |>
        dplyr::slice_max(order_by = abs(correlation), n = 30)
    
    edges <- tibble::tibble(
        from             = plot_data$bact_label,
        to               = plot_data$fung_label,
        weight           = plot_data$correlation,
        interaction_type = plot_data$interaction_type
    )
    
    nodes <- tibble::tibble(
        name   = c(unique(edges$from), unique(edges$to)),
        domain = c(rep("Bacteria", length(unique(edges$from))),
                   rep("Fungi",    length(unique(edges$to))))
    )
    
    cross_graph <- igraph::graph_from_data_frame(
        d = edges, vertices = nodes, directed = FALSE)
    igraph::V(cross_graph)$degree <- igraph::degree(cross_graph)
    
    set.seed(42)
    layout       <- igraph::layout_with_fr(cross_graph)
    cross_net_df <- ggnetwork::ggnetwork(cross_graph, layout = layout)
    label_df     <- cross_net_df |> dplyr::distinct(name, .keep_all = TRUE)
    
    domain_colors      <- c("Bacteria" = "#0072B2", "Fungi" = "#D55E00")
    interaction_colors <- c("Positive" = "#00BA38", "Negative" = "#F8766D")
    
    ggplot(cross_net_df, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(color = interaction_type,
                       alpha = abs(weight),
                       linewidth = abs(weight)),
                   curvature = 0.1) +
        geom_nodes(aes(color = domain, size = degree)) +
        ggrepel::geom_label_repel(
            data             = label_df,
            aes(x = x, y = y, label = name, fill = domain),
            size             = labelsize,
            alpha            = 0.85,
            label.size       = 0.2,
            color            = "white",
            label.padding    = unit(0.15, "lines"),
            box.padding      = 0.5,
            max.overlaps     = 20,
            segment.alpha    = 0.75,
            segment.colour   = "black",
            segment.size     = 0.3,
            min.segment.length = 0,
            inherit.aes      = FALSE
        ) +
        scale_color_manual(values = c(domain_colors, interaction_colors)) +
        scale_fill_manual(values = domain_colors) +
        scale_linewidth_continuous(range = c(0.5, 3)) +
        scale_size_continuous(range = c(0.5, 3)) +
        scale_alpha(range = c(0.2, 0.8)) +
        theme_void() +
        theme(legend.position = "right",
              plot.title = element_text(size = 12)) +
        guides(
            color     = guide_legend(title = NULL, order = 1),
            fill      = guide_legend(title = "Domain", order = 2),
            size      = guide_legend(title = "Interaction\nStrength", order = 3),
            alpha     = "none",
            linewidth = "none"
        ) +
        labs(title    = title,
             subtitle = "Significant correlations between bacterial and fungal taxa")
}

ps_bact <- readRDS('data/bacteria_phyloseq_processed.rds')
ps_fung <- readRDS('data/fungi_phyloseq_processed.rds')

ps_bact_before <- phyloseq::subset_samples(ps_bact, soledum == "before")
ps_bact_after  <- phyloseq::subset_samples(ps_bact, soledum == "after")
ps_fung_before <- phyloseq::subset_samples(ps_fung, soledum == "before")
ps_fung_after  <- phyloseq::subset_samples(ps_fung, soledum == "after")

cross_interactions_before <- analyze_cross_domain_interactions(
    ps_bact_before, ps_fung_before,
    min_cor = 0.4, max_p = 0.05, max_taxa = 250)

cross_interactions_after <- analyze_cross_domain_interactions(
    ps_bact_after, ps_fung_after,
    min_cor = 0.4, max_p = 0.05, max_taxa = 250)

p_cross_before <- plot_cross_domain_interactions(
    cross_interactions_before,
    "Cross-domain Interactions Before Treatment", labelsize = 4)

p_cross_after <- plot_cross_domain_interactions(
    cross_interactions_after,
    "Cross-domain Interactions After Treatment", labelsize = 4)

p_cross_before + p_cross_after +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = 'A')
ggsave("plots/Fig5.pdf",
       width = 18, height = 10)
