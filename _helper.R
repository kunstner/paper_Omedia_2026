
colv2 <- c("before" = "#0072B2", "after" = "#D55E00")
# colv_resp <- c("Responder" = "#5470c6", "Non-Responder" = "#ee6666") 
colv_resp <- c("Responder" = "#2E86AB", "Non-Responder" = "#A23B72")

theme_def <- theme(axis.line = element_line(colour = "black"),
                   legend.key = element_rect(fill = "white"),
                   legend.position = "right",
                   legend.title = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.x = element_text(angle = 0, size=12),
                   strip.text.x = element_text(angle = 0, size = 12),
                   axis.text.y = element_text(size=12)
                   # axis.ticks.x=element_blank()
)


# Functions ---------------------------------------------------------------

# Create a taxa plot for a given phyloseq object
create_taxa_plot <- function(ps, rank = "Phylum", group_var = NULL, 
                             top_n = NULL, other_label = "Other") {
    
    if(!is.null(top_n) && top_n > 20) {
        print('Max 20 taxa allowed, set top_n to 20')
        top_n <- 20
    }
    
    # Ensure proper column names
    colnames(ps@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Agglomerate at the specified taxonomic rank
    ps_glom <- phyloseq::tax_glom(ps, taxrank = rank)
    
    # Transform to relative abundance
    ps_rel <- phyloseq::transform_sample_counts(ps_glom, function(x) {x/sum(x)})
    
    # Melt to long format
    ps_melt <- phyloseq::psmelt(ps_rel) %>%
        dplyr::mutate(!!as.name(rank) := gsub("_", " ", !!as.name(rank)))
    
    # Ensure group_var is in the melted data
    if(!(group_var %in% colnames(ps_melt))) {
        stop(paste("Variable", group_var, "not found in sample data"))
    }
    
    # For lower taxonomic ranks, identify top taxa
    if (!is.null(top_n)) {
        # Calculate mean abundance for each taxon across all samples
        taxa_sums <- ps_melt %>%
            dplyr::group_by(!!as.name(rank)) %>%
            dplyr::summarize(mean_abundance = mean(Abundance), .groups = "drop") %>%
            dplyr::arrange(dplyr::desc(mean_abundance))
        
        # Get top taxa names
        top_taxa <- taxa_sums %>%
            dplyr::slice_head(n = top_n) %>%
            dplyr::pull(!!as.name(rank))
        
        # Replace taxa not in top_n with "Other"
        ps_melt <- ps_melt %>%
            dplyr::mutate(TaxaGroup = ifelse(!!as.name(rank) %in% top_taxa, 
                                             !!as.name(rank), 
                                             other_label))
        
        # Calculate sum by group per sample to preserve the stacking behavior
        plot_data <- ps_melt %>%
            dplyr::group_by(Sample, TaxaGroup, !!as.name(group_var)) %>%
            dplyr::summarize(Abundance = sum(Abundance), .groups = "drop") %>%
            # Then calculate the mean across samples for each group
            dplyr::group_by(TaxaGroup, !!as.name(group_var)) %>%
            dplyr::summarize(mean = mean(Abundance), .groups = "drop") %>%
            dplyr::rename(!!as.name(rank) := TaxaGroup)
    } else {
        # Just use all taxa with their individual abundances
        plot_data <- ps_melt %>%
            dplyr::group_by(!!as.name(rank), !!as.name(group_var)) %>%
            dplyr::summarize(mean = mean(Abundance), .groups = "drop")
    }
    
    # Order taxa by abundance for the stacked bars (with "Other" at the end)
    if (other_label %in% plot_data[[rank]]) {
        # Calculate mean abundance across groups for ordering
        taxa_order <- plot_data %>%
            dplyr::filter(!!as.name(rank) != other_label) %>% 
            dplyr::group_by(!!as.name(rank)) %>%
            dplyr::summarize(mean_abund = mean(mean), .groups = "drop") %>%
            dplyr::arrange(dplyr::desc(mean_abund)) %>%
            dplyr::pull(!!as.name(rank))
        
        # Add "Other" at the end
        taxa_order <- c(taxa_order, other_label)
    } else {
        taxa_order <- plot_data %>%
            dplyr::group_by(!!as.name(rank)) %>%
            dplyr::summarize(mean_abund = mean(mean), .groups = "drop") %>%
            dplyr::arrange(dplyr::desc(mean_abund)) %>%
            dplyr::pull(!!as.name(rank))
    }
    
    # Apply the factor ordering for the stacked bars
    plot_data[[rank]] <- factor(plot_data[[rank]], levels = taxa_order)
    
    # Generate colors
    n_taxa <- length(unique(plot_data[[rank]]))
    
    if (other_label %in% plot_data[[rank]]) {
        # If we have "Other", use n-1 colors from palette plus purple for "Other"
        top_n_colors <- pals::tableau20(min(20, n_taxa-1))
        # Use a distinctive purple for "Other"
        all_cols <- c(top_n_colors, "purple4")
    } else {
        # Otherwise just use palette colors
        all_cols <- pals::tableau20(min(20, n_taxa))
    }
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = !!as.name(group_var), y = mean, fill = !!as.name(rank))) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = all_cols) +
        theme_minimal() +
        theme(
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 0, size=12),
            strip.text.x = element_text(angle = 0, size = 12),
            axis.text.y = element_text(size=12)
        ) +
        ylab("Relative Abundance\n")
    
    # Apply italic formatting to legend if applicable
    if (rank %in% c("Family", "Genus", "Species")) {
        p <- p + guides(fill = guide_legend(title = rank, ncol = 1, 
                                            label.theme = element_text(face = "italic")))
    } else {
        p <- p + guides(fill = guide_legend(title = rank, ncol = 1))
    }
    
    return(p)
}


