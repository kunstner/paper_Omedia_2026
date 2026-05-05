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

# Compare alpha diversity -------------------------------------------------

# Bacteria
shannon_with_metadata <- ps_bact %>% 
    sample_data() %>%
    as_tibble()
shannon_with_metadata %>% 
    dplyr::group_by(soledum) %>%
    dplyr::summarise(
        mean_shannon = mean(dn_est),
        se_shannon = sd(dn_est)/sqrt(n()),
        lower_ci = mean_shannon - 1.96*se_shannon,
        upper_ci = mean_shannon + 1.96*se_shannon,
        n = n()
    )

# t_test_result <- t.test(dn_est ~ soledum, data = shannon_with_metadata)
# w_test_result <- wilcox.test(dn_est ~ soledum, data = shannon_with_metadata)
# lm_result <- lm(dn_est ~ soledum + Patient, data = shannon_with_metadata)
# summary(lm_result)

b_test_table <- breakaway::betta(chats = shannon_with_metadata$dn_est,
                 ses = shannon_with_metadata$dn_error,
                 X = model.matrix(~soledum, data = shannon_with_metadata),
                 p.digits = 5)$table
b_test_result <- b_test_table %>% 
    data.frame() %>% pull(p.values)
b_test_table

# Visualization with uncertainty
p_a1 <- ggplot(shannon_with_metadata, aes(x = soledum, y = dn_est, fill = soledum)) +
    geom_violin(width = 0.75) +
    geom_boxplot(width = 0.2, fill = 'white') +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = colv2) +
    # stat_compare_means() + 
    theme_def +
    labs(
        title = "Bacteria: DivNet Shannon Diversity Estimates",
        subtitle = paste("p-value =", format(b_test_result[2], digits = 3)),
        x = "Treatment Status",
        y = "DivNet estimate of Shannon"
    ) +
    ylim(0,5)
p_a1

# Fungi
shannon_with_metadata <- ps_fung %>% 
    sample_data() %>%
    as_tibble()
shannon_with_metadata %>% 
    dplyr::group_by(soledum) %>%
    dplyr::summarise(
        mean_shannon = mean(dn_est),
        se_shannon = sd(dn_est)/sqrt(n()),
        lower_ci = mean_shannon - 1.96*se_shannon,
        upper_ci = mean_shannon + 1.96*se_shannon,
        n = n()
    )

# t_test_result <- t.test(dn_est ~ soledum, data = shannon_with_metadata)
# w_test_result <- wilcox.test(dn_est ~ soledum, data = shannon_with_metadata)
# lm_result <- lm(dn_est ~ soledum + Patient, data = shannon_with_metadata)
# summary(lm_result)

b_test_table <- breakaway::betta(chats = shannon_with_metadata$dn_est,
                                 ses = shannon_with_metadata$dn_error,
                                 X = model.matrix(~soledum, data = shannon_with_metadata),
                                 p.digits = 5)$table
b_test_result <- b_test_table %>% 
    data.frame() %>% pull(p.values)
b_test_table

# Visualization with uncertainty
p_a2 <- ggplot(shannon_with_metadata, aes(x = soledum, y = dn_est, fill = soledum)) +
    geom_violin(width = 0.75) +
    geom_boxplot(width = 0.2, fill = 'white') +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = colv2) +
    # stat_compare_means()+ 
    theme_def +
    labs(
        title = "Fungi: DivNet Shannon Diversity Estimates",
        subtitle = paste("p-value =", format(b_test_result[2], digits = 3)),
        x = "Treatment Status",
        y = "DivNet estimate of Shannon"
    ) +
    ylim(0,5)
p_a2

p_a1 + theme(legend.position = "bottom") + 
    p_a2 + theme(legend.position = "none") +
    plot_annotation(tag_levels = 'A')

# Beta diversity ----------------------------------------------------------

# Bacteria
ps_clr <- microbiome::transform(x = ps_bact, transform = "clr")
otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")

# PERMANOVA    
adonis_result <- vegan::adonis2(
    formula = ps_clr_dist ~ soledum,
    data = data.frame(ps_clr %>% sample_data),
    permutations = 9999
)
adonis_result

ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
p_beta_1 <- plot_ordination(
    physeq = ps_clr,
    ordination = ps_clr_ord, color='soledum')  +
    scale_color_manual(values = colv2) +
    geom_point(size=3) +
    # geom_text(mapping = aes(label = Group), size = 3, nudge_x = 0.25) +
    stat_ellipse(type = "t", linetype = 2, level = 0.95) +
    ggtitle("Bacteria: RDA of Aitchison distance") +
    theme_def +
    theme(legend.position = 'bottom') +
    labs(subtitle = paste("PERMANOVA p-value =", format(adonis_result$`Pr(>F)`[1], digits = 3)))
p_beta_1

# Fungi
ps_clr <- microbiome::transform(x = ps_fung, transform = "clr")
otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")

# PERMANOVA    
adonis_result <- vegan::adonis2(
    formula = ps_clr_dist ~ soledum,
    data = data.frame(ps_clr %>% sample_data),
    permutations = 9999
)
adonis_result

ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
p_beta_2 <- plot_ordination(
    physeq = ps_clr,
    ordination = ps_clr_ord, color='soledum')  +
    scale_color_manual(values = colv2) +
    geom_point(size=3) +
    # geom_text(mapping = aes(label = Group), size = 3, nudge_x = 0.25) +
    stat_ellipse(type = "t", linetype = 2, level = 0.95) +
    ggtitle("Fungi: RDA of Aitchison distance") +
    theme_def +
    theme(legend.position = 'bottom') +
    labs(subtitle = paste("PERMANOVA p-value =", format(adonis_result$`Pr(>F)`[1], digits = 3)))
p_beta_2

p_beta_1 + theme(legend.position = "bottom") + 
    p_beta_2 + theme(legend.position = "none") +
    plot_annotation(tag_levels = 'A')

p_a1 +
    labs(title = NULL, subtitle = NULL) +
    theme(legend.position = 'none') +
    p_a2 +
    labs(title = NULL, subtitle = NULL) +
    theme(legend.position = 'none') +
    p_beta_1 +
    labs(title = NULL, subtitle = NULL) +
    theme(legend.position = 'none') +
    p_beta_2 +
    labs(title = NULL, subtitle = NULL) +
    theme(legend.position = 'none') +
    plot_layout(ncol = 4, widths = c(2,2,3,3)) +
    plot_annotation(tag_levels = 'A')
ggsave(filename = "plots/Fig1.pdf", height = 6, width = 15)
