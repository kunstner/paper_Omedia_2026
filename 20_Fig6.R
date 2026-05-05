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
    as_tibble() %>% 
    dplyr::filter( !is.na(type)) 

shannon_with_metadata %>%
    dplyr::count(type, time) %>%
    tidyr::pivot_wider(names_from = time, values_from = n, values_fill = 0) %>%
    gt::gt() %>% 
    gt::tab_header(title = "Bacteria") 

shannon_with_metadata %>% 
    dplyr::group_by(type) %>%
    dplyr::summarise(
        mean_shannon = mean(dn_est),
        se_shannon = sd(dn_est)/sqrt(n()),
        lower_ci = mean_shannon - 1.96*se_shannon,
        upper_ci = mean_shannon + 1.96*se_shannon,
        n = n()
    )

shannon_before <- shannon_with_metadata %>% dplyr::filter(time == "T1")
shannon_after <- shannon_with_metadata %>% dplyr::filter(time == "T2")

# Test before/after timepoint
# t_test_before <- wilcox.test(dn_est ~ type, data = shannon_before)
# t_test_after <- wilcox.test(dn_est ~ type, data = shannon_after)
# lm_result <- lm(dn_est ~ type + Patient, data = shannon_with_metadata)
# summary(lm_result)

b_test_t1 <- breakaway::betta(
    chats = shannon_before$dn_est,
    ses = shannon_before$dn_error,
    X = model.matrix(~type, data = shannon_before),
    p.digits = 5)$table

b_test_result_t1 <- b_test_t1 %>% 
    data.frame() %>% dplyr::pull(p.values)

b_test_t2 <- breakaway::betta(
    chats = shannon_after$dn_est,
    ses = shannon_after$dn_error,
    X = model.matrix(~type, data = shannon_after),
    p.digits = 5)$table

b_test_result_t2 <- b_test_t2 %>% 
    data.frame() %>% dplyr::pull(p.values)

print("T1 (Before) betta results:")
print(b_test_t1)

print("T2 (After) betta results:")
print(b_test_t2)

# Visualization with uncertainty
p_a1 <- ggplot(shannon_with_metadata, aes(x = type, y = dn_est, fill = type)) +
    geom_violin(width = 0.75) +
    geom_boxplot(width = 0.2, fill = 'white') +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = colv_resp) +
    facet_wrap(~time) +
    # stat_compare_means() + 
    theme_def +
    theme(legend.position = 'none') +
    labs(
        title = "Bacteria: DivNet Shannon Diversity Estimates",
        subtitle = paste("before/after p-value =", 
                         format(b_test_result_t1[2], digits = 3), "/",
                         format(b_test_result_t2[2], digits = 3)),
        x = "",
        y = "DivNet estimate of Shannon"
    ) +
    ylim(0,5)
p_a1

# Fungi
shannon_with_metadata <- ps_fung %>% 
    sample_data() %>%
    as_tibble()

shannon_with_metadata %>%
    dplyr::count(type, time) %>%
    tidyr::pivot_wider(names_from = time, values_from = n, values_fill = 0) %>%
    gt::gt() %>% 
    gt::tab_header(title = "Fungi") 

shannon_with_metadata %>% 
    dplyr::group_by(type) %>%
    dplyr::summarise(
        mean_shannon = mean(dn_est),
        se_shannon = sd(dn_est)/sqrt(n()),
        lower_ci = mean_shannon - 1.96*se_shannon,
        upper_ci = mean_shannon + 1.96*se_shannon,
        n = n()
    )

shannon_before <- shannon_with_metadata %>% dplyr::filter(time == "T1")
shannon_after <- shannon_with_metadata %>% dplyr::filter(time == "T2")

# Test before/after timepoint
# t_test_before <- wilcox.test(dn_est ~ type, data = shannon_before)
# t_test_after <- wilcox.test(dn_est ~ type, data = shannon_after)
# lm_result <- lm(dn_est ~ type + Patient, data = shannon_with_metadata)
# summary(lm_result)

b_test_t1 <- breakaway::betta(
    chats = shannon_before$dn_est,
    ses = shannon_before$dn_error,
    X = model.matrix(~type, data = shannon_before),
    p.digits = 5)$table

b_test_result_t1 <- b_test_t1 %>% 
    data.frame() %>% dplyr::pull(p.values)

b_test_t2 <- breakaway::betta(
    chats = shannon_after$dn_est,
    ses = shannon_after$dn_error,
    X = model.matrix(~type, data = shannon_after),
    p.digits = 5)$table

b_test_result_t2 <- b_test_t2 %>% 
    data.frame() %>% dplyr::pull(p.values)

print("T1 (Before) betta results:")
print(b_test_t1)

print("T2 (After) betta results:")
print(b_test_t2)

# Visualization with uncertainty
p_a2 <- ggplot(shannon_with_metadata, aes(x = type, y = dn_est, fill = type)) +
    geom_violin(width = 0.75) +
    geom_boxplot(width = 0.2, fill = 'white') +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~time) +
    scale_fill_manual(values = colv_resp) +
    # stat_compare_means()+
    theme_def +
    theme(legend.position = 'none') +
    labs(
        title = "Fungi: DivNet Shannon Diversity Estimates",
        subtitle = paste("before/after p-value =", 
                         format(b_test_result_t1[2], digits = 3), "/",
                         format(b_test_result_t2[2], digits = 3)),      
        x = "",
        y = "DivNet estimate of Shannon"
    ) +
    ylim(0,5)
p_a2

p_a1 + theme(legend.position = "none") + 
    p_a2 + theme(legend.position = "none") +
    plot_annotation(tag_levels = 'A')

# Beta diversity ----------------------------------------------------------

# Bacteria
ps_clr <- microbiome::transform(x = ps_bact, transform = "clr")

# Filter by timepoint and remove samples with NA type
ps_clr_t1 <- phyloseq::subset_samples(ps_clr, time == "T1" & !is.na(type))
ps_clr_t2 <- phyloseq::subset_samples(ps_clr, time == "T2" & !is.na(type))

# T1 (Before) analysis
otu_table_clr_t1 <- otu_table(ps_clr_t1) %>% t()
ps_clr_dist_t1 <- dist(otu_table_clr_t1, method="euclidean")

# PERMANOVA for T1
adonis_result_t1 <- vegan::adonis2(
    formula = ps_clr_dist_t1 ~ type,
    data = data.frame(ps_clr_t1 %>% sample_data),
    permutations = 9999
)
print("T1 (Before) PERMANOVA:")
print(adonis_result_t1)

# T2 (After) analysis  
otu_table_clr_t2 <- otu_table(ps_clr_t2) %>% t()
ps_clr_dist_t2 <- dist(otu_table_clr_t2, method="euclidean")

# PERMANOVA for T2
adonis_result_t2 <- vegan::adonis2(
    formula = ps_clr_dist_t2 ~ type,
    data = data.frame(ps_clr_t2 %>% sample_data),
    permutations = 9999
)
print("T2 (After) PERMANOVA:")
print(adonis_result_t2)

# Ordination plots
ps_clr_ord_t1 <- phyloseq::ordinate(ps_clr_t1, "RDA", distance = "euclidean")
p_beta_1_t1 <- plot_ordination(
    physeq = ps_clr_t1,
    ordination = ps_clr_ord_t1, color='type')  +
    scale_color_manual(values = colv_resp) +
    geom_point(size=3) +
    stat_ellipse(type = "t", linetype = 2, level = 0.95) +
    ggtitle("Bacteria T1 (Before): RDA of Aitchison distance") +
    theme_def +
    theme(legend.position = 'bottom') +
    labs(subtitle = paste("PERMANOVA p-value =", format(adonis_result_t1$`Pr(>F)`[1], digits = 3)))

ps_clr_ord_t2 <- phyloseq::ordinate(ps_clr_t2, "RDA", distance = "euclidean")
p_beta_1_t2 <- plot_ordination(
    physeq = ps_clr_t2,
    ordination = ps_clr_ord_t2, color='type')  +
    scale_color_manual(values = colv_resp) +
    geom_point(size=3) +
    stat_ellipse(type = "t", linetype = 2, level = 0.95) +
    ggtitle("Bacteria T2 (After): RDA of Aitchison distance") +
    theme_def +
    theme(legend.position = 'none') +
    labs(subtitle = paste("PERMANOVA p-value =", format(adonis_result_t2$`Pr(>F)`[1], digits = 3)))

# Fungi
ps_clr <- microbiome::transform(x = ps_fung, transform = "clr")

# Filter by timepoint and remove samples with NA type
ps_clr_t1 <- phyloseq::subset_samples(ps_clr, time == "T1" & !is.na(type))
ps_clr_t2 <- phyloseq::subset_samples(ps_clr, time == "T2" & !is.na(type))

# T1 (Before) analysis
otu_table_clr_t1 <- otu_table(ps_clr_t1) %>% t()
ps_clr_dist_t1 <- dist(otu_table_clr_t1, method="euclidean")

# PERMANOVA for T1
adonis_result_t1 <- vegan::adonis2(
    formula = ps_clr_dist_t1 ~ type,
    data = data.frame(ps_clr_t1 %>% sample_data),
    permutations = 9999
)
print("Fungi T1 (Before) PERMANOVA:")
print(adonis_result_t1)

# T2 (After) analysis
otu_table_clr_t2 <- otu_table(ps_clr_t2) %>% t()
ps_clr_dist_t2 <- dist(otu_table_clr_t2, method="euclidean")

# PERMANOVA for T2
adonis_result_t2 <- vegan::adonis2(
    formula = ps_clr_dist_t2 ~ type,
    data = data.frame(ps_clr_t2 %>% sample_data),
    permutations = 9999
)
print("Fungi T2 (After) PERMANOVA:")
print(adonis_result_t2)

# Ordination plots
ps_clr_ord_t1 <- phyloseq::ordinate(ps_clr_t1, "RDA", distance = "euclidean")
p_beta_2_t1 <- plot_ordination(
    physeq = ps_clr_t1,
    ordination = ps_clr_ord_t1, color='type')  +
    scale_color_manual(values = colv_resp) +
    geom_point(size=3) +
    stat_ellipse(type = "t", linetype = 2, level = 0.95) +
    ggtitle("Fungi T1 (Before): RDA of Aitchison distance") +
    theme_def +
    theme(legend.position = 'bottom') +
    labs(subtitle = paste("PERMANOVA p-value =", format(adonis_result_t1$`Pr(>F)`[1], digits = 3)))

ps_clr_ord_t2 <- phyloseq::ordinate(ps_clr_t2, "RDA", distance = "euclidean")
p_beta_2_t2 <- plot_ordination(
    physeq = ps_clr_t2,
    ordination = ps_clr_ord_t2, color='type')  +
    scale_color_manual(values = colv_resp) +
    geom_point(size=3) +
    stat_ellipse(type = "t", linetype = 2, level = 0.95) +
    ggtitle("Fungi T2 (After): RDA of Aitchison distance") +
    theme_def +
    theme(legend.position = 'none') +
    labs(subtitle = paste("PERMANOVA p-value =", format(adonis_result_t2$`Pr(>F)`[1], digits = 3)))

p_a1 +
    labs(title = NULL, subtitle = NULL) + 
    p_beta_1_t1 + 
    labs(title = NULL, subtitle = NULL) + 
    theme(legend.position = "none") + 
    p_beta_1_t2 + 
    labs(title = NULL, subtitle = NULL) + 
    theme(legend.position = "none") +
    plot_annotation(tag_levels = 'A') +
    plot_layout(ncol = 3, widths = c(2,3,3))
ggsave(filename = "plots/Fig6.pdf", height = 6, width = 15)
