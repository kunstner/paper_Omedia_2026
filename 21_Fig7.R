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
library(ANCOMBC)

# Data --------------------------------------------------------------------

seedID <- 138
source('_helper.R')

ps_bact <- readRDS('data/bacteria_phyloseq_processed.rds')
ps_fung <- readRDS('data/fungi_phyloseq_processed.rds')

# Filter out NA type samples
ps_bact <- phyloseq::subset_samples(ps_bact, !is.na(type))
ps_fung <- phyloseq::subset_samples(ps_fung, !is.na(type))

ps_bact_t1 <- phyloseq::subset_samples(ps_bact, time == "T1")
ps_bact_t2 <- phyloseq::subset_samples(ps_bact, time == "T2")

ps_fung_t1 <- phyloseq::subset_samples(ps_fung, time == "T1")
ps_fung_t2 <- phyloseq::subset_samples(ps_fung, time == "T2")

# Differential abundance analysis by timepoint ----------------------------

# Bacteria T1 (Before)
ancombc_result_bact_t1 <- ANCOMBC::ancombc2(
    data = ps_bact_t1,
    tax_level = "species",
    fix_formula = "type",
    rand_formula = NULL
)

res_df_bact_t1 <- ancombc_result_bact_t1$res %>% 
    dplyr::mutate(
        taxon_clean = stringr::str_remove(taxon, "^[^_]+_"),
        taxon = taxon_clean
    ) %>%
    dplyr::select(-taxon_clean)

p_da_bact_t1 <- res_df_bact_t1 %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon)) %>% 
    dplyr::filter(passed_ss_typeResponder == TRUE & p_typeResponder < 0.001) %>%
    dplyr::filter(abs(lfc_typeResponder) > 0.5) %>% 
    dplyr::arrange(desc(abs(lfc_typeResponder))) %>%
    dplyr::mutate(
        taxon = factor(taxon, levels = taxon[order(lfc_typeResponder)]),
        Direction = ifelse(lfc_typeResponder > 0, "Higher in Responder", "Higher in Non-Responder")
    ) %>%
    ggplot(aes(x = lfc_typeResponder, y = taxon, fill = Direction)) +
    geom_col() +
    geom_errorbar(
        aes(xmin = lfc_typeResponder - se_typeResponder, xmax = lfc_typeResponder + se_typeResponder),
        width = 0.2
    ) +
    scale_fill_manual(values = c("Higher in Responder" = "#2E86AB", "Higher in Non-Responder" = "#A23B72")) +
    labs(
        title = 'Bacteria T1 (Before): Differential abundance (p < 0.001)',
        x = "Log2 Fold Change (Responder/Non-Responder)",
        y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(angle = 0, size=10, face = 'italic')) +
    xlim(-7.5, 7.5)
p_da_bact_t1

# Bacteria T2 (After)
ancombc_result_bact_t2 <- ANCOMBC::ancombc2(
    data = ps_bact_t2,
    tax_level = "species",
    fix_formula = "type",
    rand_formula = NULL
)

res_df_bact_t2 <- ancombc_result_bact_t2$res %>% 
    dplyr::mutate(
        taxon_clean = stringr::str_remove(taxon, "^[^_]+_"),
        taxon = taxon_clean
    ) %>%
    dplyr::select(-taxon_clean)

p_da_bact_t2 <- res_df_bact_t2 %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon)) %>% 
    dplyr::filter(passed_ss_typeResponder == TRUE & p_typeResponder < 0.001) %>%
    dplyr::filter(abs(lfc_typeResponder) > 0.5) %>% 
    dplyr::arrange(desc(abs(lfc_typeResponder))) %>%
    dplyr::mutate(
        taxon = factor(taxon, levels = taxon[order(lfc_typeResponder)]),
        Direction = ifelse(lfc_typeResponder > 0, "Higher in Responder", "Higher in Non-Responder")
    ) %>%
    ggplot(aes(x = lfc_typeResponder, y = taxon, fill = Direction)) +
    geom_col() +
    geom_errorbar(
        aes(xmin = lfc_typeResponder - se_typeResponder, xmax = lfc_typeResponder + se_typeResponder),
        width = 0.2
    ) +
    scale_fill_manual(values = c("Higher in Responder" = "#2E86AB", "Higher in Non-Responder" = "#A23B72")) +
    labs(
        title = 'Bacteria T2 (After): Differential abundance (p < 0.001)',  # Fixed title
        x = "Log2 Fold Change (Responder/Non-Responder)",
        y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(angle = 0, size=10, face = 'italic')) +
    xlim(-7.5, 7.5)
p_da_bact_t2

p_da_bact_t1 + coord_flip() +
    labs(title = NULL, subtitle = NULL) +
    theme(
        axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "plain", angle = 0, hjust = 1)) + 
    p_da_bact_t2 + coord_flip() +
    labs(title = NULL, subtitle = NULL) +
    theme(
        axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "plain", angle = 0, hjust = 1)) +
    plot_layout(ncol = 2, widths = c(12,21), guides = "collect") + 
    plot_annotation(tag_levels = 'A') &
    theme(legend.position = "bottom")   # <- apply once
ggsave(filename = 'plots/Fig7.pdf', width = 15, height = 6)
# ggsave(filename = 'plots/Fig7.jpg', width = 15, height = 6, units = "in", dpi = 300)
# ggsave(filename = 'plots/Fig7.png', width = 15, height = 6, units = "in", dpi = 300)
