library(pheatmap)
library(tidyverse)
library(readxl)
library(stringr)
library(knitr)

result_list <- readRDS("specific_snp_enrichment_pfdr_11_28_18.RData")
results_df <- bind_rows(result_list, .id = "lead_variant")
annotated_tissue <- read_excel("roadmap_instatehub_wilson.xlsx")
annotated_tissue <- annotated_tissue %>%
  filter(is.na(include)) %>%
  mutate_all(toupper) %>%
  dplyr::select(ab, annotation)

all_samples <- tbl_df(unique(results_df$sample))

results_df <- results_df %>%
  mutate(a1_fg = 0.5 + fg.success,
         b1_fg = 0.5 + fg.total - fg.success,
         a1_bg = 0.5 + bg.success,
         b1_bg = 0.5 + bg.total - bg.success,
         avg_success_fg = fg.success/fg.total,
         avg_success_bg = bg.success/bg.total)

# active enhancer tissue samples
ae_samples <- results_df %>%
  ungroup() %>%
  filter(state == "Active Enhancer") %>%
  dplyr::select(sample) %>%
  distinct()

# active promoter tissue samples
ap_samples <- results_df %>%
  ungroup() %>%
  filter(state == "Active Promoter") %>%
  dplyr::select(sample) %>%
  distinct()

table(all_samples$value %in% annotated_tissue$ab)

results_df <- results_df %>%
  left_join(., annotated_tissue, by = c("sample" = "ab"))

promoters <- results_df %>%
  ungroup() %>%
  filter(state == "Active Promoter") %>%
  dplyr::select(lead_variant, sample, fdr_sig) %>%
  spread(sample, fdr_sig) %>%
  base::as.data.frame()

rownames(promoters) <- promoters$lead_variant
promoters[promoters == TRUE] <- 1
promoters[promoters == FALSE]  <- 0
promoters <- promoters %>%
  dplyr::select(-lead_variant) %>%
  as.matrix()

col_groups <- results_df %>%
  ungroup() %>%
  distinct(sample, annotation) %>%
  arrange(sample) %>%
  dplyr::select(annotation) %>%
  base::as.data.frame()

row.names(col_groups) <- colnames(promoters)

pheatmap(promoters, cluster_cols = FALSE, main = "Active Promoter", annotation_col = col_groups, legend = FALSE,
         annotation_names_col = FALSE) 


# active enhancer
# active promoter tissue samples
ae_samples <- results_df %>%
  ungroup() %>%
  filter(state == "Active Enhancer") %>%
  dplyr::select(sample) %>%
  distinct()

enhancers <- results_df %>%
  ungroup() %>%
  filter(state == "Active Enhancer") %>%
  dplyr::select(lead_variant, sample, fdr_sig) %>%
  spread(sample, fdr_sig) %>%
  base::as.data.frame()

rownames(enhancers) <- enhancers$lead_variant
enhancers[enhancers == TRUE] <- 1
enhancers[enhancers == FALSE]  <- 0
enhancers <- enhancers %>%
  dplyr::select(-lead_variant) %>%
  as.matrix()

col_groups <- results_df %>%
  ungroup() %>%
  distinct(sample, annotation) %>%
  arrange(sample) %>%
  dplyr::select(annotation) %>%
  base::as.data.frame()

row.names(col_groups) <- colnames(enhancers)

pheatmap(enhancers, cluster_cols = FALSE, main = "Active Enhancers", annotation_col = col_groups, legend = FALSE,
         annotation_names_col = FALSE) 

results_df %>%
  ungroup() %>%
  filter(state %in% c("Active Promoter", "Active Enhancer")) %>%
  filter(fdr_sig == TRUE) %>%
  dplyr::select(lead_variant, state, sample, fdr_sig, annotation) %>%
  group_by(lead_variant, state) %>%
  summarise(n_sig = sum(fdr_sig),
            tissue = base::paste(base::unique(annotation), collapse = ", ")) %>%
  arrange(n_sig) %>%
  ungroup() %>%
  spread(state, tissue)

results_df %>%
  ungroup() %>%
  filter(state %in% c("Active Promoter", "Active Enhancer")) %>%
  filter(fdr_sig == TRUE) %>%
  dplyr::select(lead_variant, state, sample, fdr_sig, annotation) %>%
  group_by(lead_variant) %>%
  mutate(n_sig = length(unique(annotation))) %>%
  ungroup() %>%
  group_by(lead_variant, state) %>%
  summarise(n_sig = as.numeric(paste(unique(n_sig))),
            tissue = base::paste(base::unique(annotation), collapse = ", ")) %>%
  arrange(n_sig) %>%
  ungroup() %>%
  spread(state, tissue) %>%
  arrange(desc(n_sig)) #%>%
  # write_tsv("significant_tissue.tsv")


results_df %>%
  ungroup() %>%
  filter(state %in% c("Active Promoter", "Active Enhancer")) %>%
  filter(fdr_sig == TRUE) %>%
  group_by(annotation, state) %>%
  summarise(n = length(unique(lead_variant))) %>%
  ungroup() %>%
  spread(state, n) %>%
  arrange(desc(`Active Enhancer`)) #%>%
  #write_tsv("enrihment_by_annotation.tsv")

