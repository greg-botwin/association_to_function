#!/usr/bin/env Rscript
## github packages
library(StatePaintR)
library(funciVar)

## bioconductor packages
library(GenomicRanges)
library(VariantAnnotation)

##cran packages
library(tidyverse)
library(readxl)
library(parallel)

## helper functions
source("helpers.R")


#roadmap segemnts query
#roadmap_segmentation <- GetSegmentations(files = segmentation_urls)
#roadmap_segmentations<- unlist(roadmap_segmentation)
#saveRDS(roadmap_segmentations, file = "roadmap_74_tissue_segmentations.RData")

rsegs <- readRDS("roadmap_74_tissue_segmentations.RData")

# sample population data
my.samples <- read.delim("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel", stringsAsFactors = FALSE)
my.samples <- my.samples %>%
  dplyr::select(sample, pop, super_pop, gender)

# Human Annotate Segmentation
mcols(rsegs[rsegs$state %in% c("EAR", "AR", "ARC", "EARC"), ])$state <- "Active Enhancer"
mcols(rsegs[rsegs$state %in% c("EWR", "EWRC"), ])$state <- "Weak Enhancer"
mcols(rsegs[rsegs$state %in% c("PAR", "PARC"), ])$state <- "Active Promoter"
mcols(rsegs[rsegs$state %in% c("PWR", "PWRC"), ])$state <- "Weak Promoter"
mcols(rsegs[rsegs$state %in% c("HET", "SCR", "TRS", "RPS"), ])$state <- "Other"

## query parameters
ibd_snps_by_hd_annotated <- read_excel("data/nature22969-s2.xlsx", sheet = 3)

un_annotated <- ibd_snps_by_hd_annotated %>%
  dplyr::filter(tier2 == "No") %>%
  dplyr::filter(Annotated == FALSE) 

ibd_snps_by_hd <- read_excel("data/nature22969-s2.xlsx", sheet = 1)

loci_to_annotate <- ibd_snps_by_hd %>%
  dplyr::semi_join(., un_annotated, by = c("HD" = "HD", "chr" = "chr", "signal" = "signal")) %>%
  dplyr::select(variant.lead, chr, position.lead) %>%
  dplyr::mutate(window_size = 1e6)

good_snps <- c("rs6426833", "rs6702254", "rs12463658", "rs17229679", 
               "rs9941524", "rs4676408", "rs7711427", "rs7725339",
               "rs62408223", "rs4946717", "rs28701841", "rs1836846", 
               "rs9494844", "rs7468800", "rs10995271", "rs10761648", 
               "rs7915475", "rs630923", "rs79157249", "rs11614178", 
               "rs184788345", "rs72796367", "rs138425259", "rs4807570",
               "rs145530718", "rs2019262", "rs7517810", "rs77981966", "rs11679753", 
               "rs17225380", "rs56167332", "rs3749925", "rs7997823",
               "rs56083426")

# "rs796819847" is an indel
# "chr9:117571294" is not in 1000genome by position or name although 117571293 is, possible indel TNFSF15
# "rs146029108" is an indel and can't calcualte ld
# "rs12722504" is an indel 
# "rs148319899" is in low complexity site similar to rs561645224, screwing up LD calc
# "rs7307562" is in low complextity site similar to rs561645224, screwing up LD calc
# "chr20:43258079" is not in 1000genome by position or name, 43258081 is indel 
# investigate bad snps
bad_snps <- c("rs796819847", "chr9:117571294", "rs146029108",
                 "rs12722504", "rs148319899", "rs7307562",
                 "rs56083426", "chr20:43258079")

good_loci_to_annotate <- loci_to_annotate %>%
  dplyr::filter(variant.lead %in% good_snps)

result_list <- mcmapply(annotate_loci, variant = good_loci_to_annotate$variant.lead,
       chr = good_loci_to_annotate$chr,
       position = good_loci_to_annotate$position.lead,
       window_size = good_loci_to_annotate$window_size,
       mc.cores = 10, SIMPLIFY = FALSE)

result_list <- lapply(result_list, add_q_value)

saveRDS(result_list, file = "specific_snp_enrichment_pfdr_11_28_18.RData")

for (i in names(result_list)) {
  n_fg <- result_list[[paste0(i)]]$n_fg
  n_bg <- result_list[[paste0(i)]]$n_bg
  
  result_list[[paste0(i)]]%>%
    filter(state %in% c("Active Enhancer", "Active Promoter", "Regulatory")) %>%
    ggplot(aes(x = sample, y = difference, color = as.factor(fdr_sig))) +
    geom_point(fill = NA) +
    geom_errorbar(aes(ymin=lower, ymax=upper)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          legend.position="bottom") +
    scale_color_discrete(name="FDR Significant Diff") +
    facet_wrap(~state, ncol = 1) +
    labs(title = paste0("Difference in Beta Distributions Between Foreground and Background SNPs (95% Credible Interval) for ", i, 
                        "\n", n_fg, " SNPs in Foreground and ", n_bg, " SNPs in Background"))
  
  ggsave(filename = paste0("plots/", i, "_func_enrichment_active_reg.png"), width = 11, height = 8.5, units = c("in"))
  
}
