IBD Sub-Type Associated SNPs Chromatin Functional Enrichment
================
Translational Genomics Group
29 November, 2018

``` r
## github packages
library(StatePaintR)
library(funciVar)

##cran packages
library(tidyverse)
library(readxl)

## bioconductor packages
library(GenomicRanges)
library(VariantAnnotation)

## helper functions
source("helpers.R")
```

Genome Segmentations
--------------------

StateHub was browsed for relevant segementations.

All NIH Roadmap Consortium tracks annotated with the following marks excluding cells lines and stem cell lines were chosen for analysis. These were largely similar to the tracks chosen in the BFG Parkinson's Disease candidate SNP analysis.

All segments were annotated with at least the following 6 histone marks:

-   H3K27ac
-   H3K27me3
-   H3K36me3
-   H3K4me1
-   H3K4me3
-   H3K9me3

The Focused Poised Promoter Model (TrackHub ID5813b67f46e0fb06b493ceb0).

``` r
# statehub url for poised model + roadmap
statehub.roadmap.aws <- "http://s3-us-west-2.amazonaws.com/statehub-trackhub/tracks/5813b67f46e0fb06b493ceb0/hg19/Roadmap/"

# 89 Roadmap Tissue Tracks Anoted with >=6 Marks
tracks <- c("adrl.glnd.fet.7mark.segmentation.bed",
            "bld.cd14.mono.7mark.segmentation.bed",
            "bld.cd14.pc.7mark.segmentation.bed",
            "bld.cd19.ppc.7mark.segmentation.bed",
            # "bld.cd3.ppc.7mark.segmentation.bed", error loading removing from analysis
            "bld.cd56.pc.7mark.segmentation.bed",
            "bld.mob.cd34.pc.f.7mark.segmentation.bed",
            "brn.nha.7mark.segmentation.bed",
            "brst.hmec.7mark.segmentation.bed",
            "gi.l.int.fet.7mark.segmentation.bed",
            "gi.s.int.7mark.segmentation.bed",
            "gi.s.int.fet.7mark.segmentation.bed",
            "gi.stmc.fet.7mark.segmentation.bed",
            "gi.stmc.gast.7mark.segmentation.bed",
            "lng.nhlf.7mark.segmentation.bed",
            "mus.hsmm.7mark.segmentation.bed",
            "mus.hsmmt.7mark.segmentation.bed",
            "mus.leg.fet.7mark.segmentation.bed",
            "mus.psoas.7mark.segmentation.bed",
            "mus.trnk.fet.7mark.segmentation.bed",
            "ovry.7mark.segmentation.bed",
            "panc.7mark.segmentation.bed",
            "plcnt.fet.7mark.segmentation.bed",
            "skin.nhdfad.7mark.segmentation.bed",
            "skin.nhek.7mark.segmentation.bed",
            "skin.pen.frsk.fib.01.7mark.segmentation.bed",
            "skin.pen.frsk.fib.02.7mark.segmentation.bed",
            "skin.pen.frsk.mel.01.7mark.segmentation.bed",
            "thym.fet.7mark.segmentation.bed",
            "vas.huvec.7mark.segmentation.bed",
            "bld.cd4.cd25.cd127m.tregpc.6mark.segmentation.bed", 
            "bld.cd4.cd25i.cd127.tmempc.6mark.segmentation.bed",
            "bld.cd4.cd25m.cd45ra.npc.6mark.segmentation.bed",
            "bld.cd4.cd25m.cd45ro.mpc.6mark.segmentation.bed", 
            "bld.cd4.cd25m.il17m.pl.tpc.6mark.segmentation.bed",
            "bld.cd4.cd25m.il17p.pl.tpc.6mark.segmentation.bed",
            "bld.cd4.cd25m.tpc.6mark.segmentation.bed",
            "bld.cd4.mpc.6mark.segmentation.bed",
            "bld.cd4.npc.6mark.segmentation.bed",
            "bld.cd8.mpc.6mark.segmentation.bed",
            "bld.cd8.npc.6mark.segmentation.bed",
            "bld.per.monuc.pc.6mark.segmentation.bed",
            "bone.osteo.6mark.segmentation.bed",
            "brn.ang.gyr.6mark.segmentation.bed",
            "brn.ant.caud.6mark.segmentation.bed",
            "brn.cing.gyr.6mark.segmentation.bed",
            "brn.dl.prfrntl.crtx.6mark.segmentation.bed",
            "brn.hipp.mid.6mark.segmentation.bed",
            "brn.inf.tmp.6mark.segmentation.bed",
            "brn.sub.nig.6mark.segmentation.bed",
            "fat.adip.nuc.6mark.segmentation.bed",
            "gi.cln.muc.6mark.segmentation.bed",
            "gi.cln.sig.6mark.segmentation.bed",
            "gi.cln.sm.mus.6mark.segmentation.bed",
            "gi.duo.sm.mus.6mark.segmentation.bed",
            "gi.eso.6mark.segmentation.bed",
            "gi.rect.muc.29.6mark.segmentation.bed",
            "gi.rect.muc.31.6mark.segmentation.bed",
            "gi.rect.sm.mus.6mark.segmentation.bed",
            "gi.stmc.mus.6mark.segmentation.bed",
            "hrt.atr.r.6mark.segmentation.bed",
            "hrt.vent.l.6mark.segmentation.bed",
            "hrt.vnt.r.6mark.segmentation.bed",
            "liv.adlt.6mark.segmentation.bed",
            "lng.6mark.segmentation.bed",
            "mus.sklt.f.6mark.segmentation.bed",
            "panc.islt.6mark.segmentation.bed",
            "plcnt.amn.6mark.segmentation.bed",
            "skin.pen.frsk.ker.03.6mark.segmentation.bed",
            "skin.pen.frsk.mel.03.6mark.segmentation.bed",
            "spln.6mark.segmentation.bed",
            "strm.chon.mrw.dr.msc.6mark.segmentation.bed",
            "strm.mrw.msc.6mark.segmentation.bed",
            "thym.6mark.segmentation.bed",
            "vas.aor.6mark.segmentation.bed"
            )

segmentation_urls<- unlist(lapply(tracks, paste_urls))

#roadmap_segmentation <- GetSegmentations(files = segmentation_urls)
#rsegs <- unlist(roadmap_segmentation)
#saveRDS(rsegs, file = "roadmap_74_tissue_segmentations.RData")
```

Designating Foreground and Background SNPs from IBD Fine Mapping Study
----------------------------------------------------------------------

Genetic data from 67k subejcts all gentoyped on iChip was obtained. The iChip array can be brokedn down into 187 high density regions. Previous IBD genome wide associations were mapped to a high density region, imputation was performed, and fine mapping was performed.

"We applied three complementary Bayesian fine-mapping methods that used different priors and model selection strategies to identify independent association signals within a region, and to assign a posterior probability of causality to each variant (Supplementary Methods and Extended Data Fig. 2). For each independent signal detected by each method, we sorted all variants by the posterior probability of association, and added variants to the ‘credible set’ of associated variants until the sum of their posterior probability exceeded 95%: that is, the credible set contained the minimum list of DNA variants that were &gt;95% likely to contain the causal variant (Fig. 1). These sets ranged in size from 1 to &gt;400 variants."

### Total SNPs

``` r
## Get Background and Goreground SNPs from Fine Mapping Paper
all_snps <- read_tsv("data/filtered_imputed_snps_per_block_with_dups.tsv")

### Need to add Chr to blocks
ichip_hd_blocks <- read_tsv("data/ichip_regions_mt2_HG19-2.txt", 
                            col_names = c("Chr", " Start", "End", "block", "Disease"))

ichip_hd_blocks <- ichip_hd_blocks %>%
  dplyr::select(Chr, block)

all_snps <- left_join(all_snps, ichip_hd_blocks, by = "block")

### Remove Dups Because of extending the region
all_snps <- all_snps %>%
  dplyr::select(rs_id, position, Chr, block) %>%
  distinct(rs_id, position, Chr, .keep_all = TRUE)

all_snps
```

#### Fix All Snp Issues

``` r
### Get list of possible FG snps
fg <- read_excel("data/nature22969-s2.xlsx", sheet = 2)
fg_snps <- fg %>%
  dplyr::select(variant, chr, position, HD, signal) %>%
  distinct()

problem_snps <- fg_snps %>%
  filter(!variant %in% all_snps$rs_id)

# list of snps i can't uniquely assign 
non_assignable <- c("imm_5_40408209", "rs59418409", "rs75900575", "rs2823259")

problem_snps_fixed <- problem_snps %>% 
  filter(!variant %in% non_assignable) %>%
  left_join(., all_snps, by = c("chr" = "Chr", "position" = "position")) %>%
  dplyr::select(-block) %>%
  mutate(variant = rs_id) %>%
  dplyr::select(-rs_id)

fg_snps <- fg_snps %>%
  filter(!variant %in% problem_snps$variant) %>%
  bind_rows(., problem_snps_fixed) %>%
  mutate(chr = paste0("chr", chr))

table(fg_snps$variant %in% all_snps$rs_id)
```

CD Specific Loci
----------------

``` r
ibd_snps_by_hd <- read_excel("data/nature22969-s2.xlsx", sheet = 1)

# restrict to 139 loci that passed fine mapping criteria
ibd_snps_by_hd <- ibd_snps_by_hd %>%
  filter(tier2 == "No")

# number of signals assinged to IBD sub-type
table(ibd_snps_by_hd$trait.reassigned)

cd_signals <- ibd_snps_by_hd %>%
  filter(trait.reassigned == "CD") %>%
  dplyr::select(HD, signal)

uc_signals <- ibd_snps_by_hd %>%
  filter(trait.reassigned == "UC") %>%
  dplyr::select(HD, signal)

ibd_signals <- ibd_snps_by_hd %>%
  filter(trait.reassigned == "IBD") %>%
  dplyr::select(HD, signal)
```

``` r
cd_fg_snps <- inner_join(fg_snps, cd_signals, by = c("HD" = "HD", "signal" = "signal"))


# create fg object

cd_fg <- GRanges(seqnames = cd_fg_snps$chr,
              ranges = IRanges(start = cd_fg_snps$position,
                               end = cd_fg_snps$position,
                               width = 1,
                               names = cd_fg_snps$variant))

cd_fg
```

The background SNPs are all SNPs in the implicated HD regions. IN HD regions, where there are multiple signals, the background SNPs will consist of SNPs that were assigned to UC and IBD.

``` r
cd_bg_snps <- all_snps %>%
  #filter(!rs_id %in% cd_fg_snps$variant) %>%
  separate(block, into = c("prefix", "HD"), sep = 2) %>%
  dplyr::select(-prefix)%>%
  filter(HD %in% cd_fg_snps$HD) %>%
  mutate(Chr = paste0("chr", Chr))

cd_bg <- GRanges(seqnames = cd_bg_snps$Chr,
              ranges = IRanges(start = cd_bg_snps$position,
                               end = cd_bg_snps$position,
                               width = 1,
                               names = cd_bg_snps$rs_id))

cd_bg

cd_snps <- list(fg = cd_fg, bg = cd_bg)
```

UC Specifci Loci
----------------

``` r
uc_fg_snps <- inner_join(fg_snps, uc_signals, by = c("HD" = "HD", "signal" = "signal"))

# create fg object

uc_fg <- GRanges(seqnames = uc_fg_snps$chr,
              ranges = IRanges(start = uc_fg_snps$position,
                               end = uc_fg_snps$position,
                               width = 1,
                               names = uc_fg_snps$variant))

uc_fg
```

The background SNPs are all SNPs in the implicated HD regions. IN HD regions, where there are multiple signals, the background SNPs will consist of SNPs that were assigned to CD and IBD.

``` r
uc_bg_snps <- all_snps %>%
  #filter(!rs_id %in% uc_fg_snps$variant) %>%
  separate(block, into = c("prefix", "HD"), sep = 2) %>%
  dplyr::select(-prefix)%>%
  filter(HD %in% uc_fg_snps$HD) %>%
  mutate(Chr = paste0("chr", Chr))

uc_bg <- GRanges(seqnames = uc_bg_snps$Chr,
              ranges = IRanges(start = uc_bg_snps$position,
                               end = uc_bg_snps$position,
                               width = 1,
                               names = uc_bg_snps$rs_id))

uc_bg

uc_snps <- list(fg = uc_fg, bg = uc_bg)
```

### IBD Shared Loci

``` r
ibd_fg_snps <- inner_join(fg_snps, ibd_signals, by = c("HD" = "HD", "signal" = "signal"))


# create fg object

ibd_fg <- GRanges(seqnames = ibd_fg_snps$chr,
              ranges = IRanges(start = ibd_fg_snps$position,
                               end = ibd_fg_snps$position,
                               width = 1,
                               names = ibd_fg_snps$variant))

ibd_fg
```

The background SNPs are all SNPs in the implicated HD regions. IN HD regions, where there are multiple signals, the background SNPs will consist of SNPs that were assigned to CD and UC, but not shared in IBD.

``` r
ibd_bg_snps <- all_snps %>%
  #filter(!rs_id %in% ibd_fg_snps$variant) %>%
  separate(block, into = c("prefix", "HD"), sep = 2) %>%
  dplyr::select(-prefix)%>%
  filter(HD %in% ibd_fg_snps$HD) %>%
  mutate(Chr = paste0("chr", Chr))

ibd_bg <- GRanges(seqnames = ibd_bg_snps$Chr,
              ranges = IRanges(start = ibd_bg_snps$position,
                               end = ibd_bg_snps$position,
                               width = 1,
                               names = ibd_bg_snps$rs_id))

ibd_bg

ibd_snps <- list(fg = ibd_fg, bg = ibd_bg)
```

### Load Segs

``` r
rsegs <- readRDS("roadmap_74_tissue_segmentations.RData")
# Human Annotate Segmentation
mcols(rsegs[rsegs$state %in% c("EAR", "AR", "ARC", "EARC"), ])$state <- "Active Enhancer"
mcols(rsegs[rsegs$state %in% c("EWR", "EWRC"), ])$state <- "Weak Enhancer"
mcols(rsegs[rsegs$state %in% c("PAR", "PARC"), ])$state <- "Active Promoter"
mcols(rsegs[rsegs$state %in% c("PWR", "PWRC"), ])$state <- "Weak Promoter"
mcols(rsegs[rsegs$state %in% c("HET", "SCR", "TRS", "RPS"), ])$state <- "Other"
```

### Perform CD Specific Enrichment

``` r
# perform enrichment

cd_enrichment <-CalculateEnrichment(variants = cd_snps,
                                   features = rsegs,
                                   feature.type = "segmentations", 
                                   CI = 0.95, 
                                   return.overlaps = TRUE)

cd_enrichment <- tbl_df(cd_enrichment$enrichment) %>%
  mutate(n_fg = length(cd_snps$fg)) %>%
  mutate(n_bg = length(cd_snps$bg)) %>%
  mutate(sub_type = "CD") %>% 
  mutate(PEP = 1 - probability) %>%
  group_by(state) %>%
  arrange(PEP) %>%
  mutate(qvalue = cummean(PEP)) %>%
  mutate(fdr_sig = if_else(qvalue < 0.01, TRUE, FALSE))
```

### Perform UC Specific Enrichment

``` r
uc_enrichment <-CalculateEnrichment(variants = uc_snps,
                                   features = rsegs,
                                   feature.type = "segmentations", 
                                   CI = 0.95, 
                                   return.overlaps = TRUE)

uc_enrichment <- tbl_df(uc_enrichment$enrichment) %>%
  mutate(n_fg = length(uc_snps$fg)) %>%
  mutate(n_bg = length(uc_snps$bg)) %>%
  mutate(sub_type = "UC") %>%
  mutate(PEP = 1 - probability) %>%
  group_by(state) %>%
  arrange(PEP) %>%
  mutate(qvalue = cummean(PEP)) %>%
  mutate(fdr_sig = if_else(qvalue < 0.01, TRUE, FALSE))
```

### Perform Shared IBD Loci Enrichment

``` r
ibd_enrichment <-CalculateEnrichment(variants = ibd_snps,
                                   features = rsegs,
                                   feature.type = "segmentations", 
                                   CI = 0.95, 
                                   return.overlaps = TRUE)

ibd_enrichment <- tbl_df(ibd_enrichment$enrichment) %>%
  mutate(n_fg = length(ibd_snps$fg)) %>%
  mutate(n_bg = length(ibd_snps$bg)) %>%
  mutate(sub_type = "Shared IBD SNPs") %>%
  mutate(PEP = 1 - probability) %>%
  group_by(state) %>%
  arrange(PEP) %>%
  mutate(qvalue = cummean(PEP)) %>%
  mutate(fdr_sig = if_else(qvalue < 0.01, TRUE, FALSE))
```

Plot Everything Maybe
---------------------

``` r
sub_type_results <- bind_rows(cd_enrichment, uc_enrichment, ibd_enrichment)

sub_type_results %>%
  filter(state == "Active Enhancer") %>%
  ggplot(aes(x = sample, y = difference, color = as.factor(fdr_sig))) +
  geom_point(fill = NA) +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        legend.position="bottom") +
  scale_color_discrete(name="Significant Diff") +
  facet_wrap(~sub_type, ncol = 1) +
  labs(title = "Difference in Beta Distributions Between Foreground and Background SNPs (95% CI) in Active Enhancer Regions")

ggsave(filename = "plots/ibd_sub_type_func_enrichment_act_enhancer.png", width = 11, height = 8.5, units = c("in"))

sub_type_results %>%
  filter(state == "Active Promoter") %>%
  ggplot(aes(x = sample, y = difference, color = as.factor(fdr_sig))) +
  geom_point(fill = NA) +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        legend.position="bottom") +
  scale_color_discrete(name="Significant Diff") +
  facet_wrap(~sub_type, ncol = 1) +
  labs(title = "Difference in Beta Distributions Between Foreground and Background SNPs (95% CI) in Active Promoter Regions")

ggsave(filename = "plots/ibd_sub_type_func_enrichment_act_promoter.png", width = 11, height = 8.5, units = c("in"))

sub_type_results %>%
  filter(state == "Weak Enhancer") %>%
  ggplot(aes(x = sample, y = difference, color = as.factor(fdr_sig))) +
  geom_point(fill = NA) +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        legend.position="bottom") +
  scale_color_discrete(name="Significant Diff") +
  facet_wrap(~sub_type, ncol = 1) +
  labs(title = "Difference in Beta Distributions Between Foreground and Background SNPs (95% CI) in Weak Enhancer Regions")

ggsave(filename = "plots/ibd_sub_type_func_enrichment_weak_enhancer.png", width = 11, height = 8.5, units = c("in"))

sub_type_results %>%
  filter(state == "Weak Promoter") %>%
  ggplot(aes(x = sample, y = difference, color = as.factor(fdr_sig))) +
  geom_point(fill = NA) +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        legend.position="bottom") +
  scale_color_discrete(name="Significant Diff") +
  facet_wrap(~sub_type, ncol = 1) +
  labs(title = "Difference in Beta Distributions Between Foreground and Background SNPs (95% CI) in Weak Promoter Regions")

ggsave(filename = "plots/ibd_sub_type_func_enrichment_weak_promoter.png", width = 11, height = 8.5, units = c("in"))
```

``` r
table(sub_type_results$state)
```
