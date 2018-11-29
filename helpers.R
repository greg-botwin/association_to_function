# helper functions
## function to perform enrichments in parallel and avoiding loading in large file
track_enrichment <- function(segmentation_urls, variants) {
  # get segmentations
  #Import State
  roadmap_segmentation <- GetSegmentations(files = segmentation_urls)
  seg <- unlist(roadmap_segmentation)
  
  # Human Annotate Segmentation
  # I am guessing/copying with this assignment
  mcols(seg[seg$state %in% c("EAR", "AR", "ARC", "EARC"), ])$state <- "Active Enhancer"
  mcols(seg[seg$state %in% c("EWR", "EWRC"), ])$state <- "Weak Enhancer"
  mcols(seg[seg$state %in% c("PAR", "PARC"), ])$state <- "Active Promoter"
  mcols(seg[seg$state %in% c("PWR", "PWRC"), ])$state <- "Weak Promoter"
  mcols(seg[seg$state %in% c("HET", "SCR", "TRS", "RPS"), ])$state <- "Other"
  
  # perform enrichment
  enrichibd <- CalculateEnrichment(variants = variants,
                                   features = seg,
                                   feature.type = "segmentations", 
                                   CI = 0.95, 
                                   return.overlaps = TRUE)
  tbl_df(enrichibd$enrichment)
  
}

# paste urls
paste_urls<-  function(tracks) {
  paste0(statehub.roadmap.aws, tracks)
}

## function to annotate candidate loci

annotate_loci <- function(variant, chr, position, window_size) {
  index <- variant
  chr_pos <- paste0(chr, ":", position, "-", position)
  window_plus_minus <- window_size/2
  window <- GRanges(chr_pos) + window_plus_minus
  
  
  # get applicable vcf file
  vcf_file <- paste0("/mnt/share6/SHARED_DATASETS/1000_genomes/ALL.chr", chr, 
                     ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
  
  # get variants in window
  vcf_snps <- GetVariantsInWindow(file = vcf_file, position = window,
                                  genome = "hg19", type = "vcf")
  
  # add pop data or l calc
  vcf_snps<- SetPopulation(vcf = vcf_snps, sample_sheet = my.samples)
  
  # calc ld
  vcf_snps <- CalcLD(vcf = vcf_snps, index = index, population = "EUR")
  
  # split fg bg based on ld
  vcf_snps <- SplitVcfLd(vcf = vcf_snps, ld = c(metric = "R.squared",
                                                cutoff = 0.8,
                                                maf = 0.01))
  
  ## perform loci specific enrichment
  enrichment <- CalculateEnrichment(variants = vcf_snps,
                                    features = rsegs,
                                    feature.type = "segmentations", 
                                    CI = 0.95, 
                                    return.overlaps = TRUE)
  
  x <- tbl_df(enrichment$enrichment)
  x <- x %>%
    mutate(n_fg = length(rowRanges(vcf_snps$fg))) %>%
    mutate(n_bg = length(rowRanges(vcf_snps$bg))) %>%
    mutate(index = index)
    return(x)
  
}

# adds q value and sig if < 0.05
add_q_value <- function(result_df) {
  result_df %>%
    mutate(PEP = 1 - probability) %>%
    group_by(state) %>%
    arrange(PEP) %>%
    mutate(qvalue = cummean(PEP)) %>%
    mutate(fdr_sig = if_else(qvalue < 0.01, TRUE, FALSE))
}
