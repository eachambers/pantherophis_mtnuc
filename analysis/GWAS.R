#' Helper function to calculate relevant GWAS statistics
#' Adapted from M. Matz code: https://github.com/z0on/Multivariate_GWAS
#' 
#' @param gwas an RData object outputted from GWAS analysis
#'
#' @return gwas df with position in Mb and signed log p-value columns
gwas_stats <- function(gwas) {
  # Identify which Z-scores fall above 0
  sign = as.numeric(gwas$zscore > 0) # returns 1s and 0s for each SNP
  sign[sign == 0] = -1 # switches occurrences of 0s to -1s
  gwas$signed.logp = gwas$logp*sign # assign signs to log p-values based on Z-score results
  gwas$pos.Mb = gwas$pos/1e+6 # convert bp to Mb
  
  return(gwas)
}

#' Function to visualize outlier SNPs occurring in certain regions of the genome (i.e., N-mts)
#'
#' @param gwas resulting object from GWAS analysis (GWAS_analysis.R script)
#' @param regions regions we're interested in detecting outliers within (same format as gff)
#' @param sig significance threshold (alpha value) for detecting outliers (defaults to 0.05)
#' @param ld_blocks whether to also display SNPs occurring within user-specified distance from specified region (defaults to FALSE)
#' @param ld_dist if `ld_blocks = TRUE`, number of nucs to specify distance from specified region in either direction
#' @param manhattan whether to output Manhattan plot with SNPs colorized according to whether they are in specified region
#' @param color_block if `manhattan = TRUE`, whether to highlight regions that are specified regions on Manhattan plot
#' @param subdivide if `manhattan = TRUE`, whether to subdivide Manhattan plot (defaults to FALSE)
#' @param subdivide_length if `subdivide = TRUE`, what length to display per Manhattan plot (in Mb; defaults to 10)
#'
#' @return table with results
gwas_vis <- function(gwas, regions, sig = 0.05, ld_blocks = FALSE, ld_dist = NULL, manhattan = TRUE, color_block = TRUE, subdivide = FALSE, subdivide_length = 10) {
  # Gather relevant stats
  gwas <- gwas_stats(gwas)
  
  # Extract SNPs above significance threshold
  sig_snps <- gwas %>% 
    dplyr::filter(signed.logp > -log10(sig))

  join <-
    fuzzyjoin::fuzzy_inner_join(sig_snps, regions, by = c("chrom" = "chrom",
                                               "pos" = "start",
                                               "pos" = "end"),
                   match_fun = list(`==`, `>=`, `<=`))
  
  # Make table with results from significant SNPs only
  gwas_table(join)
  
  if(manhattan){
    gwas_manhattan(gwas = gwas, regions = regions, sig = sig, color_block = color_block, subdivide = TRUE, subdivide_length = 10)
  }
  return(join)
}

#' Build table of GWAS results
#' Adapted from rda_table() function in algatr package
#' 
#' @param gwas resulting object from GWAS analysis (RDA_GWAS.R script)
#' @param sig significance threshold (defaults to 0.05)
#' @param sig_only whether to only include loci with p-values less than `sig` (defaults to TRUE)
#' @param top whether to only include only keep the top variable for each snp in the table by the strength of the correlation (defaults to FALSE)
#' @param order whether to order by the magnitude of the correlation (defaults to FALSE)
#' @param nrow number of rows to display (defaults to displaying all rows)
#' @param digits number of digits to include (defaults to 2)
#'
#' @return prints table
gwas_table <- function(gwas, sig = 0.05, sig_only = TRUE, top = FALSE, order = FALSE, nrow = NULL, digits = 2) {
  
  gwas <- gwas_stats(gwas)
  
  if(sig_only) gwas <- gwas %>% 
      dplyr::filter(signed.logp > -log10(sig))
                              
  if(order) gwas <- gwas[order(abs(gwas$signed.logp), decreasing = TRUE),]
  if(top) gwas <- gwas %>%
      dplyr::group_by(snp) %>%
      dplyr::filter(abs(r) == max(abs(r)))
  if(!is.null(nrow)) {
    if(nrow > nrow(gwas)) nrow <- nrow(gwas)
    gwas <- gwas[1:nrow, ]
  }
  
  gwas <- gwas %>% dplyr::as_tibble()
  if(!is.null(digits)) gwas <- gwas %>% dplyr::mutate(dplyr::across(-c(var, snp), round, digits))
  
  # Get min/max values
  d <- max(abs(min(gwas$signed.logp)), abs(max(gwas$signed.logp)))
  
  suppressWarnings(
    # tbl <- 
    gwas %>%
      dplyr::select(chrom, pos, zscore, r2, signed.logp) %>% 
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric(signed.logp, trim = TRUE, domain = c(-d,d)))
  
  tbl
}

#' Build Manhattan plot of outliers within nmts
#'
#' @inheritParams gwas_vis function
#'
#' @return
gwas_manhattan <- function(gwas, regions, sig, color_block, subdivide = TRUE, subdivide_length = 10) {
  # Extract relevant values from which there are SNPs
  nmt_sub <-
    nmt_sub %>% 
    mutate(start_mb = start/1e+6,
           end_mb = end/1e+6) %>% 
    filter(start < max_pos)
  
  # Extract relevant values for plotting
  pmax <- max(gwas$signed.logp) + 0.05*(max(gwas$signed.logp))
  max_pos <- max(gwas$pos)
  
  # Build plot --------------------------------------------------------------
  p <-
    gwas %>% 
    filter(pos.Mb <= 10) %>% 
    ggplot(aes(x = pos.Mb, y = signed.logp)) +
    geom_point(size = 0.25, alpha = 0.25, col="grey") +
    ylab("-log10(p-value)") +
    xlab("position (Mb)") +
    theme_cowplot() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.6) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))

  p +
    geom_rect(data = nmt_divide, aes(xmin = nmt_divide$start_mb, xmax = nmt_divide$end_mb, ymin = 0, ymax = pmax),
              fill = "lightgreen", alpha = 0.5)
}
