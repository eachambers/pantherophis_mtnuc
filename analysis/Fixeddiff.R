#' Process sequence data
#'
#' @param vcf_file full path to vcf file
#' @param fasta_file full path to fasta file
#' @param mitotypes object with mitotypes (two cols, INDV and mitotype)
#'
#' @return df with each SNP in its own col, individual names, and mitotype
#' @export
process_seq <- function(vcf_file, fasta_file, mitotypes) {
  vcf <- vcfR::read.vcfR(vcf_file)
  
  # Get ordering of inds from vcf
  names <- colnames(vcf@gt[,-1])
  names <- as.data.frame(names) %>% 
    separate(names, into = c("INDV", "sorted", "dedup", "filetype"), sep = "_") %>% 
    dplyr::select(-sorted, -dedup, -filetype)
  
  # Get chromosome and position IDs from vcf
  chroms <- as.data.frame(vcf@fix[,1]) %>% 
    dplyr::rename(CHROM = 1)
  pos <- as.data.frame(vcf@fix[,2]) %>% 
    dplyr::rename(POS = 1)
  loci <- cbind(chroms, pos) %>% 
    unite(locus_name, c("CHROM", "POS"), sep = "_")
  
  sequences = seqinr::read.fasta(fasta_file)
  results <- do.call(rbind, sequences)
  results <- toupper(results)
  
  # Add in loci names, ind names, and mitotypes
  colnames(results) <- loci$locus_name
  results <- cbind(names, results)
  results <- left_join(results, mitotypes)
  
  return(results)
}

#' Calculate pure allele dictionary (i.e., fixed differences between reference samples)
#'
#' @param dataset either cont or nmts for control or NMT gene datasets, respectively; object output from `process_seq()` function
#' @param dataset_name either "cont" or "nmts"; must correspond to `dataset`
#' @param save_file whether to save data (as .rda)
#' @param output_path if `save_file = TRUE`, path to save output files
#'
#' @return pure allele dictionary with no ambiguous bases included, `pure_allele_dict_noambig`
#' @export
allele_dict <- function(dataset, dataset_name, save_file = TRUE, output_path) {
  dataset %>%
    # Only consider two reference samples, H21189 (slowinskii) and O42736 (emoryi)
    filter(INDV == "H21189" | INDV == "O42736") %>%
    # Make into long-form dataset, ignoring sample ID and mitotype columns
    pivot_longer(names_to = "locus", values_to = "value", cols = -c(INDV, mitotype)) %>% # 77,102 rows
    # Remove missing values
    filter(value != 'N') %>%                                    # 76,440 rows
    dplyr::count(mitotype, locus, value) %>%
    group_by(mitotype, locus) %>%
    mutate(frac = n/sum(n)) %>%
    ungroup() %>%
    filter(frac >= 1) %>%
    dplyr::select(mitotype, locus, value) %>%
    pivot_wider(names_from = mitotype, values_from = value) %>%
    filter(emoryi != slowinskii) -> pure_allele_dict
  
  # Remove ambiguous bases from the pure allele dictionary
  amb_bases = c("M", "W", "K", "R", "Y", "S")
  # If you want to remove ambiguous bases:
  pure_allele_dict_noambig <-
    pure_allele_dict %>% 
    dplyr::filter(!emoryi %in% amb_bases) %>% 
    dplyr::filter(!slowinskii %in% amb_bases)
  
  # Add a heterozygous state to the dictionary
  pure_allele_dict_noambig <-
    pure_allele_dict_noambig %>% 
    mutate(het = case_when(
      emoryi == "A" & slowinskii == "C" | emoryi == "C" & slowinskii == "A" ~ "M",
      emoryi == "A" & slowinskii == "T" | emoryi == "T" & slowinskii == "A" ~ "W",
      emoryi == "A" & slowinskii == "G" | emoryi == "G" & slowinskii == "A" ~ "R",
      emoryi == "C" & slowinskii == "T" | emoryi == "T" & slowinskii == "C" ~ "Y",
      emoryi == "C" & slowinskii == "G" | emoryi == "G" & slowinskii == "C" ~ "S",
      emoryi == "G" & slowinskii == "T" | emoryi == "T" & slowinskii == "G" ~ "K"))
  
  if (save_file) save(pure_allele_dict_noambig, file = paste0(output_path, "pure_allele_dict_", dataset_name, ".rda"))
  print(paste0("There were ", nrow(pure_allele_dict_noambig), " fixed differences found"))
  
  return(pure_allele_dict_noambig)
}

#' Calculate per-locus allele frequencies
#'
#' @param dataset either cont or nmts for control or NMT gene datasets, respectively; object output from `process_seq()` function
#' @param dataset_name either "cont" or "nmts"; must correspond to `dataset`
#' @param pure_allele_dict object output from `allele_dict()` function of fixed diffs between reference samples
#' @param threshold threshold number of individuals to be considered assigned as a match or mismatch
#' @param save_file whether to save data (as .rda) and summary stats (as .txt)
#' @param output_path if `save_file = TRUE`, path to save output files
#'
#' @return
#' @export
freq_data <- function(dataset, dataset_name, pure_allele_dict, threshold, save_file = TRUE, output_path) {
  dataset %>% 
    pivot_longer(names_to = "locus", values_to = "value", cols = -c(INDV, mitotype)) %>%
    filter(value != 'N') %>% 
    inner_join(pure_allele_dict, by = "locus") %>% 
    group_by(mitotype, locus) %>% 
    summarize(frac_emoryi = sum(value == emoryi)/n(),
              frac_slow = sum(value == slowinskii)/n(),
              frac_het = sum(value == het)/n()) -> freqs
  
  summary <-
    freqs %>% 
    summarize(emoryi_match = sum(mitotype == "emoryi" & frac_emoryi >= threshold),
              emoryi_mismatch = sum(mitotype == "emoryi" & frac_slow >= threshold),
              emory_het = sum(mitotype == "emoryi" & frac_het >= threshold),
              slow_match = sum(mitotype == "slowinskii" & frac_slow >= threshold),
              slow_mismatch = sum(mitotype == "slowinskii" & frac_emoryi >= threshold),
              slow_het = sum(mitotype == "slowinskii" & frac_het >= threshold),
              dataset = dataset_name,
              threshold = threshold)
  
  if (save_file) {
    save(freqs, file = paste0(output_path, "freq_data_", dataset_name, "_", threshold, ".rda"))
    write_tsv(summary, file = paste0(output_path, "summary_", dataset_name, "_", threshold, ".txt"))
  }
  
  return(list(freqs = freq_data, summary = summary))
}