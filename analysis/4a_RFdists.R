#' Read in trees
#'
#' @param data_path path to all data files
#' @param prefix prefix for species, mito, and minor trees (must have suffixes "_speciestree.txt", "_mitotree.txt", and "_minortree.txt")
#' @param nmt_gt_file file name for nmt gene trees (rooted)
#' @param cont_gt_file file name for glycolysis gene trees (rooted)
#'
#' @return (multi)phylo objects for each tree read in
#' @export
read_files <- function(data_path, nmt_gt_file, cont_gt_file = NULL, prefix, mito = FALSE){
  nmt_gt <- ape::read.tree(paste0(data_path, nmt_gt_file))
  cont_gt <- ape::read.tree(paste0(data_path, cont_gt_file))
  
  spptree <- ape::read.tree(paste0(data_path, prefix, "_speciestree.txt"))
  mitotree <- ape::read.tree(paste0(data_path, prefix, "_mitotree.txt"))
  minortree <- ape::read.tree(paste0(data_path, prefix, "_minortree.txt"))
  
  return(list(nmt_gt = nmt_gt, cont_gt = cont_gt, spptree = spptree, mitotree = mitotree, minortree = minortree))
}

#' Calculate RF distances; assumes trees are rooted
#' 
#' @param comp_gt gene trees to calculate RF distances from
#' @param comp_to tree to calculate RF distances to
#' @param weighted whether to perform weighted RF distance comparison (requires branch lengths)
#' @param save_tips whether to save tip labels to returned object
#'
#' @return list of RF distances
#' @export
calc_rfs <- function(comp_gt, comp_to, weighted = FALSE, save_tips = FALSE){
  gt_length <- length(comp_gt)
  dists <- integer(length(gt_length))

  for(i in 1:gt_length) {
    tree <- comp_gt[i]
    if (!weighted) dists[i] <- 1 - phangorn::RF.dist(comp_to, tree, normalize = TRUE, check.labels = TRUE, rooted = TRUE)
    if (weighted) dists[i] <- 1 - phangorn::wRF.dist(comp_to, tree, normalize = TRUE, check.labels = TRUE, rooted = TRUE)
  }
  
  if(!save_tips) return(dists)

  if(save_tips) {
    for(i in 1:gt_length) {
      tree <- comp_gt[i]
      if (!weighted) dists[i] <- 1 - phangorn::RF.dist(comp_to, tree, normalize = TRUE, check.labels = TRUE, rooted = TRUE)
      if (weighted) dists[i] <- 1 - phangorn::wRF.dist(comp_to, tree, normalize = TRUE, check.labels = TRUE, rooted = TRUE)
    }
    gts <- 1:length(comp_gt)
    tips <- 
      gts %>% 
      lapply(function(x) {
        tiplabs <- comp_gt[[x]]$tip.label
        return(tiplabs)
      }) %>% 
      dplyr::bind_cols()
    colnames(tips) <- gts
    
    results <- list("dists" = dists, "tiplabs" = tips)
    return(results)
  }
}

#' Calculate average bootstrap support values for gene trees
#'
#' @param gt_dat multiphylo object with gene trees
#'
#' @return
#' @export
calc_avg_bs <- function(gt_dat){
  gt_length <- length(gt_dat)
  avg_bs <- integer(length(gt_length))
  
  for(i in 1:gt_length){
    avg_bs[i] <- mean(as.numeric(gt_dat[[i]]$node.label), na.rm = T)
    if(is.nan(avg_bs[i])) avg_bs[i] <- 0
  }
  
  avg_bs <- as.data.frame(avg_bs) %>% 
    tibble::rownames_to_column(var = "tree_no")
  
  return(avg_bs)
}

#' Helper function to turn calculation objects into tibbles and add order column for joining
#'
#' @param dat output from calc_dists or calc_avg_bs functions; assumes these objects are named appropriately
#'
#' @return
#' @export
turn_tibble_helper <- function(dat){
  results <- tibble::as_tibble(dat) %>% 
    tibble::rownames_to_column(var = "order")
  
  return(results)
}
