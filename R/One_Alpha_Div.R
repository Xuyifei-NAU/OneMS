#' A convenient function for calculating alpha diversity
#'
#' A convenient function for calculating alpha diversity
#'
#' @param otu The OTUs table, which row are OTUs and column are samples
#' @param tree Tree file
#' @param rrarefy Logical, whether to rrarefy
#' @param melt Logical, whether to melt the result according to the alpha div index
#'
#' @return A data.frame
#' @export
#'
#' @examples One_Alpha_Div(Protist$POTUs)
#'
One_Alpha_Div <- function(otu,tree=NULL,rrarefy=T,melt=F){
  ## packages
  require(vegan)
  require(picante) # pd()
  require(reshape2)


  if(!(is.matrix(otu) | is.data.frame(otu))){
    stop('otu table must be matrix or data.frame')
  } else{
    otu <- as.data.frame(otu)
  }


  if(!(is.logical(rrarefy)&is.logical(melt))){
    stop('rrarefy and melt must be logical')
  }

  ## rrarefy
  if (rrarefy==T){
    otu_t <- t(otu)
    ## Calculate the rrarefy threshold
    raremin <- min(rowSums(otu_t))
    ## rrarefy
    set.seed(888)
    otu_rare <- rrarefy(otu_t,raremin)
  } else{
    otu_rare <- t(otu)
  }

  richness <- t(estimateR(otu_rare))[,c(1)]
  chao1 <- t(estimateR(otu_rare))[,c(2)]
  ace <- t(estimateR(otu_rare))[,c(4)]
  shannon_index <- diversity(otu_rare, index = 'shannon', base = exp(1))
  shannon_diversity <- exp(1)^shannon_index
  simpson <- diversity(otu_rare, index = 'simpson')
  pielou <- diversity(otu_rare, index = 'shannon', base = exp(1)) / log(estimateR(otu_rare)[1, ], exp(1))
  goods_coverage <- 1 - rowSums(otu_rare == 1) / rowSums(otu_rare)

  if (is.null(tree)){
    alpha_diversity <- data.frame(richness,chao1,ace,shannon_index,shannon_diversity,simpson,pielou,goods_coverage)
  } else{
    pd <- pd(otu_rare, tree, include.root = FALSE)
    alpha_diversity <- data.frame(richness,chao1,ace,shannon_index,shannon_diversity,simpson,pielou,goods_coverage,pd)
  }

  if (melt==T){
    alpha_diversity$SampleId <- rownames(alpha_diversity)
    alpha_diversity <- reshape2::melt(alpha_diversity,id.vars = 'SampleId')
  }
  return(alpha_diversity)
}
