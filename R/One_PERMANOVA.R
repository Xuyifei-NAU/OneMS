#' One step for PERMANOVA and pairwise comparison
#'
#' Conveniently carried out PERMANOVA and pairwise comparison, the result would be writeen in directory  "Adonis result"
#'
#' @param otu The OTUs table, which row are OTUs and column are samples
#' @param metadata The metadata file, a data.frame, should have column of sampleid and group information
#' @param id The sample id, which could also be found in OTUs table
#' @param gpid The column name of the group id
#' @param dis_method The method to construct the distance matrix, including "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis".
#' @param perm The number of permutations times
#' @param p.adj The method to adjust the p-value, including 'bonferroni','holm','hochberg','hommel','fdr','BH','BY','none'
#' @param prefix The prefix of the output file
#'
#' @return A list of PERMANOVA result
#' @export
#'
#' @examples One_PERMANOVA(otu=Protist$POTUs,metadata=Protist$PDesign,id='row.names',gpid = 'trt',dis_method='bray',perm=999,p.adj='BH')

One_PERMANOVA <- function(otu,metadata,id='row.names',gpid = 'trt',dis_method='bray',perm=999,p.adj='BH',prefix=NULL)
{

  dis_list <- c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis")
  p.adj_list <- c('bonferroni','holm','hochberg','hommel','fdr','BH','BY','none' )

  ## Regenerate the experimental design table
  if (!(id %in% c('row.names',colnames(metadata)) & gpid %in% c('row.names',colnames(metadata)))){
    stop('"id" and "gpid" must be one of the cols or rownames of metadata')
  } else if (id %in% 'row.names'){
    group <- data.frame(ID=rownames(metadata),gp=metadata[,gpid])
  } else {
    group <- data.frame(ID=metadata[,id],gp=metadata[,gpid])
  }


  if (!(identical(colnames(otu) ,group[,"ID"]))){
    stop("OTU table and metadata table must have same sampleid")
  } else if(!(dis_method %in% dis_list)){
    stop("dis_method not found")
  } else if(!(p.adj %in% p.adj_list)){
    stop("p.adj method not found")
  } else if (!(is.numeric(perm))){
    stop("perm must be numeric")
  }


  require(vegan)
  otu <- data.frame(t(otu))

## Adonis analysis
  adonis_result <- adonis(otu~group[,2], group, distance = dis_method, permutations =perm)

  otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
  otuput <- cbind(rownames(otuput), otuput)
  names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')

 if(!(dir.exists('Adonis result'))){
   dir.create('Adonis result')
 }

  if(!(file.exists(paste0('Adonis result/',prefix,'Adonis result all.txt')))){
    write.table(otuput, file = paste0('Adonis result/',prefix,'Adonis result all.txt'), row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }

  ## Pairwise comparison

  group_name <-  unique(group[,2])

  adonis_result_two <- NULL
  for (i in 1:(length(group_name) - 1)) {
    for (j in (i + 1):length(group_name)) {
      group_ij <- subset(group, group[,2] %in% c(group_name[i], group_name[j]))
      otu_ij <- otu[group_ij[,1], ]
      adonis_result_otu_ij <- adonis(otu_ij~group_ij[,2], group_ij, permutations = perm, distance = dis_method)	#指定距离测度，基于 999 次置换
      adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), dis_method, unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
    }
  }
  adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
  names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
  subset(group, group$Type =='BS')
  # Add p-value adjust
  adonis_result_two$'Pr (>F)' <- as.numeric(adonis_result_two$'Pr (>F)')
  adonis_result_two <- cbind(adonis_result_two,p.adjust(adonis_result_two$'Pr (>F)', method = p.adj))
  names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)',paste0('P_adj_',p.adj))

  ## Write the result
  if(!(file.exists(paste0('Adonis result/',prefix,'Adonis result two.txt')))){
    write.table(adonis_result_two, file = paste0('Adonis result/',prefix,'Adonis result two.txt'), row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }


  PERMANOVA.result_all <- otuput
  PERMANOVA.result_two <- adonis_result_two
  result <- list(PERMANOVA.result_all,PERMANOVA.result_two)
  return(result)
}




