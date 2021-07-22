#' One step for ANOSIM and pairwise comparison
#'
#' Conveniently carried out ANOSIM and pairwise comparison, the result would be writeen in directory  "ANOSIM result"
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
#' @return A list of ANOSIM result
#' @export
#'
#' @examples One_ANOSIM(otu=Protist$POTUs,metadata=Protist$PDesign,id='row.names',gpid = 'trt',dis_method='bray',perm=999,p.adj='BH')
#'
#'
One_ANOSIM <- function(otu,metadata,id='row.names',gpid = 'trt',dis_method='bray',perm=999,p.adj='BH',prefix=NULL)
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

  anosim_result <- anosim(otu, group[,2], distance = 'bray', permutations = perm)

  output <- NULL
  output <- cbind(anosim_result$permutations,dis_method,anosim_result$signif,anosim_result$statistic)
  anosim_result$signif	#p value
  anosim_result$statistic	#R value
  colnames(output) <- c('permutations','distance','P','R')

  if(!(dir.exists('Anosim result'))){
    dir.create('Anosim result')
  }

  if(!(file.exists(paste0('Anosim result/',prefix,'Anosim result all.txt')))){
    write.table(output, file = paste0('Anosim result/',prefix,'Anosim result all.txt'), row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }

  ## Pairwise comparison

  group_name <- unique(group[,2])

  anosim_result_two <- NULL
  for (i in 1:(length(group_name) - 1)) {
    for (j in (i + 1):length(group_name)) {
      group_ij <- subset(group, group[,2] %in% c(group_name[i], group_name[j]))
      otu_ij <- otu[group_ij[,1], ]

      anosim_result_otu_ij <- anosim(otu_ij, group_ij[,2], permutations = perm, distance = dis_method)	#Bray-Curtis 距离测度，基于 999 次置换
      anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = '/'), dis_method, anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif))
    }
  }

  anosim_result_two <- data.frame(anosim_result_two, stringsAsFactors = FALSE)
  names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value')



  # Add p-value adjust
  anosim_result_two$P_value <- as.numeric(anosim_result_two$P_value)
  anosim_result_two <- cbind(anosim_result_two,p.adjust(anosim_result_two$P_value, method = p.adj))
  names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value',paste0('P_adj_',p.adj))


  if(!(file.exists(paste0('Anosim result/',prefix,'Anosim result two.txt')))){
    write.table(anosim_result_two, file = paste0('Anosim result/',prefix,'Anosim result two.txt'), row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }

  result <- list(anosim_result,anosim_result_two)
  return(result)
}
