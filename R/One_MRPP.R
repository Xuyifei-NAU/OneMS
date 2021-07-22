#' One step for MRPP and pairwise comparison
#'
#' Conveniently carried out MRPP and pairwise comparison, the result would be writeen in directory  "MRPP result"
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
#' @return A list of MRPP result
#' @export
#'
#' @examples One_MRPP(otu=Protist$POTUs,metadata=Protist$PDesign,id='row.names',gpid = 'trt',dis_method='bray',perm=999,p.adj='BH')
One_MRPP <- function(otu,metadata,id='row.names',gpid = 'trt',dis_method='bray',perm=999,p.adj='BH',prefix=NULL)
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

  library(vegan)
  otu <- data.frame(t(otu))

  # MRPP analysis
  mrpp_result <- mrpp(otu, group[,2], distance = dis_method, permutations =perm)

  if(!(dir.exists('MRPP result'))){
    dir.create('MRPP result')
  }

  if(!(file.exists(paste0('MRPP result/',prefix,'MRPP result all.txt')))){
    write.table(data.frame(Group = 'all', Distance = dis_method, A = mrpp_result$A,
                           Observe_delta = mrpp_result$delta, Expect_delta = mrpp_result$E.delta, P_value = mrpp_result$Pvalue)
                , file = paste0('MRPP result/',prefix,'MRPP result all.txt'), row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }



  ## Pairwise comparison
  group_name <- unique(group[,2])

  mrpp_result_two <- NULL
  for (i in 1:(length(group_name) - 1)) {
    for (j in (i + 1):length(group_name)) {
      group_ij <- subset(group, group[,2] %in% c(group_name[i], group_name[j]))
      otu_ij <- otu[group_ij[,1], ]
      mrpp_result_otu_ij <- mrpp(otu_ij, group_ij[,2], permutations = perm, distance = dis_method)	#Bray-Curtis ???????ȣ????? 999 ???û?
      mrpp_result_two <- rbind(mrpp_result_two, c(paste(group_name[i], group_name[j], sep = '/'), dis_method, mrpp_result_otu_ij$A, mrpp_result_otu_ij$delta, mrpp_result_otu_ij$E.delta, mrpp_result_otu_ij$Pvalue))
    }
  }
  mrpp_result_two <- data.frame(mrpp_result_two, stringsAsFactors = FALSE)
  names(mrpp_result_two) <- c('group', 'distance', 'A', 'Observe_delta', 'Expect_delta', 'P_value')

  # 加上校正的p值
  mrpp_result_two$P_value <- as.numeric(mrpp_result_two$P_value)
  mrpp_result_two <- cbind(mrpp_result_two,p.adjust(mrpp_result_two$'P_value', method = p.adj))
  names(mrpp_result_two) <- c('group', 'distance', 'A', 'Observe_delta', 'Expect_delta', 'P_value',paste0('P_adj_',p.adj))

  result <- list(mrpp_result,mrpp_result_two)

  ## Write the result
  if(!(file.exists(paste0('MRPP result/',prefix,'MRPP result two.txt')))){
    write.table(mrpp_result_two, file = paste0('MRPP result/',prefix,'MRPP result two.txt'), row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }

  return(result)
}



