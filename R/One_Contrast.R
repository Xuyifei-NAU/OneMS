# DESeq2 ------------------------------------------------------------------
#' A convenient way to do one-way DESeq2 analysis
#'
#' For loop was used to conduct the multi-group pairwise DESeq2 analysis
#'
#' @param otu The OTUs table, which row are OTUs and column are samples
#' @param design The metadata file, a data.frame, should have column of sampleid and group information
#' @param treatment The column name of the treatment id
#' @param sig Numeric, the threshold of significant p-value
#' @param padj Method for p adjust
#' @param prefix Prefix for files
#'
#' @return A data.frame and files
#' @export
#'
#' @examples otu <- Protist$POTUs
#' @examples design <- Protist$PDesign
#' @examples One_DESeq2(otu,design = design,treatment = 'precipitation',prefix = 'prec')
#' @examples One_DESeq2(otu,design = design,treatment = 'temperature',prefix = 'temp')
#'
#'
One_DESeq2 <- function(otu,design,treatment,sig=0.05,padj='BH',prefix=NULL)
{
  require(DESeq2)
  p.adj_list <- c('bonferroni','holm','hochberg','hommel','fdr','BH','BY','none' )

  if(!(treatment %in% colnames(design))){
    stop("'treatment' must one the column names of design table")
  }
  if(!(padj %in% p.adj_list)){
    stop("padj method not found")
  }
  if(!(is.numeric(sig))){
    stop("sig must be numeric")
  }



  idx <-  colnames(otu) %in%  rownames(design)
  count <- otu[,idx]
  subdesign <- design[colnames(count),]

  ## reconstruct subdesign table
  subdesign <- data.frame(row.names = rownames(subdesign),group=subdesign[,treatment])

  # Build the DESeq object
  count <- as.matrix(count)
  subdesign <- as.matrix(subdesign)
  dds_group <- DESeqDataSetFromMatrix(countData = count,colData = subdesign,design = ~ group)
  ## DESeq function for analysis
  dds <- DESeq(dds_group)

  ## Get treatment info
  trts =unique(subdesign[,which(colnames(subdesign)=='group')])
  ## Pairwise comparsion
  for (i in 1:(length(trts)-1))
  {
    for (j in (i+1):length(trts))
    {
      if (j != i)
      {
        result <- results(dds,contrast = c('group',trts[i],trts[j]),alpha = sig,pAdjustMethod=padj)
        result = result[order(result$padj),] # Sort by p value

        if(!(dir.exists('DESeq2 result'))){
          dir.create('DESeq2 result')
        }

        if(!file.exists(paste0('DESeq2 result/','DESeq2_',prefix,'_',trts[i],'-',trts[j],".txt"))){
          write.table(data.frame ("OUTID"= rownames(result), result),file=paste0('DESeq2 result/','DESeq2_',prefix,'_',trts[i],'-',trts[j],".txt"),sep="\t",quote=F,row.names=F)
        }

        if (i==1&j==2)
        {
          result$compare <- paste0(trts[i],'-',trts[j])
          aa <- as.data.frame(result)
        }
        else
        {
          result$compare <- paste0(trts[i],'-',trts[j])
          aa <- rbind(aa,as.data.frame(result))
        }
      }
    }
  }
  return(aa)
}




#  edgeR ------------------------------------------------------------------


#' A convenient way to do one-way edgeR analysis
#'
#' For loop was used to conduct the multi-group pairwise edgeR analysis
#'
#' @param otu The OTUs table, which row are OTUs and column are samples
#' @param design The metadata file, a data.frame, should have column of sampleid and group information
#' @param treatment The column name of the treatment id
#' @param sig Numeric, the threshold of significant p-value
#' @param padj Method for p adjust
#' @param prefix Prefix for files
#'
#' @return A data.frame and files
#' @export
#'
#' @examples otu <- Protist$POTUs
#' @examples design <- Protist$PDesign
#' @examples One_EdgeR(otu,design = design,treatment = 'precipitation',prefix = 'prec')
#' @examples One_EdgeR(otu,design = design,treatment = 'temperature',prefix = 'temp')
#'
#'

One_EdgeR <- function(otu,design,treatment,sig=0.05,padj='BH',prefix=NULL)
{
  require(edgeR)
  require(dplyr)

  p.adj_list <- c('bonferroni','holm','hochberg','hommel','fdr','BH','BY','none' )

  if(!(treatment %in% colnames(design))){
    stop("'treatment' must one the column names of design table")
  }
  if(!(padj %in% p.adj_list)){
    stop("padj method not found")
  }
  if(!(is.numeric(sig))){
    stop("sig must be numeric")
  }



  idx <-  colnames(otu) %in%  rownames(design)
  count <- otu[,idx]
  subdesign <- design[colnames(count),]
  ## reconstruct subdesign table
  subdesign <- data.frame(row.names = rownames(subdesign),group=subdesign[,treatment])



  # DGEList
  d_group <- DGEList(counts = count,group = subdesign[,'group']) %>% calcNormFactors()
  # Organize and design the grouping matrix
  design.group <- model.matrix(~0+d_group$samples$group)
  colnames(design.group) <-levels(d_group$samples$group)

  ## Analysis
  fit <- estimateGLMCommonDisp(d_group,design.group) %>% estimateTagwiseDisp() %>% glmFit(design.group)

  design_matrix=design.group
  trts <- colnames(design_matrix)

  for (i in 1:(length(trts)-1))
  {
    for (j in (i+1):length(trts))
    {
      if (j != i)
      {
        contrast_list <- makeContrasts(contrasts = paste0(trts[i],'-',trts[j]),levels = design_matrix)
        LRD <- glmLRT(fit,contrast =contrast_list )
        res <- decideTestsDGE(LRD,p.value = sig,adjust.method = padj,lfc = 0)
        res_final <- as.data.frame(LRD$table)
        res_final <- merge(res_final,as.data.frame(res),by='row.names')
        res_final <- cbind(res_final,padj=p.adjust(res_final$PValue,method = padj))
        res_final$level <- ifelse(res_final[,6]==1,'enriched',ifelse(res_final[,6]==-1,'depleted','nosig'))

        if(!(dir.exists('EdgeR result'))){
          dir.create('EdgeR result')
        }

        if(!file.exists(paste0('EdgeR result/','EdgeR',prefix,'_',trts[i],'-',trts[j],".txt"))){
          write.table(data.frame (res_final[,-6]),file= paste0('EdgeR result/','EdgeR_',prefix,'_',trts[i],'-',trts[j],".txt"),sep="\t",quote=F,row.names=F)
        }

        if (i==1&j==2)
        {
          res_final$compare <- paste0(trts[i],'-',trts[j])
          bb <- as.data.frame(res_final[,-6])
        }
        else
        {
          res_final$compare <- paste0(trts[i],'-',trts[j])
          bb <- rbind(bb,as.data.frame(res_final[,-6]))
        }
      }
    }
  }
  return(bb)
}




