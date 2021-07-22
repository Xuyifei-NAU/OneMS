#' An easy way to transform data for Veen and Upset plot
#'
#' Transform OTUs table according to the design table
#'
#' @param otu The OTUs table, which row are OTUs and column are samples
#' @param design The design file, a data.frame, should have column of sampleid and group information
#' @param group The column name of the group
#'
#' @return A data.frame
#' @export
#'
#' @examples otu <- Protist$POTUs
#' @examples design <- Protist$PDesign
#' @examples group='trt'
#' @examples One_Format2Veen(otu,design,'trt')
One_Format2Veen <- function(otu,design,group){

  if(!(group %in% colnames(design))){
    stop('"group" must be one of the cols or rownames of metadata')
  }

  idx <-  colnames(otu) %in%  rownames(design)
  otu <- otu[,idx]
  design <- design[colnames(otu),]

  subdesign <- data.frame(row.names = rownames(design),group=design[[group]])

  ### Merge the OTU table with the design table
  data_merge <- merge(subdesign,t(otu),by='row.names')[,-1]

  ### Calculate the mean by group
  data_mean <- aggregate(data_merge[,-1],by=list(data_merge[,1]),mean)
  ### Transform format
  data_long <- as.data.frame(do.call(rbind,data_mean))
  colnames(data_long) <-data_long[1,]
  data_long <- data_long[-1,]
  ### Convert a non-zero number to 1
  data_long[data_long>0] <- 1

  data_long <- as.data.frame(lapply(data_long, as.numeric),row.names = rownames(data_long))
  ## Pick the ID of the Core OTUs
  # core_id <- rownames(data_long)[rowSums(data_long)==ncol(data_long)]
  return(data_long)
}
