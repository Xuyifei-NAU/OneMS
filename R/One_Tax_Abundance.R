#' A convenient way to integrate OTUs tables
#'
#' You can generate counts table or relative abundance table on other ranks, not just OTUs levels
#'
#' @param tax TAXs table, which row are OTUs and column are ranks
#' @param otu OTUs table, which row are OTUs and column are samples
#' @param tar A column name of the TAXs table
#' @param all Logical, if all ranks would be caculated
#' @param norm Logical, if the counts need to be normlized
#'
#' @return A data.frame or list
#' @export
#'
#' @examples tax <- Protist$PTAXs
#' @examples otu <- Protist$POTUs
#' @examples One_Tax_Abundance(tax=tax,otu=otu,tar='phylum',norm=TRUE) ## specified target rank with normalization
#' @examples One_Tax_Abundance(tax=tax,otu=otu,all=TRUE)[['phylum']] ## All ranks
#' @examples One_Tax_Abundance(tax=tax,otu=otu,all=TRUE,norm = TRUE )[['phylum']] ## All ranks with normalization
#'
#'
One_Tax_Abundance <- function(tax,otu,tar=NULL,all=F,norm=F){

  require(dplyr)
  colnm <- colnames(tax)

  if(!(nrow(tax)==nrow(otu))){
    stop("OTUs table must have same rows with TAXs table!")
  } else if(!(is.null(tar))){
    if (!(tar %in% colnm)){
      stop("tar must be one of the column names of TAXs table!")
    }
  }

  if(!(is.logical(all) & is.logical(norm))){
    stop("all and norm must be logical!")
  }

  all_tab <- merge(tax,otu,by='row.names') ## Merge OTU and TAX
  start <- dim(tax)[2]+2 # OTUs table start position
  end <- dim(all_tab)[2] # End position of the OTUs table

  if (all==FALSE){
    if (is.null(tar)) tar <- names(tax)[2]

    select_tab <- data.frame(all_tab[[tar]],all_tab[,start:end]) ## Merge data set
    na.omit(select_tab) # Deletes the row on which the NA value resides
    select_tab[,1] <- factor(select_tab[,1])
    res <- aggregate(select_tab[,-1],by=list(select_tab[,1]),sum)
    names(res)[1] <- tar # Rename

    rownames(res) <- res[,1]
    res <- as.data.frame(res)[,-1]

    if (norm==TRUE){
      res <- t(t(res)/colSums(res)) %>% as.data.frame()
    }

    return(res)
  } else {
    class <- names(tax)
    res_list <- list()
    for (tar in class){
      select_tab <- data.frame(all_tab[[tar]],all_tab[,start:end]) ## Merge data set
      na.omit(select_tab) # Deletes the row on which the NA value resides
      select_tab[,1] <- factor(select_tab[,1])
      res <- aggregate(select_tab[,-1],by=list(select_tab[,1]),sum)
      names(res)[1] <- tar # Rename

      rownames(res) <- res[,1]
      res <- as.data.frame(res)[,-1]

      if (norm==TRUE){
        res <- t(t(res)/colSums(res)) %>% as.data.frame()
      }

      res_list[[tar]] <- res
    }
    return(res_list)
  }
}






#' Integrating abundance according to the group information
#'
#'You can generate counts table or relative abundance table according to the treatments
#'
#' @param otu OTUs table, which row are OTUs(or other ranks) and column are samples
#' @param metadata The metadata file, a data.frame, should have column of sampleid and group information
#' @param id The sample id, which could also be found in OTUs table
#' @param group The column name of the group
#' @param norm Logical, if the counts need to be normlized
#'
#' @return A data.frame with 3 columns
#' @export
#'
#' @examples otu <- Protist$POTUs
#' @examples design <- Protist$PDesign
#' @examples One_Group_Abundance(otu,metadata=design,id='row.names',group = 'temperature',norm = FALSE)
#' @examples One_Group_Abundance(otu,metadata=design,id='row.names',group = 'temperature',norm = TRUE)
#' @examples tax <- Protist$PTAXs
#' @examples res_phy <- One_Tax_Abundance(tax=tax,otu=otu,tar='phylum',norm=TRUE)
#' @examples One_Group_Abundance(otu=res_phy,metadata=design,id='row.names',group = 'temperature',norm = FALSE)
#'
#'
One_Group_Abundance <- function(otu,metadata,id='row.names',group=NULL,norm=T){
  require(dplyr)

  if(!(is.matrix(otu) | is.data.frame(otu))){
    stop('otu table must be matrix or data.frame')
  } else{
    otu <- as.data.frame(otu)
  }

  if((id %in% colnames(metadata) | id=='row.names')& group %in% colnames(metadata) ){
    ## Unified group Name
    pos <- which(names(metadata)==group)
    names(metadata)[pos] <- 'temp_name'

    if (id=='row.names'){
      test_id <- rownames(metadata)
      } else {
      test_id <- metadata[[id]]
      }
    } else{
    stop("'id' and 'group' must one the column name of metadata table")
  }

  if(!(identical(colnames(otu),test_id))){
    stop('OTU table and metadata table must have same sampleid')
  }

  if (norm==TRUE){
    otu <-  t(t(otu)/colSums(otu)) %>% as.data.frame()
  }

  ## Change the OTUs table format
  otu$OTUID <- rownames(otu)
  otu_melt <- reshape2::melt(otu,id.vars='OTUID' ,value.name = 'abundance',variable.name='sampleid')
  data_merge <- merge(otu_melt,metadata,by.x='sampleid',by.y = id,all = T)

  ## Abundance is calculated according to group
  data_merge$OTUID <- factor(data_merge$OTUID)
  data_merge$temp_name <- factor(data_merge$temp_name)
  data <- data_merge %>%
    group_by(OTUID,temp_name) %>%
    mutate(mean_abundance=mean(abundance)) %>% as.data.frame()

  ## Sorting result
  res <- data.frame(OTUID=data[["OTUID"]],group=data[["temp_name"]],mean_abundance=data[["mean_abundance"]])%>% distinct()
  names(res)[2] <- group
  return(res)
  }




#' A method for screening OTU tables according to abundance
#'
#' OTUs tables can be filtered based on abundance thresholds or number thresholds
#'
#' @param otu OTUs table, which row are OTUs(or other ranks) and column are samples
#' @param input_format The value of the OTUs table, counts or relative (abundance)
#' @param mt The threshold of mean relative abundance
#' @param num The numbers of OTUs you want to have in the result
#'
#' @return A data.frame with relative abundance
#' @export
#'
#' @examples otu <- Protist$POTUs
#' @examples One_Filter_Aabundance(otu,mt=0.001)
#' @examples One_Filter_Aabundance(otu,num=200)
#' @examples One_Filter_Aabundance(otu,mt=0.01,num = 3)
#' @examples One_Filter_Aabundance(otu,mt=0.01,num = 40)
#' @examples ## Filter Abundance table with other ranks
#' @examples tax <- Protist$PTAXs
#' @examples res_phy <- One_Tax_Abundance(tax=tax,otu=otu,tar='phylum',norm=TRUE)
#' @examples One_Filter_Aabundance(res_phy,mt=0.01,num = 20)
#' @examples
#'
One_Filter_Aabundance <- function(otu,input_format='count',mt=NULL,num=NULL){

  require(dplyr)

  if(is.null(mt) & is.null(num)){
    stop('One of the mt or num must be assigned!')
  }

  ## Caculate the relative abundance
  if (input_format=='count'){
    re_ab <-  t(t(otu)/colSums(otu))
  } else if (input_format=='relative'){
    re_ab <- otu
  } else {
    stop('input_format must be set as "count" or "relative"')
  }

  ab_all <- data.frame(abundance=rowMeans(re_ab),ID=rownames(re_ab))
  ab_all <- ab_all [order(ab_all$abundance,decreasing = T),]

  ## The extraction is performed according to the abundance threshold by default
  if (!(is.null(mt))){
    ## Extraction is performed according to abundance
    filter_ID <- ab_all[ab_all$abundance>=mt,][,2]
    ## Consider the setting of num
    if(!is.null(num)){
      len <- length(filter_ID)
      len2 <- ifelse(num<len,num,len)
      filter_ID <- filter_ID[1:len2]
    }

    } else {
    ## Extract according to the number
    filter_ID <-ab_all [order(ab_all$abundance,decreasing = T),][1:num,2]
  }

  ## Extract the table by ID
  ID <- rownames(otu) %in% filter_ID
  final_table <- re_ab[ID,] %>% as.data.frame()
  return(final_table)
}












