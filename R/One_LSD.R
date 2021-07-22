#' One step for multi-group LSD Test
#'
#'A great way to test the significance if your data possess multi-groups
#'
#' @param data Whole data.frame
#' @param group Column containing grouping information
#' @param compare Column containing the difference analysis group
#' @param value Column containing the value
#'
#' @return A data.frame containing mean, std, results of LSD
#' @export
#'
#' @examples data(ARGs)
#' @examples ONE_LSD(data=ARGs,group='gene',compare='trt',value='value')
#'
ONE_LSD <- function(data=ARGs,group='gene',compare='trt',value='value'){

  colnm <- colnames(data)

  if(!is.data.frame(data)){
    message('Error in input: "data" must be format of data.frame ')
  } else if(!(group %in% colnm & compare %in% colnm & value %in% colnm )){
    message('Error in input: "group", "compare", "value" must be one of the cols of data')
  } else if(!(is.numeric(data[,value]))){
    message('Error in input: "value" must be numeric')
  } else{
    a <- data.frame(stringsAsFactors = F)
    type <- unique(data[,group])
    for (i in type)
    {
      # sub_dat <- subset(data,group == i)
      sub_dat <- data[data[,group]==i,]
      # fit <- aov(value ~ compare,sub_dat)
      fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
      out <- agricolae::LSD.test(fit,'sub_dat[, compare]',p.adj='BH')

      out$groups$type <- i
      out$groups$compare <- rownames(out$groups)

      a <- rbind(a,merge(out$means[,1:2], out$groups,by='sub_dat[, value]'))
      # print(a)
    }
    names(a) <- c('mean','std','lsd',group,compare)
    return(a)
  }
}
