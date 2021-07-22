#' A great way to perform multi-group Wilcox Test
#'
#' @param data Whole data.frame
#' @param group Column containing grouping information
#' @param compare Column containing the difference analysis group
#' @param value Column containing the value
#' @param sort A way to sort the data for add letters with two modes "mean" and "median"
#'
#' @return A list of results
#' @export
#'
#' @examples One_Wilcox(data = ARGs,group = 'gene',compare = 'location',value = 'value',sort = 'mean')
#' @examples One_Wilcox(data = ARGs,group = 'gene',compare = 'location',value = 'value',sort = 'median')
#'
One_Wilcox <- function(data = ARGs,group = 'gene',compare = 'location',value = 'value',sort = 'mean'){
  require(dplyr)
  colnm <- colnames(data)

  if(!is.data.frame(data)){
    message('Error in input: "data" must be format of data.frame ')
  } else if(!(group %in% colnm & compare %in% colnm & value %in% colnm )){
    message('Error in input: "group", "compare", "value" must be one of the cols of data')
  } else if(!(is.numeric(data[,value]))){
    message('Error in input: "value" must be numeric')
  } else if(!(sort %in% c('mean','median'))){
    message('Error in input: "sort" must be mean or median')
  } else{
    aa <- data.frame(stringsAsFactors = F)
    bb <- data.frame(stringsAsFactors = F)
    type <- unique(data[,group])
    for (ele in type){

      sub_dat <- data[data[,group]==ele,]

      ## Paired comparison
      GP <- unique(sub_dat[[compare]]) %>% as.character()

      group1 <- NULL
      group2 <- NULL
      p <- NULL
      med <- data.frame(stringsAsFactors = F)
      mean <- data.frame(stringsAsFactors = F)
      std <- data.frame(stringsAsFactors = F)
      for (i in 1:(length(GP) - 1)) {
        for (j in (i + 1):length(GP)) {
          group1 <- c(group1, GP[i])
          group2 <- c(group2, GP[j])
          group_ij <- subset(sub_dat, sub_dat[[compare]] %in% c(GP[i], GP[j]))
          pos <- which(colnames(group_ij) %in% compare)
          group_ij[,pos] <- factor(group_ij[,pos], levels = c(GP[i], GP[j]))

          wilcox_test <- wilcox.test(group_ij[[value]]~group_ij[[compare]],  alternative = 'two.sided', conf.level = 0.95)
          p <- c(p, wilcox_test$p.value)

          med <- rbind(med,aggregate(group_ij[[value]],by=list(group_ij[[compare]]),FUN=median))
          mean <- rbind(mean,aggregate(group_ij[[value]],by=list(group_ij[[compare]]),FUN=mean))
          std <- rbind(std,aggregate(group_ij[[value]],by=list(group_ij[[compare]]),FUN=sd))
        }
      }
      ## Sorting result
      result <- data.frame(group1, group2,  p)
      result$padj <- p.adjust(result$p, method = 'BH')   #推荐加上 p 值校正，这里使用 Benjamini 方法校正 p 值
      result$sig <- ifelse(result$p <= 0.01,'**',ifelse(result$p <= 0.05,'*', 'no sig'))
      result$sig.adj <- ifelse(result$padj <= 0.01,'**',ifelse(result$padj <= 0.05,'*', 'no sig'))

      result$type <- ele

      ## Add mean and median
      mean <- distinct(mean);rownames(mean) <- mean$Group.1
      med <- distinct(med);rownames(med) <- med$Group.1
      std <- distinct(std);rownames(std) <- std$Group.1

      result$mean1 <- mean[result$group1,][,2]
      result$mean2 <- mean[result$group2,][,2]
      result$median1 <- med[result$group1,][,2]
      result$median2 <- med[result$group2,][,2]
      result$std1 <- std[result$group1,][,2]
      result$std2 <- std[result$group2,][,2]

      ### Add letters
      m = 1; n = 2; l = 1;pp=0.05
      df <- data.frame(trt=c(result$group1,result$group2),
                       mean=c(result$mean1,result$mean2),
                       median=c(result$median1,result$median2),
                       std=c(result$std1,result$std2)) %>% distinct()
      df <- df[order(df[,sort],decreasing = T),]

      trt=unique(df$trt)

      stat_p <- matrix(rep(1,length(trt)*length(trt)),ncol = length(trt)) %>% data.frame(row.names = trt)
      names(stat_p) <- trt
      for (a in trt){
        for (b in trt){
          if(sum(result$group1==a&result$group2==b)){
            stat_p[a,b] <- result$padj[result$group1==a&result$group2==b]
          } else if(sum(result$group1==b&result$group2==a)){
            stat_p[a,b] <- result$padj[result$group1==b&result$group2==a]
          }
        }
      }

      df[,'sig'] <- letters[l]

      if (n < length(trt)){
        while (n < length(trt)) {
          for (n in (m+1):length(trt)) {
            if (stat_p[trt[m],trt[n]] < pp) {
              df[n,'sig'] <- letters[l+1]
              if (n - m != 1) {
                for (x in (m+1):(n-1)) {
                  if (stat_p[trt[x],trt[n]] < pp) df[x,'sig'] <- letters[l]
                  else df[x,'sig'] <- paste(letters[l], letters[l+1], sep = '')
                }
              }
              l = l + 1
              m = n
              break
            }
          }
        }
        df[(m:n),'sig'] <- letters[l]
      }
      ## Pairwise comparisons only have two treatments
      df[2,'sig'] <- ifelse(stat_p[1,2]<0.05,'b','a')

      df[,group] <- ele
      names(df)[1] <- compare
      ## Merge the results
      aa <- rbind(aa,result)
      bb <- rbind(bb,df)
    }
    names(aa) <- c('compare1','compare2','p','padj','sig','sigadj',group,'mean1','mean2','median1','median2','std1','std2')
    re <- list(aa,bb)
    return(re)
  }
}
