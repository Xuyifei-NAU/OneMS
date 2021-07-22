#' One step for multi-group Tukey HSD Test
#'
#'A great way to test the significance if your data possess multi-groups
#'
#' @param data Whole data.frame
#' @param group Column containing grouping information
#' @param compare Column containing the difference analysis group
#' @param value Column containing the value
#'
#' @return A data.frame containing mean, std, results of Tukey HSD Test
#' @export
#'
#' @examples data(ARGs)
#' @examples One_Tukey_HSD1(ARGs,group='gene',compare = 'location',value = 'value')
#' @examples One_Tukey_HSD2(ARGs,group='gene',compare = 'location',value = 'value')
#'
#'
One_Tukey_HSD1 <- function(data=ARGs,group='gene',compare='trt',value='value'){

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
      g1=compare
      sub_dat <- data[data[,group]==i,]

      ## Rename it for later use
      names(sub_dat)[names(sub_dat)==compare] <- 'g1'
      names(sub_dat)[names(sub_dat)==value] <- 'value'
      sub_dat$g1 <- factor(sub_dat$g1)

      fit <- aov(value ~ g1,data = sub_dat )
      Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
      options(warn = -1)
      tuk <- multcomp::cld(multcomp::glht(fit, alternative = 'two.sided', linfct = multcomp::mcp(g1 = 'Tukey')), decreasing = TRUE)
      Tukey.labels <- data.frame(Letters=tuk$mcletters$Letters, stringsAsFactors = FALSE)
      ## Extract the alphabetic grouping row name as the group name
      Tukey.labels$compare = rownames(Tukey.labels)
      Tukey.labels$type <- i

      mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                       aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
      )
      names(mean_sd) <- c('compare','std','mean')

      a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))
    }

    names(a) <- c(compare,'std','mean','Letters',group)
    return(a)
  }
}




One_Tukey_HSD2 <- function(data=ARGs,group='gene',compare='trt',value='value'){

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
      g1=compare
      sub_dat <- data[data[,group]==i,]

      ## Rename it for later use
      names(sub_dat)[names(sub_dat)==compare] <- 'g1'
      names(sub_dat)[names(sub_dat)==value] <- 'value'
      sub_dat$g1 <- factor(sub_dat$g1)

      fit <- aov(value ~ g1,data = sub_dat )
      Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
      options(warn = -1)
      tuk <- multcompView::multcompLetters2(value ~ g1, Tukey_HSD$g1[,"p adj"], sub_dat)

      Tukey.labels <- data.frame(tuk['Letters'], stringsAsFactors = FALSE)
      ## Extract the alphabetic grouping row name as the group name
      Tukey.labels$compare = rownames(Tukey.labels)
      Tukey.labels$type <- i

      mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                       aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
      )
      names(mean_sd) <- c('compare','std','mean')

      a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))
    }

    names(a) <- c(compare,'std','mean','Letters',group)
    return(a)
  }
}


