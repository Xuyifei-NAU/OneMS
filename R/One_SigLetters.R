#' Add letters according to the p-values
#'
#'A good way to add letters to results of  multiple comparisons
#'
#' @param df A data.frame with group and mean/median information
#' @param id The column of group information
#' @param sort  The column of mean/median values
#' @param pvalue  The p-value matrix, both n*n matrix and n*3 matrix can be accepted
#'
#' @return A data.frame with significance letter result
#' @export
#'
#' @examples df <- data.frame(id=c('A','B','C','D'),mean=c(0.6,0.4,0.7,0.1))
#' @examples pvalue <- data.frame(cp1 = c('A','A','A','B','B','C'),cp2=c('B','C','D','C','D','D'),p =c('0.09','0.1','0.04','0.03','0.045','0.01'))
#' @examples One_SigLetters(df,id="id",sort='mean',pvalue = pvalue)
#'
One_SigLetters <- function(df,id,sort=NULL,pvalue){
  require(dplyr)
  colnm <- colnames(df)

  if(!is.data.frame(df)){
    message('Error in input: "df" must be format of data.frame ')
  } else if(!(id %in% colnm)){
    message('Error in input: "id" must be one of the cols of df')
  } else if(!((sum(rownames(pvalue) %in% colnames(pvalue))==nrow(pvalue)) | ncol(pvalue)==3)){
    message('Error in input: "pvalue" must be n*n or n*3')
  } else {

    ## if sort was assigned, the df would be sorted according to the value
    if (!is.null(sort)){
      if ( sort %in% colnm){
        df <- df[order(df[,sort],decreasing = T),]
      } else {
        message("Warning: 'sort' was not performed, for the col 'sort' couldn't be found")
      }
    }

    compare=df[,id]

    ## transform 3 column p-value matrix
    if(sum(rownames(pvalue) %in% colnames(pvalue))==nrow(pvalue)){
      stat_p <- pvalue
    } else {
      pvalue[,3] <- as.numeric(pvalue[,3])
      stat_p <- matrix(rep(1,length(compare)*length(compare)),ncol = length(compare)) %>% data.frame(row.names = compare)
      names(stat_p) <- compare
      for (a in compare){
        for (b in compare){
          if(sum(pvalue[,1]==a & pvalue[,2]==b)){
            stat_p[a,b] <- pvalue[,3][pvalue[,1]==a & pvalue[,2]==b]
          } else if(sum(pvalue[,1]==b & pvalue[,2]==a)){
            stat_p[a,b] <- pvalue[,3][pvalue[,1]==b & pvalue[,2]==a]
          }
        }
      }
    }

    ## Initialize the result data.frame
    m = 1; n = 2; l = 1;pp=0.05
    df[,'sig'] <- letters[l]
    ## 比较p值，标字母
    if (n < length(compare)){
      while (n < length(compare)) {
        for (n in (m+1):length(compare)) {
          if (stat_p[compare[m],compare[n]] < pp) {
            df[n,'sig'] <- letters[l+1]
            if (n - m != 1) {
              for (x in (m+1):(n-1)) {
                if (stat_p[compare[x],compare[n]] < pp) df[x,'sig'] <- letters[l]
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

    ## result
    return(df)
  }

}
