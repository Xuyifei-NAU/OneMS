etwd('~/Test_Data/18S/lefse/')
# 1. 读取OTU表
otutab = read.table('~/Test_Data/18S/lefse/otu_table.txt', header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# 2. 读取物种注释
taxonomy = read.table("~/Test_Data/18S/lefse/taxonomy.txt", header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F)
# 3. 读取样本分组信息
metadata = read.table("~/Test_Data/18S/lefse/metadata.txt", header=T, row.names=1, sep="\t", comment.char="")

##
library(do)

## 方式1 把 _unclassified 的字段删除 ----
# taxonomy2 <- Replace(taxonomy,from='_unclassified',to='')
# for (i in colnames(taxonomy2)){
#   taxonomy2[[i]] <-  paste0(tolower(strsplit(i,"")[[1]][1]),'_',taxonomy2[[i]]) ## 添加物种分类的前缀
# }

## 方式2 把含有 _unclassified 的单元格删除 ----
taxonomy3 <- taxonomy
for (i in colnames(taxonomy3)){
  taxonomy3[[i]] <-  paste0(tolower(strsplit(i,"")[[1]][1]),'_',taxonomy3[[i]])
}

taxonomy3 <- Replace0(taxonomy3,from=c('\\S+_unclassified','')) ## 将_unclassified的单元格全部替换为“”
taxonomy3 <- taxonomy3[!apply(taxonomy3 == "", 1, all),] # 将全空的行删除

colSums(taxonomy3=='')/nrow(taxonomy3) ## 计算“”的占比
# Kingdom     Phylum      Class      Order     Family      Genus    Species
# 0.00000000 0.04601006 0.09273904 0.12149533 0.22501797 0.46441409 0.56146657

source('~/MyScript/One_format2lefse.R')
# One_format2lefse(otutab, taxonomy3, metadata, thre = 0.01, groupID='temperature',sub_group = 'precipitation',rank_start = 'Kingdom',rank = 'Family', output='temp/LEfSe_family.txt')



One_format2lefse <- function (otutab, taxonomy, metadata, thre = 0.01, groupID = "Group", sub_group='None',rank_start='None',rank='None',
                              output = "LEfSe.txt")
{
  idx = rownames(otutab) %in% rownames(taxonomy)
  otutab = otutab[idx, ]
  tax = taxonomy[rownames(otutab), ]


  p_list = c("dplyr",'tidyverse','stringr')
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
  norm = t(t(otutab)/colSums(otutab, na = T)) * 100
  idx = rowMeans(norm) > thre
  HA = norm[idx, ]
  tax = tax[rownames(HA), ]
  tax$Phylum = paste(tax$Kingdom, tax$Phylum, sep = "|")
  tax$Class = paste(tax$Phylum, tax$Class, sep = "|")
  tax$Order = paste(tax$Class, tax$Order, sep = "|")
  tax$Family = paste(tax$Order, tax$Family, sep = "|")
  tax$Genus = paste(tax$Family, tax$Genus, sep = "|")
  tax$Species = paste(tax$Genus, tax$Species, sep = "|")
  grp <- tax[rownames(tax), "Kingdom", drop = F]
  merge = cbind(HA, grp)
  HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
  colnames(HA_Kingdom)[1] = "Class"
  grp <- tax[rownames(tax), "Phylum", drop = F]
  merge = cbind(HA, grp)
  HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
  colnames(HA_Phylum)[1] = "Class"
  grp <- tax[rownames(tax), "Class", drop = F]
  merge = cbind(HA, grp)
  HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
  colnames(HA_Class)[1] = "Class"
  grp <- tax[rownames(tax), "Order", drop = F]
  merge = cbind(HA, grp)
  HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
  colnames(HA_Order)[1] = "Class"
  grp <- tax[rownames(tax), "Family", drop = F]
  merge = cbind(HA, grp)
  HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
  colnames(HA_Family)[1] = "Class"
  grp <- tax[rownames(tax), "Genus", drop = F]
  merge = cbind(HA, grp)
  HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
  colnames(HA_Genus)[1] = "Class"
  all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family,
              HA_Genus)
  metadata$group = metadata[, groupID]

  ## 针对 rank 进行数据整理
  if (rank != 'None' & (rank %in% names(taxonomy)) & rank_start != 'None' & (rank_start %in% names(taxonomy))){
    class_id <- which(names(taxonomy)==rank) ## 获取rank的位置
    star_pos <- which(names(taxonomy)==rank_start) ## 获取rank起始位置

    all <- all[str_count(all$Class,"\\|") ==(class_id-1),] ## 按照rank位置拆分分类表

    all_split <- separate(all,col = 'Class',sep='\\|',into = names(taxonomy))[,star_pos:class_id] # 按照 ｜ 拆分列

    filter_id <- all_split[[rank]] != "" ## 获取最细分类级水平下，注释到的行的索引

    all_unite <- unite(all_split[filter_id,], "Class", sep = "|", remove = T)

    all <- cbind(all_unite,all[,-1][filter_id, ])

  }

  if (rank != 'None' & (rank %in% names(taxonomy)) & rank_start == 'None'){
    class_id <- which(names(taxonomy)==rank) ## 获取rank的位置

    all <- all[str_count(all$Class,"\\|") ==(class_id-1),] ## 按照rank位置拆分分类表

    all_split <- separate(all,col = 'Class',sep='\\|',into = names(taxonomy))[,class_id] # 按照 ｜ 拆分列

    all <- cbind(all_split,all[,-1])[all_split != "", ]

    names(all)[1] <- 'Class'
  }



  if(sub_group == 'None'){
    colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]],
    ]$group)
  }


  if(sub_group!='None' & (sub_group %in% names(metadata)) ){
    colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]],
    ]$group)
    all <- rbind(sub_class=c('sub_class',as.character(metadata[[sub_group]])),all)
  }

  write.table(all, file = paste(output, sep = ""), append = FALSE,
              sep = "\t", quote = F, row.names = F, col.names = T)

  return(all)
}


