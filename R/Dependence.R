p_list = c('vegan',
           'dplyr',
           'agricolae',
           'multcomp',
           'multcompView',
           'reshape2',
           'picante',
           'DESeq2',
           'edgeR',
           'Hmisc',
           'psych',
           'igraph',
           'WGCNA',
           'multtest')
for (p in p_list) {
  if (!requireNamespace(p)) {
    install.packages(p)
  }
}
