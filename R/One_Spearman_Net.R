#' One step to construct Spearman Network
#'
#' A easy way to construct Spearman Network and caculate the network features altogether
#'
#' @param normtab Standardized OTU table, which row are OTUs and column are samples
#' @param tax TAXs table, which row are OTUs and column are ranks
#' @param rank Which rank of the normtab, usually "OTU"
#' @param abundance Abundance threshold, numeric
#' @param appear The least occurrence of OTUs in the sample, numeric
#' @param num The numbers of OTUs you want to keep to construct the network, numeric
#' @param cor_type Types of correlation, could be set as 'pearson','kendall','spearman'
#' @param P.adj  Method for p adjust, could be set as 'bonferroni','holm','hochberg','hommel','fdr','BH','BY','none'
#' @param rt The threshold of correlation, numeric
#' @param pt The threshold of p-value, numeric
#' @param write Logical, whether to write result to files for Gephi or Cytoscape
#' @param prefix  Prefix for files
#'
#' @return A list of igraph object, node table, edge table, and correlation matrix
#' @export
#'
#' @examples normtab <- Protist$POTUs
#' @examples normtab <- as.data.frame(t(t(normtab)/colSums(normtab)))
#' @examples tax <- Protist$PTAXs
#' @examples One_Spearman_Net(normtab=normtab,tax=tax,rank='OTU',abundance = 0.001,appear=4,num=NULL,cor_type = 'spearman',P.adj = 'BH',rt=0.7,pt=0.05,write = TRUE)

One_Spearman_Net <- function(normtab,tax,rank='OTU',abundance = 0,appear=0,num=NULL,cor_type = 'spearman',
                        P.adj = 'BH',rt=0.7,pt=0.05,write=F,prefix=NULL){
  ## Packages
  require(Hmisc)
  require(psych)
  require(igraph)
  require(dplyr)
  require(WGCNA)
  require(multtest)

  if(!(is.numeric(abundance)&is.numeric(appear)&is.numeric(rt)&is.numeric(pt))){
    stop('"abundance,appear,rt and pt must be numeric')
  }

  if( rank!='OTU' & !(rank %in% colnames(tax)) ){
    stop('rank must be set as "OTU" or other classification rank ')
  }

  if(!(cor_type %in%  c('pearson','kendall','spearman') )){
    stop("method for cor_type not found")
  }

  p.adj_list <- c('bonferroni','holm','hochberg','hommel','fdr','BH','BY','none' )
  if(!(P.adj %in% p.adj_list)){
    stop("P.adj method not found")
  }



  idx <- rownames(normtab) %in% rownames(tax)
  normtab <- normtab[idx,]
  tax <- tax[rownames(normtab),]

  ## 1 Filter otutabs -----
  filtered_otu <- normtab[which(rowMeans(normtab) >= abundance), ]

  if (!(is.null(appear))&is.numeric(appear)){
    filtered_otu1 <- filtered_otu
    filtered_otu1[filtered_otu1>0] <- 1
    filtered_otu <- filtered_otu[which(rowSums(filtered_otu1) >= appear), ]
  }

  if (!(is.null(num))&is.numeric(num)) {
    filtered_otu$mean <- rowMeans(filtered_otu)
    filtered_otu <- filtered_otu[order(filtered_otu$mean,decreasing = T),]
    filtered_otu <- filtered_otu[1:num,!names(filtered_otu) %in% c("mean")]
    filtered_otu <- filtered_otu[,!names(filtered_otu) %in% c("mean")] # 把生成的mean删除
  }


  # 2 Calculate the correlation value and pvalue ------
  occor<-WGCNA::corAndPvalue(t(filtered_otu),method = c( "spearman"))
  # multiple test the p values
  mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc=P.adj)
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(t(filtered_otu))[2])
  ## R value
  occor.r<-occor$cor
  diag(occor.r) <- 0
  # Determine the threshold value of the interaction relationship between species, and convert the incongruent data in the correlation R matrix to 0
  occor.r[occor.p>pt|abs(occor.r)<rt] = 0
  occor.r[is.na(occor.r)]=0


  # 3 To build the network --------------------------------------------------------------
  if(sum(occor.r)==0) {
    warning('Unable to form a network')
    res <- NULL
  } else {
    g <- graph.adjacency(as.matrix(occor.r), weighted = TRUE, mode = 'undirected')
    g

    #Autocorrelation removal
    g <- simplify(g)
    #Delete the nodes with 0 degree
    g <- delete.vertices(g, names(degree(g)[degree(g) == 0]) )

    # 4 Calculate the properties of the network nodes and edges-------------------------------------------------------------
    E(g)$correlation <- E(g)$weight
    E(g)$weight <- abs(E(g)$weight)

    # negative or positive correlation
    E(g)$state <-  E(g)$correlation
    E(g)$state[E(g)$state<0] <- -1
    E(g)$state[E(g)$state>0] <- 1


    # Add tax annotation information to the nodes
    if (rank=='OTU'){
      filtered_tax <- tax[as.character(V(g)$name), ]
    } else if(rank %in% colnames(tax)){
      pos_id <- which(colnames(tax)==rank)
      delt <- ncol(tax)-pos_id
      filtered_tax <- tax[,1:pos_id] %>% dplyr::distinct() # De-redundancy of columns
      ## Add the UN
      for (i in 1:delt){
        filtered_tax[[colnames(tax)[pos_id+i]]] <- 'UN'
      }
      rownames(filtered_tax) <- filtered_tax[[rank]] # Rename the rows
      filtered_tax <- filtered_tax[as.character(V(g)$name), ] ## Sort
    } else{
      stop('"rank" must one of the column names of the tax table')
    }
    names(filtered_tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

    V(g)$Kingdom <- filtered_tax$Kingdom
    V(g)$Phylum <- filtered_tax$Phylum
    V(g)$Class <- filtered_tax$Class
    V(g)$Order <- filtered_tax$Order
    V(g)$Family <- filtered_tax$Family
    V(g)$Genus <- filtered_tax$Genus

    # Degree
    V(g)$degree <- degree(g)

    # Weighted degree
    V(g)$weight_degree <- strength(g)

    # Closeness centrality (There are two algorithms, but I don't know which one is right.)
    V(g)$closeness_centrality <- closeness(g)
    # V(g)$closeness_centrality2 <- centralization.closeness(g)$res


    # Betweenness centrality
    V(g)$betweenness_centrality <- betweenness(g)

    # Eigenvector centrality
    V(g)$eigenvector_centrality <- evcent(g)$vector

    # Edge betweenness centrality
    E(g)$betweenness_centrality <- edge.betweenness(g)

    # modularity
    fc <- cluster_fast_greedy(g)
    modularity <- modularity(g, membership(fc))

    set.seed(888)
    V(g)$modularity <- membership(cluster_fast_greedy(g))

    # Zi and Pi
    ## function
    zi.pi<-function(nodes_bulk, z.bulk, modularity_class, degree){

      z.bulk[abs(z.bulk)>0]<-1
      module<-which(colnames(nodes_bulk)==modularity_class)
      module.max<-max(nodes_bulk[,module])
      degree<-which(colnames(nodes_bulk)==degree)

      # Divide the related matrix into modules
      bulk.module<-list(NA)
      length(bulk.module)<-module.max

      for(i in 1:max(nodes_bulk[,module])){
        bulk.module[[i]]<-z.bulk[which(nodes_bulk[,module]==i),which(nodes_bulk[,module]==i)]
        bulk.module[[i]]<-as.data.frame(bulk.module[[i]])
        rownames(bulk.module[[i]])<-rownames(z.bulk)[which(nodes_bulk[,module]==i)]
        colnames(bulk.module[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
      }

      # within-module degree z
      z_bulk<-list(NA)
      length(z_bulk)<-module.max

      for(i in 1:length(z_bulk)){
        z_bulk[[i]]<-bulk.module[[i]][,1]
        z_bulk[[i]]<-as.data.frame(z_bulk[[i]])
        colnames(z_bulk[[i]])<-"z"
        rownames(z_bulk[[i]])<-rownames(bulk.module[[i]])
      }

      # Caculate z value
      for(i in 1:max(nodes_bulk[,module])){
        if(length(bulk.module[[i]])==1){
          z_bulk[[i]][,1]<-0
        }else if(sum(bulk.module[[i]])==0){
          z_bulk[[i]][,1]<-0
        }else{
          k<-rowSums(bulk.module[[i]])
          mean<-mean(k)
          sd<-sd(k)
          if (sd==0){
            z_bulk[[i]][,1]<-0
          }else{
            z_bulk[[i]][,1]<-(k-mean)/sd
          }
        }
      }

      for(i in 2:max(nodes_bulk[,module])) {
        z_bulk[[i]]<-rbind(z_bulk[[i-1]],z_bulk[[i]])
      }
      z_bulk<-z_bulk[[module.max]]

      # Split the related matrix columns by module
      bulk.module1<-list(NA)
      length(bulk.module1)<-module.max

      for(i in 1:max(nodes_bulk[,module])){
        bulk.module1[[i]]<-z.bulk[,which(nodes_bulk[,module]==i)]
        bulk.module1[[i]]<-as.data.frame(bulk.module1[[i]])
        rownames(bulk.module1[[i]])<-rownames(z.bulk)
        colnames(bulk.module1[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
      }

      # among-module connectivity c
      c_bulk<-list(NA)
      length(c_bulk)<-module.max

      for(i in 1:length(c_bulk)){
        c_bulk[[i]]<-z.bulk[,1]
        c_bulk[[i]]<-as.matrix(c_bulk[[i]])
        colnames(c_bulk[[i]])<-"c"
        rownames(c_bulk[[i]])<-rownames(z.bulk)
        c_bulk[[i]][,1]<-NA
      }

      # The square of the connections of each module of each node
      for(i in 1:max(nodes_bulk[,module])){
        c_bulk[[i]]<-rowSums(bulk.module1[[i]])
        c_bulk[[i]]<-as.matrix(c_bulk[[i]])
        c_bulk[[i]]<-c_bulk[[i]]*c_bulk[[i]]
        colnames(c_bulk[[i]])<-"c"
        rownames(c_bulk[[i]])<-rownames(z.bulk)
      }

      # Sum of squares
      for(i in 2:max(nodes_bulk[,module])){
        c_bulk[[i]]<-c_bulk[[i]]+c_bulk[[i-1]]
      }
      c_bulk<-c_bulk[[module.max]]

      c_bulk1<-1-(c_bulk/(nodes_bulk[,degree]*nodes_bulk[,degree]))
      colnames(c_bulk1)<-"c"

      # Integration of z and c values
      z_c_bulk<-c_bulk1
      z_c_bulk<-as.data.frame(z_c_bulk)
      z_c_bulk$z<-z_bulk[match(rownames(c_bulk1),rownames(z_bulk)),]
      z_c_bulk<-z_c_bulk[,c(2,1)]
      names(z_c_bulk)[1:2]<-c('within_module_connectivities','among_module_connectivities')

      z_c_bulk$nodes_id<-rownames(z_c_bulk)
      nodes_bulk$nodes_id<-rownames(nodes_bulk)
      z_c_bulk<-merge(z_c_bulk,nodes_bulk,by='nodes_id')
      z_c_bulk

    }
    # A list of node properties containing the modules that the node is divided into
    nodes_list <- data.frame(
      nodes_id = V(g)$name,
      degree = V(g)$degree,
      modularity = V(g)$modularity
    )
    rownames(nodes_list) <-  nodes_list$nodes_id
    # The order of nodes in both files should be the same
    z2 <- occor.r[V(g)$name,V(g)$name]
    nodes_list <- nodes_list[rownames(z2), ]
    # Caculate Zi and Pi
    zi_pi <- zi.pi(nodes_list, z2, degree = 'degree', modularity_class = 'modularity')
    V(g)$within_module_connectivities <- zi_pi$within_module_connectivities
    V(g)$among_module_connectivities <- zi_pi$among_module_connectivities

    ## 4.1 A table that summarizes node characteristics ----
    nodes <- data.frame(
      nodes_id = V(g)$name,
      degree = V(g)$degree,
      modularity = V(g)$modularity,
      weight_degree = V(g)$weight_degree,
      closeness_centrality = V(g)$closeness_centrality,
      betweenness_centrality = V(g)$betweenness_centrality,
      eigenvector_centrality = V(g)$eigenvector_centrality,
      within_module_connectivities = V(g)$within_module_connectivities,
      among_module_connectivities = V(g)$among_module_connectivities,
      Kingdom = V(g)$Kingdom,
      Phylum = V(g)$Phylum,
      Class = V(g)$Class,
      Order = V(g)$Order,
      Family = V(g)$Family,
      Genus = V(g)$Genus
    )
    ## 4.2 A table that summarizes the edge features ----
    edge <- data.frame(as_edgelist(g))
    edges <- data.frame(
      source = edge[[1]],
      target = edge[[2]],
      weight = E(g)$weight,
      correlation = E(g)$correlation,
      state = E(g)$state,
      betweenness_centrality = E(g)$betweenness_centrality
    )


    # 5 Network characteristic calculation ------------------------------------------------------------------

    # Number of nodes and Number of edges
    nodes_num <- length(V(g))

    edges_num <- length(E(g))

    # Average degree
    average_degree <- mean(degree(g))

    # Average weighted degree (Applicable to weighted networks only)
    #average_weight_degree <- mean(strength(igraph))

    # Connectivity
    nodes_connectivity <- vertex.connectivity(g)

    edges_connectivity <- edge.connectivity(g)

    # Average path length
    average_path_length <- average.path.length(g, directed = FALSE)

    # Diameter
    graph_diameter <- diameter(g, directed = FALSE)

    # Density
    graph_density <- graph.density(g)

    # Clustering coefficient
    clustering_coefficient <- transitivity(g)

    # Betweenness centralization
    betweenness_centralization <- centralization.betweenness(g)$centralization

    # Degree centralization
    degree_centralization <- centralization.degree(g)$centralization

    # Modularity
    modularity


    # 5.1 Select the sections to do a summary output ----
    network_character <- data.frame(
      nodes_num,    #（number of nodes）
      edges_num,    #（number of edges）
      average_degree,    #（average degree)
      nodes_connectivity,    #（vertex connectivity）
      edges_connectivity,    #（edges connectivity）
      average_path_length,    #（average path length）
      graph_diameter,    #（diameter）
      graph_density,    #（density）
      clustering_coefficient,    #（clustering coefficient）
      betweenness_centralization,    #（betweenness centralization)
      degree_centralization,    # (degree centralization)
      modularity    #（modularity）
    )

    # 6 A list of aggregated results -------------------------------------------------------------------------

    res <- list(g=g,nodes_list=nodes,edges_list=edges,network_character=network_character,cormat=occor.r)

    if(write==TRUE){
      if(!(dir.exists('Spearman Network'))){
        dir.create('Spearman Network')
      }

      if(!(file.exists(paste0('Spearman Network/','Spearman Network ',prefix,'.gml')))){
        write.graph(g, paste0('Spearman Network/','Spearman Network ',prefix,'.gml'), format = 'gml')
      }
      if(!(file.exists(paste0('Spearman Network/','Spearman Network ',prefix,'.graphml')))){
        write.graph(g, paste0('Spearman Network/','Spearman Network ',prefix,'.graphml'), format = 'graphml')
      }
    }

  }

  return(res)

}
