## @Ni, X., Geng, B., Zheng, H., Shi, J., Hu, G., & Gao, J. (2021). Accurate estimation of single-cell differentiation potency based on network topology and gene ontology information. IEEE/ACM transactions on computational biology and bioinformatics.

##############################################################################

##  description: Processes human datasets including quantile normalization, log2-transformation and other steps, and performs a gene ID conversion on PPI network.

##  usage: Process_Hs(counts, PPI)

##  arguments: 
##  counts: The scRNA-seq data matrix with rows labeling genes and columns labeling single cells.
##  PPI: The adjacency matrix of a user-given PPI network with rownames and colnames labeling a gene ID.

##  value:
##  exp: The scRNA-seq data matrix after preprocessing.
##  adj: The adjacency matrix of PPI network after gene ID conversion.

##############################################################################


Process_Hs<-
  function(counts,PPI){
    
    require("preprocessCore")
    require("AnnotationDbi")
    require("org.Hs.eg.db")
    
    counts <- counts[!duplicated(counts[,1]),]
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    
    if(max(counts) > 100){
      ncounts <- normalize.quantiles(as.matrix(counts),copy=FALSE)
      ncounts[ncounts < 1] <- 1
      ncounts <- log2(ncounts+0.1)
      rownames(ncounts) <- rownames(counts)
      colnames(ncounts) <- colnames(counts)
    }
    
    #Gene ID Conversion
    r <- rownames(PPI)
    geneIDselect <- select(org.Hs.eg.db, keys=r, columns="SYMBOL", keytype="ENTREZID")
    #geneIDselect <- select(org.Hs.eg.db, keys=r, columns="ENSEMBL", keytype="ENTREZID")
    geneIDselect <- geneIDselect[!duplicated(geneIDselect[,1]),]
    rownames(PPI) <- geneIDselect[,2]
    colnames(PPI) <- geneIDselect[,2]
    
    return(list(exp=ncounts,adj=PPI))
  }


