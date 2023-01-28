## @Ni, X., Geng, B., Zheng, H., Shi, J., Hu, G., & Gao, J. (2021). Accurate estimation of single-cell differentiation potency based on network topology and gene ontology information. IEEE/ACM transactions on computational biology and bioinformatics.

##############################################################################

##  description: Processes mouse datasets including quantile normalization, log2-transformation and other steps.

##  usage: Process_Mm(counts)

##  arguments: 
##  counts: The scRNA-seq data matrix with rows labeling genes and columns labeling single cells.

##  value:
##  exp: The scRNA-seq data matrix after preprocessing.

##############################################################################


Process_Mm<-
  function(counts){
    
    require("preprocessCore")
    require("AnnotationDbi")
    require("org.Mm.eg.db")
    
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
       
      return(ncounts)
  }


