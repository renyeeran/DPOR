##############################################################################

##  description: Gene ID conversion in processing mouse datasets.

##  usage: M2H(mdata)

##  arguments: 
##  mdata: The mouse scRNA-seq data matrix after preprocessing.

##  value:
##  hdata: The human homologene scRNA-seq data matrix.

##############################################################################


M2H<-	
  function(mdata){
    
    library(homologene)
    library(org.Mm.eg.db)
    
    genelist <- mdata[,1]
    
    id <- homologene(genelist, inTax=10090, outTax=9606)
    id <- id[!duplicated(id[,1]),]
    id_list <- id[,"10090"]
    
    idx <- match(id_list, genelist)
    hdata <- mdata[idx,]
    hdata[,1] <- id1[,"9606_ID"]
    
    return(hdata)
  }
  

  