## @Ni, X., Geng, B., Zheng, H., Shi, J., Hu, G., & Gao, J. (2021). Accurate estimation of single-cell differentiation potency based on network topology and gene ontology information. IEEE/ACM transactions on computational biology and bioinformatics.

##############################################################################

##  description: This function computes the DPOR value for each cell with the scRNA-Seq expression matrix, curvature matrix and gene-to-gene functional similarity matrix. 

##  usage: CompDPOR(exp,curve,km)

##  arguments: 
##  exp: The output value of function Process_Hs or E2S.
##  curve: The curvature matrix transformed from the output value of CompORC.py or of function E2S.
##  km: The pre-compiled pairwise Kappa similarity matrix on Gene Ontology of human genes.

##  value:
##  DPOR_nor: The normalized single-cell potency measure computing by scRNA-seq, Ollivier Curvature and Gene Ontology similarity scores.

##############################################################################


CompDPOR<-	
  function(exp, curve, km){
    
    commonID1.v <- intersect(rownames(km), rownames(curve))
    commonID.v <- intersect(commonID1.v, rownames(exp))
    match(commonID.v, rownames(curve)) -> map1.idx
    cint <- curve[map1.idx, map1.idx]
    match(commonID.v, rownames(km)) -> map2.idx
    kmint <- km[map2.idx, map2.idx]
    match(commonID.v, rownames(exp)) -> map3.idx
    expint <- exp[map3.idx, ]
    
    cr <- as.matrix(kmint*cint)
    m <- matrix(1, nrow=length(commonID.v), ncol=1)
    ORC <- cr%*%m  #summation
    
    DPOR <- t(expint)%*%ORC
    DPOR_nor <- DPOR/max(DPOR)  #normalization
    
    return(DPOR_nor)
  }

