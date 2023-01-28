
## @Ni, X., Geng, B., Zheng, H., Shi, J., Hu, G., & Gao, J. (2021). Accurate estimation of single-cell differentiation potency based on network topology and gene ontology information. IEEE/ACM transactions on computational biology and bioinformatics.

##############################################################################

##  description: This function computes the DPOR value for each cell with the scRNA-Seq expression matrix, curvature matrix and gene-to-gene functional similarity matrix. 

##  usage: CompDPOR_Mm(exp,curve,km)

##  arguments: 
##  exp: The output value of function Process_Mm.
##  curve: The curvature matrix transformed from the output value of CompORC.py.
##  km: The pre-compiled pairwise Kappa similarity matrix on Gene Ontology of human genes.

##  value:
##  DPOR_nor: The normalized single-cell potency measure computing by scRNA-seq, Ollivier Curvature and Gene Ontology similarity scores.

##############################################################################


CompDPOR_Mm<-	
  function(exp, curve, km){
    
    require("AnnotationDbi")
    require("org.Hs.eg.db")
    
    exp <- exp[!duplicated(exp[,1]),]
    rownames(exp) <- exp[,1]
    exp <- exp[,-1]
    
    curve <- curve[!duplicated(curve[,1]),]
    rownames(curve) <- curve[,1]
    curve <- curve[,-1]
    
    r = as.character(rownames(exp))
    geneIDselect1.0 <- select(org.Hs.eg.db, keys=r, columns="SYMBOL", keytype="ENTREZID")
    geneIDselect1.1 <- geneIDselect1.0[!duplicated(geneIDselect1.0[,2]),]
    geneIDselect1.2 <- geneIDselect1.1[!is.na(geneIDselect1.1[,2]),]
    
    exp = exp[!duplicated(geneIDselect1.0[,2]),]
    exp = exp[!is.na(geneIDselect1.1[,2]),]
    rownames(exp) = geneIDselect1.2[,2]
    
    curve = curve[!duplicated(geneIDselect1.0[,2]),!duplicated(geneIDselect1.0[,2])]
    curve = curve[!is.na(geneIDselect1.1[,2]),!is.na(geneIDselect1.1[,2])]
    rownames(curve) = geneIDselect1.2[,2]
    colnames(curve) = geneIDselect1.2[,2]
    
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

setwd("/Users/renyu/文稿/cellpotency-main")

load("lib/GO/hs_km.Rda")
exp <- read.csv("interface/exp/Briggs_exp.csv", header=TRUE)
curve <- read.csv("interface/curve/Briggs_curve.csv", header=TRUE)

score <- CompDPOR_Mm(exp, curve, km)


