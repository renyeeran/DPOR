# install relevant R packages:

if (!requireNamespace("BiocManager", quietly=TRUE)){
  install.packages("BiocManager")
}

package.list=c("AnnotationDbi","preprocessCore","homologene",
               "org.Hs.eg.db", "org.Mm.eg.db")
for (package in package.list){
  if (!requireNamespace(package, quietly=TRUE)){
    BiocManager::install(package)
  }
}

if (!requireNamespace("igraph", quietly=TRUE)){
  install.packages("igraph")
}

