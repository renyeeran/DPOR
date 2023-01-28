setwd("/Users/renyu/Downloads/cellpotency-main/")

# Packages installing
source("code/environment.R")

# Data processing
source("code/Process_Hs.R")
load("lib/PPI_Network/net17Jan16.rda")
counts <- read.csv("examples/GSE75748_Chu1.csv", header=TRUE) #Dataset:Chu1
res0 <- Process_Hs(counts, net17Jan16.m)

# Calculation of DPOR
source("code/DoIntegPPI.R")
source("code/CompDPOR.R")
load("lib/GO/hs_km.Rda")
res1 <- DoIntegPPI(res0[["exp"]], res0[["adj"]])
exp <- res1[["expMC"]]
adj <- res1[["adjMC"]]
## "Chu1_curve.csv" can be obtained by running CompORC.py
curve <- read.csv("examples/Chu1_curve.csv", header=TRUE)
rownames(curve) <- curve[,1]
curve <- curve[,-1]
score <- CompDPOR(exp, curve, km)
write.csv(score, "examples/Chu1_score.csv")
