setwd("/Users/renyu/Downloads/cellpotency-main/")

# Packages installing
source("code/environment.R")

# Data processing
source("code/M2H.R")
source("code/Process_Mm.R")
counts <- read.csv("examples/GSE97390_Briggs.csv", header=TRUE) # Dataset:Briggs
counts <- counts[-1,]
counts <- M2H(counts)
res0 <- Process_Mm(counts)

# Calculation of DPOR
source("code/DoIntegPPI.R")
source("code/E2S.R")
source("code/CompDPOR.R")
load("lib/PPI_Network/net17Jan16.rda")
load("lib/GO/hs_km.Rda")
res1 <- DoIntegPPI(res0, net17Jan16.m)
exp <- res1[["expMC"]]
adj <- res1[["adjMC"]]
## "Briggs_curve.csv" can be obrained by running CompORC.py
curve <- read.csv("examples/Briggs_curve.csv", header=TRUE)
rownames(curve) <- curve[,1]
curve <- curve[,-1]
res2 <- E2S(exp, curve)
score <- CompDPOR(res2[["exp"]], res2[["curve"]], km)
write.csv(score, "examples/Briggs_score.csv")

