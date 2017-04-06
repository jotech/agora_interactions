#
# Script to determine pair-wise joint exchanges for metabolic models
# Returns a data.frame with exchanged metabolites.

library(sybil)
library(foreach)
library(doParallel)
library(parallel)

source("./join_models.R")
source("./coupling.R")
source("./check_crossfeeding.R")

#sybil::SYBIL_SETTINGS("SOLVER","glpkAPI"); ok <- 5
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1

# download of agora 1.0.1 constraint models from http://vmh.uni.lu
agora <- readRDS("./agora_western-paper.RDS")
names <- sapply(agora, function(x){x@mod_desc})
comb <- combn(agora, 2)

# determine single growth rate for each model
single_growth <- sapply(agora, function(mod){
  sol <- optimizeProb(mod, retOptSol=F)
  ifelse(sol$stat==ok, sol$obj, 0)
})

no_cores <- 4 # choose core for parallel computing
cl<-makeCluster(no_cores)
registerDoParallel(no_cores)
DFinteract <- data.frame()

print(system.time(DFinteract <- foreach(i=1:ncol(comb), .combine=rbind) %dopar% {
  spec1 <- comb[1,i]; nr1 <- which(names==spec1)
  spec2 <- comb[2,i]; nr2 <- which(names==spec2)
  mod1 <- agora[[nr1]]; mod2 <- agora[[nr2]]
  obj1 <- single_growth[nr1]
  obj2 <- single_growth[nr2]
  modj <- join_models(mod1, mod2)
  coupling <- get_coupling_constraints(modj)
  modj_warm <- sysBiolAlg(modj, algorithm = "mtfEasyConstraint", easyConstraint=coupling)
  solj <- optimizeProb(modj_warm)
  if(solj$stat==ok){objj <- solj$fluxes[which(modj@obj_coef!=0)]; obj1j <- objj[1]; obj2j <- objj[2]
  }else {obj1j <- 0; obj2j <- 0}
  if(obj1j>0 & obj2j>0){ check_crossfeeding(modj, solj, spec1, spec2, modj_warm)
  }else {data.frame()}
}))
stopCluster(cl)
