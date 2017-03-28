#
# Script to determine pair-wise joint growth rates for metabolic models
# Returns a data.frame with single and joint growth rates.

library(sybil)
library(foreach)
library(doParallel)
library(parallel)

source("./join_models.R")
source("./coupling.R")
sybil::SYBIL_SETTINGS("SOLVER","glpkAPI") # if changed to another solver please also change sol$stat

# download of agora 1.0.1 constraint models from http://vmh.uni.lu
agora <- readRDS("./agora_western-paper.RDS")
names <- sapply(agora, function(x){x@mod_desc})
comb <- combn(rel_agora, 2)

# determine single growth rate for each model
single_growth <- sapply(agora, function(mod){
  sol <- optimizeProb(mod, retOptSol=F)
  ifelse(sol$stat==5, sol$obj, 0)
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
  solj <- optimizeProb(modj, algorithm=("fbaEasyConstraint"), easyConstraint=coupling, retOptSol=F)
  if(solj$stat==5){objj <- solj$fluxes[which(modj@obj_coef!=0)]; obj1j <- objj[1]; obj2j <- objj[2]
  }else {obj1j <- 0; obj2j <- 0}
  data.frame(spec1, spec2, obj1, obj2, obj1j, obj2j)
}))
stopCluster(cl)
