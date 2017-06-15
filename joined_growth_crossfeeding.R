#
# Script to determine pair-wise joint exchanges for metabolic models
# Returns a data.table with exchanged metabolites.

library(sybil)
library(foreach)
library(doParallel)
#library(parallel)
library(data.table)
library(lubridate)

source("./join_models.R")
source("./coupling.R")
source("./check_crossfeeding.R")
source("./correct_common_errors_agora2.R")

writeProcessToFile <- function(ts,i,n,file="process.txt") {
  tp <- as.numeric(Sys.time()) - ts
  tr <- tp / i * (n-i)
  tp <- seconds_to_period(tp)
  tr <- seconds_to_period(tr)
  line <- paste0(day(tp)*24+hour(tp),":",ifelse(minute(tp)>9,minute(tp),paste0(0,minute(tp))),":",
                 ifelse(round(second(tp))>9,round(second(tp)),paste0(0,round(second(tp)))),"   ",  
                 i," / ",n,"\t",
                 "(",round(i/n*100,digits=1),"%)\t",
                 day(tr)*24+hour(tr)," h ",minute(tr)," m remaining\n")
  cat(line, file = file, append = T)
}

file.create("process.txt")
file.create("paired_model.log")

#sybil::SYBIL_SETTINGS("SOLVER","glpkAPI"); ok <- 5
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1

# download of agora 1.0.1 constraint models from http://vmh.uni.lu
cat("Loading AGORA models...\n", file = "process.txt", append = T)
agora <- readRDS("./agora_western-paper.RDS")
cat("Correcting models...\n", file = "process.txt", append = T)
agora <- lapply(agora,correct_common_errors_agora)
names <- sapply(agora, function(x){x@mod_desc})
comb <- combn(names, 2)
#set.seed(1)
#ids <- sample(x = 1:ncol(comb),750)
#comb <- comb[,1:40000]

# determine single growth rate for each model
cat("Calculating single growth rates...\n", file = "process.txt", append = T)
single_growth <- sapply(agora, function(mod){
  sol <- optimizeProb(mod, retOptSol=F)
  ifelse(sol$stat==ok, sol$obj, 0)
})

no_cores <- 30 # choose core for parallel computing
registerDoParallel(no_cores)
DFinteract <- data.table()

ts <- as.numeric(Sys.time())
cat("Simulating co-growth...\n", file = "process.txt", append = T)
print(system.time(DFinteract <- foreach(i=1:ncol(comb), .combine=function(...) rbindlist(list(...)), .multicombine=T) %dopar% {
  if(i %% (no_cores*15-1) == 0) writeProcessToFile(ts,as.numeric(sub(" paired_model.log","",system("wc -l paired_model.log",intern=T),fixed=T)),ncol(comb))
  spec1 <- comb[1,i]; nr1 <- which(names==spec1)
  spec2 <- comb[2,i]; nr2 <- which(names==spec2)
  cat(paste0(Sys.time(),"\t",Sys.getpid(),"\t",i,"\t",spec1,"\t",spec2,"\n"), file = "paired_model.log", append = T)
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
  }else {data.table()}
}))

saveRDS(DFinteract, file = "DFinteract_paper_western_anaerob_SWcorrect_pt1.RDS")
