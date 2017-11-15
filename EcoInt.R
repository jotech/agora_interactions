args <- commandArgs(trailingOnly = TRUE)
# 1: Path to AGORA Model List
# 2: SWcorr / nSWcorr
# 3: minGrowth / nminGrowth
# 4: Number of cores to use

args <- c()
args[1] <- "agora_HumanMilkDietMM"
args[2] <- "SWcorr"
args[3] <- "nminGrowth"
args[4] <- "7"

pFBAcoeff <- 1e-5

library(sybil)
library(foreach)
library(doParallel)
library(data.table)
library(lubridate)
source("/home/swaschina/workspace-kiel/Resources/AGORA/AGORA/correct_common_errors_agora2.R")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPX_PARAM_THREADS"=1L))

# Scripts for community modelling
source("join_mult_models.R")
source("addMultiReact.R")
source("coupling.R")
source("get_metabolic_interchange.R")


nr.cores <- as.numeric(args[4])

process.file <- paste0("results/PairInt_",args[1],"_",args[2],"_",args[3],"_process.txt")
log.file <- paste0("results/PairInt_",args[1],"_",args[2],"_",args[3],"_log.txt")

file.create(process.file)
file.create(log.file)

## Function for logging
writeProcessToFile <- function(ts,i,n,file = process.file) {
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

cat("Loading AGORA models...\n", file = process.file, append = T)
agora <- readRDS("/home/swaschina/workspace-kiel/Resources/AGORA/AGORA/agora_HumanMilkDietMM.RDS")
if(args[2] == "SWcorr") {
  cat("Correcting models...\n", file = process.file , append = T)
  agora <- lapply(agora,correct_common_errors_agora)
}

names <- sapply(agora, function(x){x@mod_desc})
comb <- combn(names, 2)
ids <- sample(x = 1:ncol(comb),100)
#comb <- comb[,1:298378] # end: 298378
comb <- comb[,ids] # end: 298378


# determine single growth rate for each model
cat("Calculating single growth rates...\n", file = process.file, append = T)
single_growth <- sapply(agora, function(mod){
  sol <- optimizeProb(mod, retOptSol=F)
  ifelse(sol$stat==ok, sol$obj, 0)
})

registerDoParallel(nr.cores)
DFinteract <- list()

ts <- as.numeric(Sys.time())
cat("Simulating co-growth...\n", file = process.file, append = T)
# 2->15 (two lines below)
cat(system.time(DFinteract <- foreach(i=1:ncol(comb)) %dopar% {
  if(i %% (nr.cores*2-1) == 0) writeProcessToFile(ts,as.numeric(sub(paste0(" ",log.file),"",
                                                                     system(paste0("wc -l ",log.file),intern=T),fixed=T))
                                                   ,ncol(comb))
  spec1 <- comb[1,i]; nr1 <- which(names==spec1)
  spec2 <- comb[2,i]; nr2 <- which(names==spec2)
  cat(paste0(Sys.time(),"\t",Sys.getpid(),"\t",i,"\t",spec1,"\t",spec2,"\n"), file = log.file, append = T)
  mod1 <- agora[[nr1]]; mod2 <- agora[[nr2]]
  
  if(args[3]=="minGrowth") {
    mod1@lowbnd[which(mod1@obj_coef==1)] <- single_growth[spec1]; mod2@lowbnd[which(mod2@obj_coef==1)] <- single_growth[spec2] 
  }
  
  obj1 <- single_growth[nr1]
  obj2 <- single_growth[nr2]
  ag.joined <- join_mult_models(list(mod1,mod2))
  ag.joined$coupling <- get_coupling_constraints_mult(ag.joined$modj)
  
  # construct solver object
  modj_warm <- sysBiolAlg(ag.joined$modj,
                          algorithm = "mtfEasyConstraint2",
                          easyConstraint=ag.joined$coupling,
                          pFBAcoeff = pFBAcoeff)
  # Performing pFBA
  solj <- optimizeProb(modj_warm)
  
  if(solj$stat==ok){objj <- solj$fluxes[which(ag.joined$modj@obj_coef!=0)]; obj1j <- objj[1]; obj2j <- objj[2]
  }else {obj1j <- 0; obj2j <- 0}
  
  exchange <- get_metabolic_interchange(ag.joined$modj, solj)
  exchange$spec1 <- spec1
  exchange$spec2 <- spec2
  exchange <- exchange[abs(flux) > 0 | abs(o.flux) > 0]
  
  out <- list(growth = data.table(spec1, spec2, obj1, obj2, obj1j, obj2j),
              exchange = exchange)
  
  out
}), file = process.file, append = T)

DFgrowth <- rbindlist(lapply(DFinteract, FUN = function(x) x$growth))
DFmets <- rbindlist(lapply(DFinteract, FUN = function(x) x$exchange))

