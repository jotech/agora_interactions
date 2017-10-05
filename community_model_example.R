#
# An AGORA-Community model example
#
library(sybil)
library(data.table)

source("./join_mult_models.R")
source("./coupling.R")
source("/home/swaschina/workspace-kiel/Resources/AGORA/AGORA/correct_common_errors_agora2.R")

sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
agora <- readRDS("/home/swaschina/workspace-kiel/Resources/AGORA/AGORA/agora_western-paper.RDS")
agora <- lapply(agora,correct_common_errors_agora)

# Join 10 models
system.time({
  t.new <- join_mult_models(agora[sample(1:773,10)])
  #t.new <- join_mult_models(agora[c(158,449)])
  cat("Getting coupling constraints...\n")
  coupling <- get_coupling_constraints_mult(t.new$modj)
  cat("Applying coupling constraints...\n")
  modj_warm <- sysBiolAlg(t.new$modj, algorithm = "mtfEasyConstraint", easyConstraint=coupling)
  cat("Performing pFBA...\n")
  solj1 <- optimizeProb(modj_warm)})


# What is exchanged between community members?
get_metabolic_interchange(t.new$modj, solj1)