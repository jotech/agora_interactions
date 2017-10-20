#source("./join_mult_models.R")
#source("./coupling.R")
#source("./get_metabolic_interchange.R")

#
# simulate the growth of an AGORA-Community with fixed biomass ratios between community members
# spec.ratio - A data.frame or data.table with two coloums: 1. spec (Agora model name), 2. ratio
#
simulate_agora_commmunity <- function(agora.list, spec.ratio, pFBAcoeff = 1e-6) {
  require(data.table)
  
  ind <- which(names(agora.list) %in% spec.ratio$spec)

  # construct joined model 
  ag.joined <- join_mult_models(agora[ind])
  
  # introduce new objective reaction with fixed biomass ratios
  cat("Introducing new objective reaction with fixed biomass ratios...\n")
  ag.joined$model.IDs <- merge(ag.joined$model.IDs,spec.ratio,by.x="model.name",by.y="spec")
  bm.mets <- paste0(ag.joined$model.IDs$model.id,"_biomass[c]")
  ag.joined$modj <- addReact(model = ag.joined$modj, id = "EX_BIOMASS", 
                             met = bm.mets, 
                             Scoef = -ag.joined$model.IDs$ratio, 
                             reversible = F, 
                             reactName = "joined Biomass with fixed ratios")
  ag.joined$modj@obj_coef <- rep(0, length(ag.joined$modj@react_id))
  ag.joined$modj@obj_coef[grep("EX_BIOMASS", ag.joined$modj@react_id)] <- 1
  # block individual Biomass outflow reactions
  ag.joined$modj@uppbnd[grep("M[0-9]+_EX_biomass",ag.joined$modj@react_id)] <- 0
  
  # Perform Optimization
  cat("Getting coupling constraints...\n")
  coupling <- get_coupling_constraints_mult(ag.joined$modj)
  
  cat("Applying coupling constraints...\n")
  #modj_warm <- sysBiolAlg(ag.joined$modj, algorithm = "mtfEasyConstraint", easyConstraint=coupling)
  modj_warm <- sysBiolAlg(ag.joined$modj, algorithm = "mtfEasyConstraint2", easyConstraint=coupling, pFBAcoeff = pFBAcoeff)
  
  cat("Performing pFBA...\n")
  solj1 <- optimizeProb(modj_warm)
  
  # Get community growth
  out.gr <- solj1$fluxes[grep("EX_BIOMASS", ag.joined$modj@react_id)]
  
  # Get metabolic interchange
  met.interchange <- get_metabolic_interchange(ag.joined$modj, solj1)
  
  return(list(modj = ag.joined$modj,
              solj = solj1,
              model.IDs = ag.joined$model.IDs,
              community.growth = out.gr,
              met.interchange = met.interchange))
  
}
