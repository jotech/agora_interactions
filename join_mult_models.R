#source("./addMultiReact.R")
#
# Joins multiple metabolic models (imported into R with sybilSBML)
#
join_mult_models <- function(model.list) {
  require(data.table)
  n <- length(model.list)
  
  model.IDs <- data.table(model.name = unlist(lapply(model.list,FUN = function(x) x@mod_desc)),
                          model.id = as.character(1:n))
  model.IDs$model.id <- paste0("M",model.IDs$model.id)
  
  # get a table of all exchange reactions and their lower bounds
  cat("Extracting all metabolite exchange (i.e. \"EX\") reactions...\n")
  ex.rxns <- data.table(model.name = character(0L), react_id = character(0L), lb = double(0L))
  for(i in 1:n) {
    ex.ind <- grep("^EX_",model.list[[i]]@react_id)
    tmp <- data.table(model.name = model.list[[i]]@mod_desc,
                      react_id = model.list[[i]]@react_id[ex.ind],
                      lb = model.list[[i]]@lowbnd[ex.ind])
    ex.rxns <- rbind(ex.rxns,tmp)
  }
  ex.rxns[,model.name := NULL]
  setkey(ex.rxns,NULL)
  ex.rxns <- unique(ex.rxns)
  if(any(duplicated(ex.rxns$react_id))) {
    stop(paste0("unequal lower bounds for at least one exchange reaction. "))
  } 
  ex.rxns[,met := gsub("^EX_","",react_id)]
  ex.rxns <- ex.rxns[met != "biomass(e)"]
  
  # rename reactions and metabolites
  cat("Renaming reaction-, metabolite, and compartment IDs...\n")
  for(i in 1:n){
    model.list[[i]]@react_id    <- paste0(model.IDs$model.id[i],"_", model.list[[i]]@react_id)
    model.list[[i]]@met_id      <- paste0(model.IDs$model.id[i],"_", model.list[[i]]@met_id)
    model.list[[i]]@mod_compart <- paste0(model.IDs$model.id[i],"_", model.list[[i]]@mod_compart)
  }
  
  
  # initiating joined/community model.
  compNames <- character(0)
  for(i in 1:n)
    compNames <- c(compNames,model.list[[i]]@mod_compart)
    
  modj <- modelorg(id = "joined.mod",
                   name = paste0("Number of species in community: ",n),
                   compartment = compNames)
  
  # add reactions from single models to joined/community model
  cat("Adding reactions from individual species to joined/community model...\n")
  for(i in 1:n) {
    cat(paste0("\r",i,"/",n))
    n.mets <- length(modj@met_comp)
    max.c <- ifelse(n.mets>0,max(modj@met_comp),0)
    model.list[[i]]@met_comp <- model.list[[i]]@met_comp + as.integer(max.c)
    modj <- addMultiReact(modj, ids=model.list[[i]]@react_id, src = model.list[[i]])
  }
  
  
  # add new external exchanges
  cat("\rAdding new external exchanges...\n")
  mod_compart(modj) <- c(mod_compart(modj), "e")
  Nex <- nrow(ex.rxns)
  modj <- addMultiReact(model=modj, ids=ex.rxns$react_id, mets=ex.rxns$met, Scoefs = rep(-1, Nex), 
                        reversible = TRUE, lb=ex.rxns$lb, metComp = rep("e", Nex))
  mod_name(modj) <- "A SW - joined model."

  # setting up original exchange interactions to interact with new common "e" compartment
  # e.g. "M1_ac(e) <->"  ==> "M1_ac(e) <-> ac(e)" + removing lower bnd (new: -1000)
  cat("Modifying original exchange reactions to interact with new common \"e\" compartment...\n")
  r.ind <- grep("^M[0-9]+_EX_",modj@react_id)
  tmp.ind <- grep("^M[0-9]+_EX_biomass",modj@react_id)
  r.ind <- r.ind[!(r.ind %in% tmp.ind)]
  tmp.mets <- gsub("^M[0-9]+_EX_","",modj@react_id[r.ind])
  m.ind.e  <- match(tmp.mets, modj@met_id)
  modj@lowbnd[r.ind] <- rep(-1000,length(r.ind))
  for(i in 1:length(r.ind))
    modj@S[m.ind.e[i],r.ind[i]] <- 1
  
  return(list(model.IDs=model.IDs,
              modj = modj,
              ex.rxns = ex.rxns))
}