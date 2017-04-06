#
# Joins two metabolic models (imported into R with sybilSBML)
#


library(sybil)
if (! exists("addMultiReact") ){
    source("./addMultiReact.R") # uses an extension of sybil to add several reactions at once
}

join_models <- function(mod1, mod2, id1="X", id2="Y"){
  ex1 <- grep("EX_", mod1@react_id)
  ex2 <- grep("EX_", mod2@react_id)
  ex  <- union(mod1@react_id[ex1], mod2@react_id[ex2])
  # get exchange metabolites and lower bounds
  ex_met <- vector(); ex_lb <- vector()
  for(x in ex){
    pos <- which(mod1@react_id %in% x)
    pos2<- which(mod2@react_id %in% x)
    # add lower bounds only one time
    if(length(pos)>0 & length(pos2)>0){ # exchange is in mod1 and mod2
      ex_met <- c(ex_met, mod1@met_id[which(mod1@S[,pos]!=0)])
      lb_min <- min(mod1@lowbnd[pos], mod1@lowbnd[pos])
      ex_lb  <- c(ex_lb, lb_min)
    } else if(length(pos)>0) { # exchange is only in mod1
      ex_met <- c(ex_met, mod1@met_id[which(mod1@S[,pos]!=0)])
      ex_lb  <- c(ex_lb, mod1@lowbnd[pos]) 
    } else{ # exchange is only in mod2
      pos <- which(mod2@react_id %in% x)
      ex_met <- c(ex_met, mod2@met_id[which(mod2@S[,pos]!=0)])
      ex_lb  <- c(ex_lb, mod2@lowbnd[pos])
    } 
  }
  
  # rename reactions and metabolites
  mod1@react_id <- paste0(id1,"_", mod1@react_id)
  mod1@met_id   <- paste0(id1,"_", mod1@met_id)
  mod2@react_id <- paste0(id2,"_", mod2@react_id)
  mod2@met_id   <- paste0(id2,"_", mod2@met_id)
  
  # set exchanges to work in inner (pseudo) compartment
  new_mets <- list(); new_Scoefs <- list()
  for(i in seq_along(ex1)){
    j <- ex1[i]
    id <- mod1@react_id[j]
    stoich <- which(mod1@S[,j]!=0)
    met<- mod1@met_id[stoich]
    Scoef <- mod1@S[stoich,j]
    new_mets[[i]]   <- c(met, gsub(paste0("^",id1,"_"), "", met))
    new_Scoefs[[i]] <- c(Scoef, -Scoef)
  }
  mod1 <- addMultiReact(mod1, ids=mod1@react_id[ex1], mets=new_mets, Scoefs=new_Scoefs)
  mod1@lowbnd[ex1] <- -1000
  
  new_mets <- list(); new_Scoefs <- list()
  for(i in seq_along(ex2)){
    j <- ex2[i]
    id <- mod2@react_id[j]
    stoich <- which(mod2@S[,j]!=0)
    met<- mod2@met_id[stoich]
    Scoef <- mod2@S[stoich,j]
    new_mets[[i]] <- c(met, gsub(paste0("^",id2,"_"), "", met))
    new_Scoefs[[i]] <- c(Scoef, -Scoef)
  }
  mod2 <- addMultiReact(mod2, ids=mod2@react_id[ex2], mets=new_mets, Scoefs=new_Scoefs)
  mod2@lowbnd[ex2] <- -1000
  
  # add reactions to one model
  modj <- addMultiReact(mod1, ids=mod2@react_id, src = mod2)
  
  # add new external exchanges
  modj <- addMultiReact(model=modj, ids=ex, mets=ex_met, Scoefs = rep(-1, length(ex)), 
                        reversible = TRUE, lb=ex_lb)
  return(modj)  
}
