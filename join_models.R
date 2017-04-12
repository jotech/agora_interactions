#
# Joins two metabolic models (imported into R with sybilSBML)
#


library(sybil)
if (! exists("addMultiReact") ){
    source("./addMultiReact.R") # uses an extension of sybil to add several reactions at once
}

join_models <- function(mod1, mod2, id1="X", id2="Y"){
  allex1 <- findExchReact(mod1); allex1 <- allex1[grep("EX_", allex1@react_id),]
  allex2 <- findExchReact(mod2); allex2 <- allex2[grep("EX_", allex2@react_id),]
  ex1 <- allex1@react_pos
  ex2 <- allex2@react_pos
  ex  <- union(mod1@react_id[ex1], mod2@react_id[ex2])
  ex_met <- gsub("^EX_","",ex) # metabolite names in new pseudo compartment
  Ncomp1 <- length(mod_compart(mod1))
  Ncomp2 <- length(mod_compart(mod2))

  # check for common ground
  ComEx <- intersect(allex1@react_id, allex2@react_id)
  if(length(ComEx) ==0) stop("No common exchange reactions")

  # get exchange metabolites and lower bounds
  avail1  <- match(ex, allex1@react_id)
  avail2  <- match(ex, allex2@react_id)
  ex_met1 <- setNames(sapply(avail1, function(x){allex1@met_pos[x]}), ex)
  ex_met2 <- setNames(sapply(avail2, function(x){allex2@met_pos[x]}),ex)
  ex_lb1  <- sapply(avail1, function(x){allex1@lowbnd[x]})
  ex_lb2  <- sapply(avail2, function(x){allex2@lowbnd[x]})
  ex_lb   <- pmin(ex_lb1, ex_lb2, na.rm = T)

  # rename reactions and metabolites
  mod1@react_id <- paste0(id1,"_", mod1@react_id)
  mod1@met_id   <- paste0(id1,"_", mod1@met_id)
  mod2@react_id <- paste0(id2,"_", mod2@react_id)
  mod2@met_id   <- paste0(id2,"_", mod2@met_id)
  
  # set exchanges to work in inner (pseudo) compartment
  CompPs <- max(Ncomp1, Ncomp2)+1L # new unique compartment number
  mod1@mod_compart <- paste0(id1, "_",mod1@mod_compart)
  mod1@mod_compart <- c(mod1@mod_compart, "pseudo")
  new_mets <- list(); new_Scoefs <- list()
  for(i in seq_along(ex1)){
    j <- ex1[i]
    id <- gsub(paste0("^",id1,"_"), "", mod1@react_id[j])
    inMet  <- met_id(mod1)[unname(ex_met1[id])]
    exMet  <- ex_met[match(id, ex)] # new metabolite in common compartment
    stoich <- which(mod1@S[,j]!=0)
    Scoef <- abs(mod1@S[stoich,j])
    new_mets[[i]]   <- c(inMet, exMet)
    new_Scoefs[[i]] <- c(-Scoef, Scoef)
  }
  mod1 <- addMultiReact(mod1, ids=mod1@react_id[ex1], mets=new_mets, Scoefs=new_Scoefs)
  met_comp(mod1)[which(is.na(met_comp(mod1)))] <- CompPs
  mod1@lowbnd[ex1] <- -1000
  
  mod2@mod_compart <- paste0(id2, "_",mod2@mod_compart)
  mod2@mod_compart <- c(mod2@mod_compart, "pseudo")
  new_mets <- list(); new_Scoefs <- list()
  for(i in seq_along(ex2)){
    j <- ex2[i]
    id <- gsub(paste0("^",id2,"_"), "", mod2@react_id[j])
    inMet  <- met_id(mod2)[unname(ex_met2[id])]
    exMet  <- ex_met[match(id, ex)] # new metabolite in common compartment
    stoich <- which(mod2@S[,j]!=0)
    Scoef <- abs(mod2@S[stoich,j])
    new_mets[[i]]   <- c(inMet, exMet)
    new_Scoefs[[i]] <- c(-Scoef, Scoef)
  }
  mod2 <- addMultiReact(mod2, ids=mod2@react_id[ex2], mets=new_mets, Scoefs=new_Scoefs)
  met_comp(mod2)[which(is.na(met_comp(mod2)))] <- CompPs
  mod2@lowbnd[ex2] <- -1000
  
  # add reactions to one model
  if(Ncomp1+1 == CompPs){
    mod2@met_comp <- match(mod2@met_comp, c(rep(0, Ncomp1), CompPs, 1:Ncomp2))
    compall <- c(mod_compart(mod1), setdiff(mod_compart(mod2), "pseudo"))
  }else{
    mod1@met_comp <- match(mod1@met_comp, c(rep(0, Ncomp2), CompPs, 1:Ncomp1))  
    compall <- c(mod_compart(mod2), setdiff(mod_compart(mod1), "pseudo"))
  }
  modj <- addMultiReact(mod1, ids=mod2@react_id, src = mod2)
  mod_compart(modj) <- compall
  
  # add new external exchanges
  mod_compart(modj) <- c(mod_compart(modj), "e")
  Nex <- length(ex)
  modj <- addMultiReact(model=modj, ids=ex, mets=unname(ex_met), Scoefs = rep(-1, Nex), 
                        reversible = TRUE, lb=ex_lb, metComp = length(modj@mod_compart))
  mod_name(modj) <- paste("joint model:",mod_name(mod1), "+",mod_name(mod2))
  return(modj)  
}