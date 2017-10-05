# assumes irreversible model
get_coupling_constraints <- function(modj, id1="X", id2="Y"){
  coupling <- list(react=list(), x=list(), lb=vector(), ub=vector(), rtype=vector())
  intern_rea <- grep("EX_", modj@react_id, invert = T)
  allObj <- which(modj@obj_coef!=0) # get objectives
  obj1 <- allObj[grep(paste0("^",id1,"_"), react_id(modj)[allObj])]
  obj2 <- allObj[grep(paste0("^",id2,"_"), react_id(modj)[allObj])]
  for(j in intern_rea){
    # get biomass for coupling
    if(length(grep(paste0("^",id1,"_"), modj@react_id[j])) > 0)
      obj <- obj1
    if(length(grep(paste0("^",id2,"_"), modj@react_id[j])) > 0)
      obj <- obj2
    reversible <- modj@lowbnd[j] < 0
    # add coupling backwards direction
    if(reversible) coupling <- add_coupling(coupling, react=c(j, obj), x=c(1,400), lb=0.01, ub=NA, rtype="L")
    # add coupling forward direction
    coupling <- add_coupling(coupling, react=c(j, obj), x=c(1,-400), lb=NA, ub=0.01, rtype="U")
  }
  return(coupling)
}

add_coupling <- function(coupling, react, x, lb, ub, rtype){
  index <- length(coupling$react) + 1
  coupling$react[[index]] <- react
  coupling$x[[index]]     <- x
  coupling$lb[index]      <- lb
  coupling$ub[index]      <- ub
  coupling$rtype[index]   <- rtype
  return(coupling)
}

get_coupling_constraints_mult <- function(modj) {
  coupling <- list(react=list(), x=list(), lb=vector(), ub=vector(), rtype=vector())
  intern_rea  <- grep("EX_", modj@react_id, invert = T)
  allObj.ind  <- which(modj@obj_coef!=0)
  allObj.name <- react_id(modj)[allObj.ind]
  for(j in intern_rea){
    # get corresponding objective reaction for reaction j
    m.id <- gsub("_.*", "", react_id(modj)[j])
    obj.ind <- grep(paste0("^",m.id,"_"),allObj.name)
    obj <- allObj.ind[obj.ind]
    
    reversible <- modj@lowbnd[j] < 0
    # add coupling backwards direction
    if(reversible) coupling <- add_coupling(coupling, react=c(j, obj), x=c(1,400), lb=0.01, ub=NA, rtype="L")
    # add coupling forward direction
    coupling <- add_coupling(coupling, react=c(j, obj), x=c(1,-400), lb=NA, ub=0.01, rtype="U")
  }
  return(coupling)
}
