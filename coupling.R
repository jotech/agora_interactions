# assumes irreversible model
get_coupling_constraints <- function(modj, id1="X", id2="Y"){
  coupling <- list(react=list(), x=list(), lb=vector(), ub=vector(), rtype=vector())
  intern_rea <- grep("EX_", modj@react_id, invert = T)
  for(j in intern_rea){
    # get biomass for coupling
    if(length(grep(paste0("^",id1,"_"), modj@react_id[j])) > 0)
      obj <- grep(paste0("^", id1, "_biomass"), modj@react_id)
    if(length(grep(paste0("^",id2,"_"), modj@react_id[j])) > 0)
      obj <- grep(paste0("^", id2, "_biomass"), modj@react_id)
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
