# addMultiReact
#
# extension for sybil to allow to add several reaction from one model to another
# 
# TODO: add handling of gprAssoc and subSystem

addMultiReact <- function(model, 
                          ids,
                          src = NA, # if defined then all other parameters will be taken from this model
                          mets = NA,   # vector needed
                          Scoefs = NA, # vector needed
                          reversible = FALSE,
                          lb = 0,
                          ub = SYBIL_SETTINGS("MAXIMUM"),
                          obj = 0,
                          subSystem = NA,
                          gprAssoc = NA,
                          reactName = NA,
                          metName = NA,
                          metComp = NA) {

  
  # ------------------------------------------------------------------------ #
  # check arguments
  # ------------------------------------------------------------------------ #
  
  if (!is(model, "modelorg")) {
    stop("needs an object of class modelorg!")
  }
  
  stopifnot(checkVersion(model))
  
  if (!is(src, "modelorg")){
    Nids <- length(ids)
    if (length(mets) != Nids | length(Scoefs) != Nids){
      stop("all arguments have to be provided for each reaction")
    }
  }
  
  check_mets <- sapply(1:length(ids), function(i){length(mets[i]) == length(Scoefs[i])})
  if (any(is.na(check_mets) | !all(check_mets))) {
    stop("in each reaction arguments 'met' and 'Scoef' must have the same length")
  }

  
  # ------------------------------------------------------------------------ #
  # get main parameters from source model
  # ------------------------------------------------------------------------ #
  
  if(is(src, "modelorg")){
    reactInd  <- match(ids, react_id(src))
    reactName <- react_name(src)[reactInd]
    Crev      <- react_rev(src)[reactInd]
    lb        <- lowbnd(src)[reactInd]
    ub        <- uppbnd(src)[reactInd]
    obj       <- obj_coef(src)[reactInd]
    subSystem <- NA #subSys(src)[reactInd,] # subsystems not implemented yet
    gprAssoc  <- NA # gene association not implemented yet
    mets      <- sapply(reactInd, function(r){met_id(src)[which(S(src)[,r]!=0)]})
    Scoefs    <- sapply(reactInd, function(r){S(src)[,r][which(S(src)[,r]!=0)]})
    met       <- unique(unlist(mets))
    metInd    <- match(met, met_id(src))
    metName   <- met_name(src)[metInd]
    metComp   <- met_comp(src)[metInd]
  }else{
    met       <- unique(unlist(mets))
    nR        <- length(ids)
    if(length(obj)==1 & nR>1) obj <- rep(obj, nR)
    if(length(lb)==1 & nR>1)  lb  <- rep(lb,  nR)
    if(length(ub)==1 & nR>1)  ub  <- rep(ub,  nR)
    if(length(reversible)==1 & nR>1) reversible <- rep(reversible, nR)
    Crev      <- ifelse( ( (ub > 0) && (lb < 0) ) && (!isTRUE(reversible)), TRUE, reversible)
  }
  
  
  # ------------------------------------------------------------------------ #
  # check, if we need to add columns and/or rows
  # ------------------------------------------------------------------------ #
  
  # reactions
  colInd    <- match(ids, react_id(model))
  newR      <- which(is.na(colInd))
  nCols     <- react_num(model)
  nNewCols  <- length(newR)
  addCol    <- FALSE
  
  for(i in seq(along = newR)){
    addCol <- TRUE
    nCols <- nCols + 1
    colInd[newR[i]] <- nCols
  }
  
  
  # metabolites
  
  rowInd  <- match(met, met_id(model))
  
  newM    <- which(is.na(rowInd))
  nRows   <- met_num(model)          # number of rows in the model
  nNewRows<- length(newM)            # number of new rows
  addRow  <- FALSE
  
  for (i in seq(along = newM)) {
    addRow   <- TRUE
    nRows    <- nRows + 1
    rowInd[newM[i]] <- nRows
  }
  
  
  if ( (isTRUE(addCol)) || (isTRUE(addRow)) ) {    
    
    # -------------------------------------------------------------------- #
    # make a new model
    # -------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------- #
    # data structures
    
    newmet_num      <- met_num(model)
    newmet_id       <- met_id(model)
    newmet_name     <- met_name(model)
    newmet_comp     <- met_comp(model)
    newmet_single   <- met_single(model)
    newmet_de       <- met_de(model)
    
    newreact_num    <- react_num(model)
    newreact_rev    <- react_rev(model)
    newreact_id     <- react_id(model)
    newreact_name   <- react_name(model)
    newreact_single <- react_single(model)
    newreact_de     <- react_de(model)
    newlowbnd       <- lowbnd(model)
    newuppbnd       <- uppbnd(model)
    newobj_coef     <- obj_coef(model)
    
    newgprRules     <- gprRules(model)
    newgenes        <- genes(model)
    newgpr          <- gpr(model)
    newallGenes     <- allGenes(model)
    newrxnGeneMat   <- rxnGeneMat(model)
    newsubSys       <- subSys(model)
    
    newS            <- S(model)
    
    newMetAttr <- met_attr(model)
    newReactAttr <- react_attr(model)
    newCompAttr <- comp_attr(model)
    newModAttr <- mod_attr(model)
    
    
    if (isTRUE(addRow)) {
      
      # new number of metabolites
      newmet_num  <- nRows
      
      # new metabolite id's
      newmet_id   <- append(met_id(model), met[newM])
      
      # new metabolite names
      if (any(is.na(metName))) {
        newmet_name <- append(met_name(model), met[newM])
      }
      else {
        newmet_name <- append(met_name(model), metName[newM])
      }
      
      # new metabolite compartments
      if (any(is.na(metComp))) {
        newmet_comp <- append(met_comp(model), rep(NA, nNewRows))
      }
      else {
        if (is(metComp, "numeric")) {
          newmet_comp <- append(met_comp(model), metComp[newM])
        }
        else {
          newmet_comp <- append(met_comp(model),
                                match(metComp[newM],
                                      mod_compart(model)))
        }
      }
      
      # singleton and dead end metabolites (not checked!)
      newmet_single <- append(met_single(model), rep(NA, nNewRows))
      newmet_de     <- append(met_de(model),     rep(NA, nNewRows))
      
      # new rows in stoichiometric matrix
      newRows <- Matrix::Matrix(0,
                                nrow = nNewRows,
                                ncol = react_num(model))
      newS <- Matrix::rBind(newS, newRows)
      
      # new met attrs
      if(ncol(newMetAttr) > 0){
        newMetAttr[nrow(newMetAttr)+1:nNewRows, ] <- NA
      }
    }
    
    if (isTRUE(addCol)) {                        # we add at least one column
      # new number of reactions
      newreact_num  <- nCols
      
      # new reaction ids
      newreact_id   <- append(react_id(model), ids[newR])
      
      # new reaction names
      if (any(is.na(reactName))) {
        newreact_name <- append(react_name(model), ids[newR])
      }
      else {
        newreact_name <- append(react_name(model), reactName[newR])
      }
      
      # reaction contains singleton or dead end metabolites (not checked!)
      newreact_single <- append(react_single(model), rep(NA, nNewCols))
      newreact_de     <- append(react_de(model),     rep(NA, nNewCols))
      
      # reversibility, lower and upper bounds, objective coefficient
      newreact_rev <- append(react_rev(model), Crev[newR])
      newlowbnd    <- append(lowbnd(model),    lb[newR])
      newuppbnd    <- append(uppbnd(model),    ub[newR])
      newobj_coef  <- append(obj_coef(model),  obj[newR])
      
      # new columns in stoichiometric matrix
      newColumns <- Matrix::Matrix(0,
                                ncol = nNewCols,
                                nrow = nrow(newS))
      
      newS <- Matrix::cBind(newS, newColumns)
      
      # new react Attr
      if(ncol(newReactAttr) > 0){
        newReactAttr[nrow(newReactAttr)+1:nNewCols, ] <- NA
      }
      
      # subsystems
      if (any(is.na(subSystem))) {
        ss <- subSys(model)
        if(ncol(ss)==0){ # if no subSys defined, rbind (see else) failed
          dim(ss) <- c(nrow(ss)+nNewCols, ncol(ss)) # new reactions are rows in this matrix
          newsubSys <- ss
        }
        else {
          newsubSys <- rbind(ss, rep(rep(FALSE, ncol(subSys(model))), nNewCols))
        }
      }
      # else {
      #   if (is(subSystem, "logical")) {
      #     newsubSys <- rBind(subSys(model), subSystem)
      #   }
      #   else {
      #     nSubsRow  <- colnames(subSys(model)) %in% subSystem
      #     newsubSys <- rBind(subSys(model), nSubsRow)
      #   }
      # }
      
      
      # gpr association
      if (ncol(rxnGeneMat(model)) > 0) {
        newrxnGeneMat   <- rbind(rxnGeneMat(model),
                                 rep(rep(FALSE, ncol(rxnGeneMat(model))), nNewCols))
      }
      else { #if (nrow(rxnGeneMat(model)) > 0) {
        newrxnGeneMat <- rxnGeneMat(model)
        dim(newrxnGeneMat) <- c(nrow(newrxnGeneMat)+nNewCols,
                                ncol(newrxnGeneMat))
      }
      # do above else always.
      
      if ( (is.na(gprAssoc)) || (gprAssoc == "") ) {
        if ((length(gprRules(model)) > 0)) {
          newgprRules     <- append(gprRules(model), rep("", nNewCols))
          newgenes        <- append(genes(model), rep("", nNewCols))
          newgpr          <- append(gpr(model), rep("", nNewCols))
        }
      }
      # else {
      #   gene_rule <- .parseBoolean(gprAssoc)
      #   
      #   geneInd <- match(gene_rule$gene, allGenes(model))
      #   
      #   # indices of new genes
      #   new_gene <- which(is.na(geneInd))
      #   
      #   # if we have new gene(s), add a column in rxnGeneMat and
      #   # gene name(s) to allGenes
      #   if (length(new_gene) > 0) {
      #     newallGenes <- append(allGenes(model),
      #                           gene_rule[["gene"]][new_gene])
      #     
      #     # update geneInd
      #     geneInd <- match(gene_rule[["gene"]], newallGenes)
      #     
      #     # if we have an empty modelorg object, we need to
      #     # initialize rxnGeneMat
      #     if (ncol(newrxnGeneMat) == 0) {
      #       newrxnGeneMat <- Matrix::Matrix(FALSE,
      #                                       nCols, max(geneInd))
      #     }
      #     else {
      #       for (i in seq(along = gene_rule[["gene"]][new_gene])) {
      #         newrxnGeneMat <- cBind(newrxnGeneMat,
      #                                rep(FALSE, nrow(newrxnGeneMat)))
      #       }
      #     }
      #   }
      #   
      #   # rxnGeneMat
      #   newrxnGeneMat[nCols, geneInd] <- TRUE
      #   
      #   # new rule
      #   newgpr <- append(gpr(model), gprAssoc)
      #   
      #   # genes per reaction
      #   newgenes <- append(genes(model), list(gene_rule$gene))
      #   newrule  <- gene_rule$rule
      #   
      #   # not needed for modelorg version 2.0
      #   #                for (j in 1 : length(geneInd)) {
      #   #                    pat  <- paste("x(", j, ")", sep = "")
      #   #                    repl <- paste("x[", geneInd[j], "]", sep = "")
      #   #    
      #   #                    newrule <- gsub(pat, repl, newrule, fixed = TRUE)
      #   #                }
      #   
      #   newgprRules <- append(gprRules(model), newrule)
      # }
    }
    
    # values for stoichiometric matrix
    newS[ , colInd] <- sapply(1:length(colInd), function(i){
      newCol <- rep(0, length = nrow(newS))
      curRInd <- colInd[i]
      #curMInd <- metInd[match(mets[[i]], met)]
      curMInd <- rowInd[match(mets[[i]], met)]
      newCol[curMInd] <- Scoefs[[i]]
      newCol
    })
    

    # -------------------------------------------------------------------- #
    # new model
    # -------------------------------------------------------------------- #
    
    if (is(model, "modelorg_irrev")) {
      mod_out <- modelorg_irrev(mod_id(model), mod_name(model))
      irrev(mod_out)     <- TRUE
      matchrev(mod_out)  <- append(matchrev(model), 0L)
      
      revReactId <- as.integer(max(irrev2rev(model))+1)
      irrev2rev(mod_out) <- append(irrev2rev(model), revReactId)
      rev2irrev(mod_out) <- rbind(rev2irrev(model), c(nCols, nCols))
    } else {
      mod_out <- modelorg(mod_id(model), mod_name(model))
    }
    
    mod_desc(mod_out)    <- mod_desc(model)
    mod_compart(mod_out) <- mod_compart(model)
    
    
    met_num(mod_out)      <- as.integer(newmet_num)
    met_id(mod_out)       <- newmet_id
    met_name(mod_out)     <- newmet_name
    met_comp(mod_out)     <- as.integer(newmet_comp)
    met_single(mod_out)   <- newmet_single
    met_de(mod_out)       <- newmet_de
    
    react_num(mod_out)    <- as.integer(newreact_num)
    react_rev(mod_out)    <- newreact_rev
    react_id(mod_out)     <- newreact_id
    react_name(mod_out)   <- newreact_name
    react_single(mod_out) <- newreact_single
    react_de(mod_out)     <- newreact_de
    lowbnd(mod_out)       <- newlowbnd
    uppbnd(mod_out)       <- newuppbnd
    obj_coef(mod_out)     <- newobj_coef
    
    gprRules(mod_out)     <- newgprRules
    genes(mod_out)        <- newgenes
    gpr(mod_out)          <- newgpr
    allGenes(mod_out)     <- newallGenes
    rxnGeneMat(mod_out)   <- newrxnGeneMat
    subSys(mod_out)       <- newsubSys
    
    S(mod_out)            <- newS
    
    react_attr(mod_out) <- newReactAttr
    met_attr(mod_out) <- newMetAttr
    comp_attr(mod_out) <- newCompAttr
    mod_attr(mod_out) <- newModAttr
    
    
  } else{
    stop("Nothing to change")
  }
  
  check <- validObject(mod_out, test = TRUE)
  
  if (check != TRUE) {
    msg <- paste("Validity check failed:", check, sep = "\n    ")
    warning(msg)
  }
  
  return(mod_out)
  
}
