check_crossfeeding <- function(modj, solj, spec1, spec2, modj_warm, id1="X", id2="Y"){
  idx1    <- grep(paste0(id1,"_EX_"), react_id(modj))
  idx2    <- grep(paste0(id2,"_EX_"), react_id(modj))
  ex1     <- gsub(paste0("^",id1,"_"), "", react_id(modj)[idx1])
  ex2     <- gsub(paste0("^",id2,"_"), "", react_id(modj)[idx2])
  flux1   <- solj$fluxes[idx1]
  flux2   <- solj$fluxes[idx2]
  df1  <- data.table(ex1, idx1, flux1)
  df2  <- data.table(ex2, idx2, flux2)
  dfj  <- merge(df1, df2, by.x = "ex1", by.y = "ex2")
  if(nrow(dfj)>0){  
    cross<- dfj[which( dfj$flux1*dfj$flux2 != abs(dfj$flux1*dfj$flux2)), ]
    len  <- nrow(cross)
    if(len>0){  
      cross$spec1<- spec1
      cross$spec2<- spec2
      return(subset(cross, select = c(-idx1, -idx2)))
    }else{return(data.table())}
  }
}