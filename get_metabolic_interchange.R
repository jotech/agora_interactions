#
# Get a quantitive measure for the metabolic interchange of individual metabolites
#
get_metabolic_interchange <- function(modj, solj) {
  int.ext.rxns <- grep("^M[0-9]+_EX_",modj@react_id)
  ext.ext.rxns <- grep("^EX_",modj@react_id)
  
  dt <- data.table(rxn = modj@react_id[int.ext.rxns],
                   flux = solj$fluxes[int.ext.rxns])
  dt[, rxn := gsub("M[0-9]+_","",rxn)]
  dt[, flux := abs(flux)]
  dt <- dt[,.(flux = sum(flux)),by="rxn"]
  
  # external exchange
  dt2 <- data.table(rxn = modj@react_id[ext.ext.rxns],
                    o.flux = solj$fluxes[ext.ext.rxns])
  
  dt <- merge(dt,dt2,by="rxn")
  dt[,flux := flux-abs(o.flux)]
  #dt <- dt[flux > 1e-6]
  
  dt[order(-flux)]
}


