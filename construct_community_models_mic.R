library(foreach)
library(doParallel)



construct_commuity_models_mic <- function(mbo, agora, min.freq=0.001) {
  #start time
  strt<-Sys.time()
  
  ag.rel <- t(t(mbo@agora.table[-1,])/colSums(mbo@agora.table[-1,]))
  ag.rel <- data.table(as.table(ag.rel))
  colnames(ag.rel) <- c("spec","sample","ratio")
  ag.rel <- ag.rel[ratio >= min.freq]
  
  rat.tab <- list()
  for(i in 1:ncol(mbo@agora.table)) {
    rat.tab[[i]] <- ag.rel[sample==colnames(mbo@agora.table)[i],.(spec,ratio)]
  }
  
  agora <<- agora
  
  cl<-makeCluster(7)
  registerDoParallel(cl)
  
  #out <- foreach(i = 1:ncol(mbo@agora.table)) %dopar% {
  out <- foreach(i = 1:30) %dopar% {
    source("/home/swaschina/workspace-kiel/2017-agora_interactions_github/agora_interactions/addMultiReact.R")
    source("/home/swaschina/workspace-kiel/2017-agora_interactions_github/agora_interactions/coupling.R")
    source("/home/swaschina/workspace-kiel/2017-agora_interactions_github/agora_interactions/join_mult_models.R")
    source("/home/swaschina/workspace-kiel/2017-agora_interactions_github/agora_interactions/simulate_agora_commmunity.R")
    simulate_agora_commmunity(agora.list = agora, spec.ratio = rat.tab[[i]])
  }
  
  names(out) <- colnames(mbo@agora.table)[1:30]
  
  cat(paste0("Elapsed time: ",Sys.time()-strt,"\n"))
  stopCluster(cl)
  
  return(out)
}