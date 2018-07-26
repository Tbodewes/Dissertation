library(parallel)
library(data.table)

computeMSE <- function(varName, dat, dag.fit, train, impute = FALSE, particles = 100){
  colNum <- which(names(dat) == varName)
  
  
  if(train){
    fits <- dag.fit[[colNum]]$fitted.values
    actuals <- dat[[colNum]]
  } else if(impute){
    #Only impute when an actual is available (otherwise it is a waste of time)
    complete.actuals <- complete.cases(dat[, colNum])
    dat.toImpute <- dat[complete.actuals,]
    actuals <- dat.toImpute[[colNum]]
    
    #Simultaneously impute missing input and predict output
    dat.toImpute[, colNum] <- as.numeric(NA)
    imputed <- impute(dag.fit, dat.toImpute, method = "bayes-lw",
                      n = particles, debug = FALSE)
    fits <- imputed[[colNum]]
  } else {
    temp <- na.omit(dat)
    fits <- predict(object = dag.fit, node = varName, data = temp[, -colNum],
                    method = "bayes-lw", n = particles)
    actuals <- temp[[colNum]]
  }
  
  RMSE <- sqrt(mean( (fits - actuals)^2 , na.rm = TRUE))
  normRMSE <- RMSE/(max(actuals, na.rm = TRUE) - min(actuals, na.rm = TRUE))
  return(normRMSE)
}

parImpute <- function(dag.fit, dat, particles){
  
  no.cores <- detectCores() - 1
  cl <- makeCluster(no.cores)
  invisible(clusterEvalQ(cl, library(bnlearn)))
  no.obs <- ceil(nrow(dat)/no.cores)
  
  dat.list <- split(dat, (as.numeric(rownames(dat))-1)%/%no.obs)
  
  impute.list <- parLapply(cl, dat.list, impute, object = dag.fit, 
                           method = "bayes-lw", n = particles)
  
  imputed <- rbindlist(impute.list)
  stopCluster(cl)
  return(imputed)
}