

#' TODO: Write tests for each layer!
#'
#' Layer 1: for a single dataset and penalty, fit DAG using NAL optimization and
#' optionally using structural EM. Extract number of queries and SHD for both
#'
#' Layer 2: for a single dataset, run layer 1 for all penalties. Collect results
#' in long dataframe with a column marking which penalty was used for each row
#'
#' Layer 3: for a given DAG, k and beta, generate 'replications' datasets and
#' run layer 2 for each. rbind results into single dataframe
#'
#' Layer 4: run layer 3 for all combinations of 1) sample
#' sizes and 2) betas. Store results in long dataframe with extra columns
#' marking each entry for 1, 2
#' 
#' Layer 5: run layer 4 for a vector of DAGs

#Collect results for single dataset and single penalty
experiment.layer1 <- function(dat, dag.true, penalty, str.em, nodes.ordered){
  if(is.null(nodes.ordered)){
    fit <- fit.dag(dat, penalty, parallel = FALSE, tabuSteps = 10)
    dag.fitted <- fit$dag
    no.queries <- fit$queries
  }
  else{
    dag.fitted <- fit.dag.ordered(nodes.ordered, max.parents = 3, dat, penalty, 
                                  parallel = FALSE)
    no.queries <- NA
    if(str.em){stop("Structural EM and ordered fitting should not be run in 
                     the same experiment")}
  }
  
  distance <- shd(dag.fitted, dag.true)
  
  if(str.em){
    dag.em <- structural.em(dat, maximize = "tabu", 
                            fit.args = list(replace.unidentifiable = TRUE))
    no.queries.em <- dag.em$learning$ntests
    distance.em <- shd(dag.em, dag.true)
  }
  else{distance.em <- NA; no.queries.em <- NA}
  
  return(list(distance, no.queries, distance.em, no.queries.em))
}

experiment.layer2 <- function(dat, dag.true, penalties, str.em, nodes.ordered){
  result.list <- lapply(penalties, experiment.layer1, dat = dat, 
                        dag.true = dag.true, str.em = str.em,
                        nodes.ordered = nodes.ordered)
  result.df <- data.frame(penalties, matrix(unlist(result.list), byrow = TRUE,
                                            nrow = length(penalties)))
  
  names(result.df) <- c("Penalty", "SHD.NAL", "Q.NAL", "SHD.EM", "Q.EM")
  return(result.df)
}

#TODO implement parallelization
experiment.layer3 <- function(dag.true, dat.true, k, beta, replications,
                              penalties, str.em, parallel, nodes.ordered){
  #If data has been provided, learn parameters from data
  #Otherwise, assume
  if(!is.null(dat.true) && class(dag.true) == "bn"){
    dag.true.fit <- bn.fit(dag.true, dat.true)
  }
  else{
    if(!class(dag.true) == "bn.fit")
    {stop("Provide either data or a bn.fit object")}
    dag.true.fit <- dag.true
    dag.true <- bn.net(dag.true)
  }
  
  generateData <- function(dag.true.fit, k, beta){
    #Generate complete dataset
    no.params <- nparams(dag.true.fit)
    n <- k*no.params
    dat.genr <- rbn(dag.true.fit, n)
    
    #Randomly remove observations
    if(beta > 0){
      p <- nnodes(dag.true.fit)
      miss <- as.data.frame(matrix(rbinom(n*p, 1, beta), nrow = n, ncol = p))
      dat.genr[miss == 1] <- NA
    }
  }
  
  #Generate datasets in list of length 'replications'
  datasets <- replicate(replications, generateData, simplify = FALSE, k = k,
                        dag.true.fit = dag.true.fit, beta = beta)
  
  result.list <- lapply(datasets, experiment.layer2, dag.true = dag.true,
                        penalties = penalties, str.em = str.em,
                        nodes.ordered = nodes.ordered)
  
  require(data.table)
  result.df <- rbindlist(result.list)
  result.df$k <- k
  result.df$beta <- beta
  return(result.df)
}

experiment.layer4 <- function(dag.true, dat.true, dag.name, k.vec, beta.vec, 
                              replications, penalties, str.em, 
                              parallel, nodes.ordered){
  input <- expand.grid(k.vec, beta.vec)
  names(input) <- c("k", "beta")
  
  result.list <- mapply(experiment.layer3, k=input$k, beta=input$beta, 
                        MoreArgs = list(dag.true = dag.true, 
                                        dat.true = dat.true, 
                                        replications = replications, 
                                        penalties = penalties,
                                        str.em = str.em, 
                                        parallel = parallel, 
                                        nodes.ordered = nodes.ordered),
                        SIMPLIFY = F)
  
  require(data.table)
  result.df <- rbindlist(result.list)
  result.df$dag <- dag.name
  return(result.df)
}

#Use mapply to simultaneously pass DAG, name and data
experiment.layer5 <- function(dag.vec, dag.names, dat.vec, k.vec, beta.vec, 
                              replications, penalties, str.em = FALSE, 
                              parallel = TRUE, nodes.ordered = NULL){
  result.list <- mapply(experiment.layer4, dag.true = dag.vec, 
                        dat.true = dat.vec, dag.name = dag.names,
                        MoreArgs = list(dag.true = dag.true, 
                                        dat.true = dat.true, 
                                        replications = replications, 
                                        penalties = penalties,
                                        str.em = str.em, 
                                        parallel = parallel, 
                                        nodes.ordered = nodes.ordered),
                        SIMPLIFY = F)
  
  require(data.table)
  result.df <- rbindlist(result.list)
  result.df$dag <- dag.name
  return(result.df)
  
}