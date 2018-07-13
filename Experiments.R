#' ** Setup computational experiment **

library(data.table)
library(dplyr)

#' Experiment is done for all combinations of 1) input DAG, 2) k (relative
#' sample size), 3) beta (probability for an entry to be missing) and 4)
#' penalty. Each of these combinations is also run on multiple generated
#' datasets. We use a layered structure to iterate over these in a natural
#' manner. Each layer collects the results from the previous layer
#'
#' Layer 1: collects output for a single given dataset, single penalty
#'
#' Layer 2: single dataset, all penalties
#'
#' Layer 3: multiple datasets for given k and beta, all penalties
#'
#' Layer 4: multiple k and beta, all penalties
#'
#' Layer full: multiple DAGs, multiple k and beta, all penalties
#'
#' Experiment can use either local search to fit for unknown node order, or the
#' optimal algorithm for known node order. In the latter case the maximum number
#' of parents is 3. Can also optionally run EM on the same datasets to compare
#' its performance with NAL optimization.


#' First (inner) layer of computational experiment
#'
#' For a single dataset and penalty, fit DAG using NAL optimization and
#' optionally using structural EM. Extract number of queries and SHD for both.
#' Can use local search or ordered fitting for NAL optimization. In case of the
#' former, each node in the true graph should have at most 3 parents
#'
#' @param dat generated dataset to use for fitting
#' @param dag.true dag to compare learned graph with
#' @param penalty penalty to be used when fitting. Positive number for a penalty
#'   of the form n^(-number), or 'bic' or 'aic'
#' @param str.em TRUE if structural EM should also be run for comparison
#' @param ordered TRUE if ordered fitting should be done (str.em should be false
#'   in this case, as no experiment is done to compare ordered fitting with EM)
#'
#' @return list containing SHD, number of queries and total running time for NAL
#'   and for EM, in that order
#'
#' @export
#'
#' @examples
experiment.layer1 <- function(dat, dag.true, penalty, str.em, ordered){
  
  #Fit DAG either using ordered or general fitting
  if(ordered){
    nodes.ordered <- node.ordering(dag.true)  
    time <- system.time(dag.fitted <- 
                          fit.dag.ordered(nodes.ordered,
                                          max.parents = 3, dat, penalty, 
                                          parallel = FALSE))[3]
    if(str.em){stop("Structural EM and ordered fitting should not be run in 
                     the same experiment")}
  }
  else{
    time <- system.time(dag.fitted <- fit.dag(dat, penalty, parallel = FALSE,
                                              tabuSteps = 10))[3]
  }
  
  #Record SHD and number of queries
  no.arcs <- narcs(dag.true)
  distance <- shd(dag.fitted, dag.true)/no.arcs
  no.queries <- dag.fitted$learning$ntests
  
  #Run structural EM and extract performance measures if necessary. Otherwise
  #set performance measures to NA
  if(str.em){
    time.em <- system.time(dag.em <- structural.em(dat, maximize = "tabu", 
                                                   fit.args = list(replace.unidentifiable = TRUE)))[3]
    distance.em <- shd(dag.em, dag.true)/no.arcs
    no.queries.em <- dag.em$learning$ntests
  }
  else{distance.em <- NA; no.queries.em <- NA; time.em <- NA}
  
  return(list(distance, no.queries, time, distance.em, no.queries.em, time.em))
}

#' Second layer of computational experiment
#'
#' Run layer 1 for different penalties
#'
#' @inheritParams experiment.layer1
#' @param penalties vector of penalties to run layer 1 for
#'
#' @return dataframe with penalties on columns and output for each penalty on
#'   rows
#'
#' @export
#'
#' @examples
experiment.layer2 <- function(dat, dag.true, penalties, str.em, ordered){
  
  if(str.em && length(penalties) > 1)
  {warning("Experiment with structural EM should be run for a single penalty.
           Program will now run EM unnecessarily for different penalties")}
  
  #Run layer 1 for each penalty and collect results in appropriately named df 
  result.list <- lapply(penalties, experiment.layer1, dat = dat, 
                        dag.true = dag.true, str.em = str.em,
                        ordered = ordered)
  result.df <- data.frame(penalties, matrix(unlist(result.list), byrow = TRUE,
                                            nrow = length(penalties)))
  
  names(result.df) <- c("penalty", "SHD.NAL", "Q.NAL", "T.NAL", "SHD.EM",
                        "Q.EM", "T.EM")
  return(result.df)
}


#' Third layer of computational experiments
#'
#' Runs layer 2 for a given number of replications, each on independently
#' generated datasets.
#'
#' @param dag.true bn or bn.fit object to be used to generate data and to
#'   compare learned DAG to. If a bn object is provided, one must also provide
#'   dat.true, so parameter learning can be done internally
#' @param dat.true optional dataframe containing true data on which synthetic
#'   data should be based
#' @param k sample sizes divided by number of parameters in dag.true
#' @param beta probability for any entry in the generated dataframes to be
#'   missing
#' @param replications number of datasets to generate
#' @param cl optional parallel cluster to be used for parallelization
#' @inheritParams experiment.layer2
#'
#' @return long dataframe with unaggrated results over replications for each
#'   penalty
#'
#' @export
#'
#' @examples
experiment.layer3 <- function(dag.true, dat.true = NULL, k, beta, 
                              replications, penalties, str.em,
                              cl = NULL, ordered){
  cat("k: ", k, "; beta: ", beta, "\n", sep ="")
  
  #If data has been provided, learn parameters from data
  #Otherwise, a bn.fit object must have been provided
  if(!is.null(dat.true) && class(dag.true)[1] == "bn"){
    dag.true.fit <- bn.fit(dag.true, dat.true)
  }
  else{
    if(!class(dag.true)[1] == "bn.fit")
    {stop("Provide either data or a bn.fit object")}
    dag.true.fit <- dag.true
    dag.true <- bn.net(dag.true)
  }
  
  #Generate k*no.params(dag.true.fit) observations, each of which is MCAR with
  #constant probability beta
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
    
    return(dat.genr)
  }
  
  #Generate datasets in list of length 'replications'
  datasets <- replicate(replications, simplify = FALSE,
                        expr = generateData(dag.true.fit = dag.true.fit,
                                            k = k, beta = beta))
  
  #Run experiment for each dataset, in parallel if a cluster has been provided
  if(!is.null(cl)){
    result.list <- parLapply(cl, datasets, experiment.layer2, dag.true = dag.true,
                             penalties = penalties, str.em = str.em,
                             ordered = ordered)
  }
  else{
    result.list <- lapply(datasets, experiment.layer2, dag.true = dag.true,
                          penalties = penalties, str.em = str.em,
                          ordered = ordered) 
  }
  
  #Process results into single dataframe
  require(data.table)
  result.df <- rbindlist(result.list)
  result.df$k <- k
  result.df$beta <- beta
  return(result.df)
}

#' Fourth layer of computational experiments
#'
#' @param dag.name character string marking DAG
#' @param k.vec numeric vector of values for k to use in layer 3
#' @param beta.vec numeric vector of values for beta to use in layer 3
#' @inheritParams experiment.layer3
#'
#' @return long dataframe with output from layer 3 for each combination of k and
#'   beta, marked by a column k and a column beta
#' @export
#'
#' @examples
experiment.layer4 <- function(dag.true, dat.true = NULL, dag.name, k.vec,  
                              beta.vec, replications, penalties, str.em, 
                              cl, ordered){
  cat("DAG:", dag.name, "\n")
  input <- expand.grid(k.vec, beta.vec)
  names(input) <- c("k", "beta")
  
  result.list <- mapply(experiment.layer3, k=input$k, beta=input$beta, 
                        MoreArgs = list(dag.true = dag.true, 
                                        dat.true = dat.true, 
                                        replications = replications, 
                                        penalties = penalties,
                                        str.em = str.em, 
                                        cl = cl, 
                                        ordered = ordered),
                        SIMPLIFY = F)
  
  require(data.table)
  result.df <- rbindlist(result.list)
  result.df$dag <- dag.name
  return(result.df)
}


#' Full computational experiment
#'
#' @param dag.vec list of bn or bn.fit objects to run experiment for
#' @param dag.names labels of DAGs in dag.vec (for use in dataframe)
#' @param dat.vec list of dataframes to be used for parameter learning in DAGs,
#'   such that datasets can be generated. Any entry can be NULL if the
#'   corresponding entry in dag.vec is a bn.fit object
#' @inheritParams experiment.layer4
#'
#' @return long dataframe with for each combination of 1) DAG, 2) k, 3) beta, 4)
#'   penalty, replications rows containing output of a single run of the
#'   computational experiment. Rows consist of SHD, number of queries and user
#'   running time for NAL optimization and optionally also for structural EM
#'
#' @export
#'
#' @examples
experiment.full <- function(dag.vec, dag.names, k.vec, beta.vec, dat.vec,
                            replications, penalties, str.em = FALSE, 
                            parallel = TRUE, ordered = FALSE){
  #Prepare parallel cluster if necessary
  if(parallel){
    cl <- makeCluster(detectCores() - 1)
    on.exit(stopCluster(cl))
    
    invisible(clusterEvalQ(cl, library(bnlearn)))
    invisible(clusterEvalQ(cl, source("Scoring.R")))
    invisible(clusterEvalQ(cl, source("Fit DAG.R")))
    invisible(clusterEvalQ(cl, source("Experiments.R")))
  }
  else{
    cl <- NULL
  }
  
  #Run layer 4 for each DAG
  result.list <- mapply(experiment.layer4, dag.true = dag.vec, 
                        dat.true = dat.vec, dag.name = dag.names,
                        MoreArgs = list(k.vec = k.vec,
                                        beta.vec = beta.vec,
                                        replications = replications, 
                                        penalties = penalties,
                                        str.em = str.em, 
                                        cl = cl, 
                                        ordered = ordered),
                        SIMPLIFY = F)
  
  #Aggregate results into single dataframe and format nicely
  require(data.table); require(dplyr)
  result.df <- rbindlist(result.list)
  setcolorder(result.df, c("dag", "k", "beta", "penalty"))
  result.df <- result.df[with(result.df, order(dag, k, beta, penalty)),]
  
  return(result.df)
}

#' Prune network to restrict size of parental sets
#'
#' Keeps only a specified number of parents for each node. The first max.parents
#' nodes from a vector sorted by node names of the parents are kept. Nodes that
#' do not have more than max.parents parents are left unaffected
#'
#' @param dag bn object to prune
#' @param max.parents maximum number of parents per node in pruned graph
#'
#' @return pruned graph
#' 
#' @export
#'
#' @examples
pruneNetwork <- function(dag, max.parents){
  for(node in nodes(dag)){
    par <- parents(dag, node)
    #If too many parents, keep the first ones based on alphabetical order
    if(length(par) > max.parents){
      parents(dag, node) <- sort(par)[1:max.parents]
    }
  }
  return(dag)
}

#' Prune bn.fit objects
#'
#' Wrapper for pruneNetwork that prunes bn.fit objects by generating a dataset
#' from the original network, then pruning the network and then learning
#' parameters for the pruned network using the generated dataset
#'
#' @param dag.fit bn.fit object to prune
#' @param max.parents maximum number of parents
#' @param k relative size of synthetic dataset. Actual size is k*number of
#'   parameters in network
#'
#' @return pruned bn.fit object
#' 
#' @export
#'
#' @examples
pruneFit <- function(dag.fit, max.parents, k = 100){
  dat <- rbn(dag.fit, k*nparams(dag.fit))
  dag.pruned <- pruneNetwork(bn.net(dag.fit, max.parents), max.parents)
  return(bn.fit(dag.pruned, dat))
}