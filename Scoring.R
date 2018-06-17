library(bnlearn)

#Computes contingency table for a discrete node and its parents in a given DAG
computeCount.node <- function(node, dag, dat){
  currentParents <- parents(dag, node)
  familydata <- na.omit(dat[,c(node, currentParents)])
  counts <- table(familydata)
  return(counts)
}

#' Computes conditional counts per node in a discrete Bayesian network
#' 
#'
#' @param dag Directed acyclic graph determining which parents to use for each 
#' node
#' @param dat Dataframe with a factor column for each variable
#' @param cl Optional parallel cluster. Must have bnlearn loaded
#' and each cluster must have access to dataframe 'dat' 
#'
#' @return List of arrays. First dimension in each array corresponds to different
#' levels of node, all other dimensions correspond to different levels of parents
#' @export
#'
#' @examples
computeCounts <- function(dag, dat, cl = NULL){
  node.names <- names(dat)
  
  #Paralellize if a cluster has been supplied, otherwise execute sequentially
  if(is.null(cl))
  {counts.list <- lapply(node.names, computeCount.node, dag = dag, dat = dat)}
  else
  {counts.list <- parLapply(cl, node.names, computeCount.node, dag = dag, 
                            dat = dat)}
  
  names(counts.list) <- node.names
  return(counts.list)
}


#' Compute node-average log-likelihood for a DBN and a dataset 
#'
#' @param dag directed acyclic graph to be used
#' @param df dataframe of factors to be used
#'
#' @return
#' @export
#'
#' @examples
computeNAL <- function(dag, dat){
  node.names <- names(dat)
  logl.nodes <- sapply(node.names, function(node)
  {computeNAL.discrete(node, dag, dat)$logl})
  logl <- sum(logl.nodes)
  
  return(logl)
}


#' Determines which parents of a given node are continuous RV
#'
#' @param node string giving node whose parents should be inspected
#' @param dag DAG determining parental set
#' @param dat dataset informing whether a variable is numeric or a factor
#'
#' @return logical vector with an entry for each parent. TRUE if parent is
#'   numeric in dat, FALSE if parent is a factor
#' @export
#'
#' @examples
parentsContinuous <- function(node, dag, dat){
  return(sapply(parents(dag, node), function(parent){
    is.numeric(dat[[parent]])}))
}

computeNAL.discrete <- function(node, dag, dat){
  
  #Check that none of the discrete node's parents are continuous. If any parent
  #is continuous, return -Inf to indicate that the DAG is invalid
  if(any(parentsContinuous(node, dag, dat))){return(-Inf)}
  
  counts <- computeCount.node(node, dag, dat)
  
  #Distinguish between nodes with and without parents. In the first case the 
  #counts are a vector, in the latter an array. These are handled differently
  #by inbuilt R functions such as colSums (which requires an array)
  familySize <- length(dim(counts))
  if(familySize == 1){
    total <- sum(counts)
    
    if(total < 1)
    {warning(c("No valid observations, result is NaN"))}
    
    mle <- counts/total
    logl <- sum(mle*log(mle))
    
    #Compute number of parameters (degrees of freedom of node log-likelihood)
    n <- sum(counts)
    no.params <- length(counts) - 1
  }
  else{
    parentIndices <- 2:familySize
    parent.counts <-  colSums(counts)
    
    #Give a warning if there are no observations for some parental configuration
    if(any(parent.counts < 1)){
      node.names <- names(dimnames(counts))
      warning(c("No valid observations for ", node.names[1] , 
                " for some configuration of parents\n", 
                encodeString(node.names[-1], quote = " ") ,
                ": result is NaN"))
      return(NaN)
    }
    mle.conditional <- sweep(counts, parentIndices, parent.counts, "/")
    
    #Replace underflown MLEs by a small positive number for numerical stability
    mle.conditional[mle.conditional < 1E-10] <- 1E-10
    term.conditional <- colSums(mle.conditional*log(mle.conditional))
    
    #Distinguish between families sized one and larger to make use of
    #rowSums, which only handles arrays of dimension two or more
    if(familySize > 2){ 
      total.counts <- sum(rowSums(parent.counts))
      mle.parents <- parent.counts/total.counts
      logl <- sum(rowSums(mle.parents*term.conditional))
    }
    else{
      total.counts <- sum(parent.counts) 
      mle.parents <- parent.counts/total.counts
      logl <- sum(mle.parents*term.conditional)
    }
    
    #Compute number of parameters
    n <- sum(rowSums(counts))
    nlevels <- dim(counts)
    no.params <- (nlevels[1]-1)*prod(nlevels[-1])
  }
  
  return(list(logl = logl, df = no.params))
}


computeNAL.gaussian <- function(node, dag, dat){
  #Obtain dataframe with only node and parents, omit rows with missing entries
  currentParents <- parents(dag, node)
  familyData <- na.omit(dat[,c(node, currentParents)])
 
  #Deal with case of no parents by filling in MLEs in normal likelihood
  if(length(currentParents) < 1){
    n <- length(familyData)
    logl <- sum(dnorm(familyData, mean(familyData), sd(familyData), log = TRUE))
    return(list(logl = logl/n, df = 2))
  }
  
  n <- nrow(familyData)
  
  #Determine which parents are continuous and which are discrete
  whichParentsContinuous <- parentsContinuous(node, dag, dat)
  continuousParents <- currentParents[whichParentsContinuous]
  discreteParents <- currentParents[!whichParentsContinuous]
  
  #Regress node on continuous parents
  f <- as.formula(paste(node, "~", paste(continuousParents, collapse = " + ")))
  
  #If there are discrete parents, do separate regression for each configuration
  if(length(discreteParents) > 0){
    dataConfigs <- split(familyData, familyData[,discreteParents])
    
    #Log-likelihood is sum of log-likelihoods for each observation. Hence it is
    #already proportional to the counts for each configuration
    weighted.logl <- sum(sapply(dataConfigs, function(dataConfig){
      return(as.vector(logLik(lm(f, data = dataConfig))))
    }))
    no.params <- length(dataConfigs)*(length(continuousParents)+2)
  } else{
    weighted.logl <- as.vector(logLik(lm(f, data = familyData)))
    no.params <- length(continuousParents) + 2
  }
  
  return(list(logl = weighted.logl/n, df = no.params))
}

computeScore.node <- function(node, dag, dat, penalty){
  
  if(is.factor(dat[[node]])){
    NAL <- computeNAL.discrete(node, dag, dat)
  }
  else{
    NAL <- computeNAL.gaussian(node, dag, dat)
  }
  logl <- NAL$logl
  no.params <- NAL$df
  
  #Configure and compute penalty
  bic = F
  aic = F
  alpha = NULL
  
  if(penalty == "bic"){bic = TRUE}
  else if(penalty == "aic"){aic = TRUE}
  else if(!is.na(is.numeric(penalty))){alpha = as.numeric(penalty)}
  else {stop("Penalty must be 'bic', 'aic' or a number")}
  
  no.nodes <- nnodes(dag)
  n <- nrow(dat)
  
  if(bic){penaltyFactor <- 0.5*log(n)/n}
  else if (aic){penaltyFactor <- 1/n}
  else {penaltyFactor <- 1/no.nodes * n^(-alpha)}

  nodeScore <-  logl - penaltyFactor*no.params
  return(nodeScore)
}


computeScore <- function(dag, dat, penalty, cl = NULL){
  
  node.names <- names(dat)
  score.nodes <- sapply(node.names, computeScore.node, dag = dag, dat = dat,
                        penalty = penalty)
  browser()
  score <- mean(score.nodes, na.rm = TRUE)*nnodes(dag)
  
  return(score)
}


