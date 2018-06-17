library(bnlearn)

#Computes contingency table for a discrete node and its parents in a given DAG
computeCount.node <- function(node, dag, dat){
  currentParents <- parents(dag, node)
  familyData <- na.omit(dat[,c(node, currentParents)])
  counts <- table(familyData)
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


computeNAL.discrete <- function(node, dag, dat){
  
  #Check that none of the discrete node's parents are continuous. If any parent
  #is continuous, return -Inf to indicate that the DAG is invalid
  parentsContinuous <- sapply(parents(dag, node), function(parent){
    is.numeric(dat[[parent]])})
  if(any(parentsContinuous)){return(-Inf)}
  
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

#TODO: replace computeScore.node by new version wherever it appears
computeScore.node <- function(node, dag, dat, penalty){
  
  if(is.factor(dat[[node]])){
    #TODO: adapt computeNAL.discrete:
    # - new signature, must compute count internally
    # - return df 
    # - must check that no parent is continuous (return -Inf logl and 0 df)
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
  score <- mean(score.nodes, na.rm = TRUE)*nnodes(dag)
  
  return(score)
}


