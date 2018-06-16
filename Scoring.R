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
computeNAL <- function(dag, df){
  counts.list <- computeCounts(dag, df)
  
  node.names <- names(df)
  logl.nodes <- sapply(counts.list, computeNAL.node)
  logl <- sum(logl.nodes)
  
  return(logl)
}

computeNAL.node <- function(counts){
  familySize <- length(dim(counts))
  if(familySize == 1){
    total <- sum(counts)
    if(total < 1)
    {warning(c("No valid observations, result is NaN"))}
    mle <- counts/total
    return(sum(mle*log(mle)))
  }
  else{
    parentIndices <- 2:familySize
    parent.counts <-  colSums(counts)
    if(any(parent.counts < 1))
    {
      node.names <- names(dimnames(counts))
      warning(c("No valid observations for ", node.names[1] , 
                " for some configuration of parents\n", 
                encodeString(node.names[-1], quote = " ") ,
                ": result is NaN"))}
    mle.conditional <- sweep(counts, parentIndices, parent.counts, "/")
    mle.conditional[mle.conditional < 1E-10] <- 1E-10
    term.conditional <- colSums(mle.conditional*log(mle.conditional))
    
    if(familySize > 2){ 
      total.counts <- sum(rowSums(parent.counts))
      mle.parents <- parent.counts/total.counts
      return(sum(rowSums(mle.parents*term.conditional)))
    }
    else{
      total.counts <- sum(parent.counts) 
      mle.parents <- parent.counts/total.counts
      return(sum(mle.parents*term.conditional))
    }
  }
}


computeScore.node <- function(counts, penalty, no.nodes = 1){
  
  logl.node <- computeNAL.node(counts)
  
  if(length(dim(counts)) > 1){
    n <- sum(rowSums(counts))
    nlevels <- dim(counts)
    no.params <- (nlevels[1]-1)*prod(nlevels[-1])
  }
  else{
    n <- sum(counts)
    no.params <- length(counts) - 1
  }
  
  #Configure and compute penalty
  bic = F
  aic = F
  alpha = NULL
  
  if(penalty == "bic"){bic = TRUE}
  else if(penalty == "aic"){aic = TRUE}
  else if(!is.na(is.numeric(penalty))){alpha = as.numeric(penalty)}
  else {stop("Penalty must be 'bic', 'aic' or a number")}
  
  if(bic){penaltyFactor <- 0.5*log(n)/n}
  else if (aic){penaltyFactor <- 1/n}
  else {penaltyFactor <- 1/no.nodes * n^(-alpha)}
  
  nodeScore <-  logl.node - penaltyFactor*no.params
  return(nodeScore)
}


computeScore <- function(dag, dat, penalty, cl = NULL){
  
  counts.list <- computeCounts(dag, dat, cl)
  
  node.names <- names(dat)
  score.nodes <- sapply(counts.list, computeScore.node, penalty)
  score <- mean(score.nodes, na.rm = TRUE)*nnodes(dag)
  
  return(score)
}


