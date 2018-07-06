library(bnlearn)


#' Computes sufficient statistics for a node in a Discrete Bayesian Network
#'
#' Makes contingency table for a multinomial node and its multinomial parents in
#' a DAG
#'
#' @param node string, name of node
#' @param dag DAG determining parental sets of nodes
#' @param dat dataframe with a factor column for each node in network
#'
#' @return Contingency table. First dimension corresponds to levels of node,
#'   further dimensions correspond to levels of parents
#'
#' @export
#'
#' @examples
computeCount.node <- function(node, dag, dat){
  currentParents <- parents(dag, node)
  familydata <- na.omit(dat[,c(node, currentParents)])
  counts <- table(familydata)
  return(counts)
}

#' Computes conditional counts per node in a discrete Bayesian network
#'
#'
#' @inheritParams computeCount.node
#' @param cl Optional parallel cluster. Must have bnlearn loaded and each
#'   cluster must have access to dataframe 'dat'
#'
#' @return List of arrays. First dimension in each array corresponds to
#'   different levels of node, all other dimensions correspond to different
#'   levels of parents
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


#' Compute node-average log-likelihood for a CGN and a given dataset
#'
#' Allows computation of likelihood over a Bayesian network in the presence of
#' missing data
#'
#' @inheritParams computeCount.node
#' @param dat dataframe to be used. Columns can be factors or continuous
#'   variables
#'
#' @return Node-average likelihood per observation (logLik will return the sum
#'   of likelihoods over observations)
#'
#' @export
#'
#' @examples
computeNAL <- function(dag, dat){
  node.names <- names(dat)
  logl.nodes <- sapply(node.names, function(node)
  {computeNAL.node(node, dag, dat)$logl})
  logl <- sum(logl.nodes)
  return(logl)
}


#' Determines which parents of a given node are continuous RV
#'
#' @param node name of node under inspection
#' @param dag DAG determining parental sets
#' @param dat dataset
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

#' Computes node-average likelihood for a single node
#'
#' @inheritParams parentsContinuous 
#'
#' @return NAL for given node
#' 
#' @export
#'
#' @examples
computeNAL.node <- function(node, dag, dat){
  if(is.factor(dat[[node]])){
    return(computeNAL.discrete(node, dag, dat))
  }
  else{
    return(computeNAL.gaussian(node, dag, dat))
  }
}


#' Computes NAL for discrete nodes in CGNs
#'
#' @inheritParams computeNAL.node
#'
#' @return NAL for given node
#' @export
#'
#' @examples
computeNAL.discrete <- function(node, dag, dat){
  
  #Check that none of the discrete node's parents are continuous. If any parent
  #is continuous, return -Inf to indicate that the DAG is invalid
  if(any(parentsContinuous(node, dag, dat)))
    {return(list(logl = -Inf, df = 0))}
  
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
    # if(any(parent.counts < 1)){
    #   node.names <- names(dimnames(counts))
    #   warning(c("No valid observations for ", node.names[1] , 
    #             " for some configuration of parents\n", 
    #             encodeString(node.names[-1], quote = " ") ,
    #             ": result is NaN"))
    #   return(list(logl = NaN, df = 0))
    # }
    
    #Compute MLEs
    mle.conditional <- sweep(counts, parentIndices, parent.counts, "/")
    
    #If there were no observations for some parental configuration, replace
    #these by a uniform distribution (maximum entropy, lowest possible score)
    mle.conditional[!is.finite(mle.conditional)] <- 1/dim(counts)[1]
    
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


#' Computes NAL for continuous nodes in CGNs
#'
#' @inheritParams computeNAL.node
#'
#' @return NAL for given node
#' @export
#'
#' @examples
computeNAL.gaussian <- function(node, dag, dat){
  #Obtain dataframe with only node and parents, omit rows with missing entries
  currentParents <- parents(dag, node)
  familyData <- na.omit(dat[,c(node, currentParents)])
  
  #If there are no parents, use unconditional MLEs 
  if(length(currentParents) == 0){
    logl <- mean(dnorm(familyData, mean(familyData), sd(familyData), 
                       log = TRUE))  
    no.params <- 2
    
  } else {
    n <- nrow(familyData)
    
    #Determine which parents are continuous and which are discrete
    whichParentsContinuous <- parentsContinuous(node, dag, dat)
    continuousParents <- currentParents[whichParentsContinuous]
    discreteParents <- currentParents[!whichParentsContinuous]
    
    #If there are no continuous parents, use conditional MLEs for each parental
    #configuration
    if(length(continuousParents) == 0){
      #Split into datasets for each configuration of discrete parents
      discreteParents <- currentParents
      dataConfigs <- split(familyData, familyData[,discreteParents])
      
      #Compute average log-likelihood for each configuration, multiplied by
      #number of observations for each configuration
      logl <- sum(sapply(dataConfigs, function(dataConfig){
        x <- dataConfig[, node]
        return(sum(dnorm(x, mean(x), sd(x), log = TRUE)))}))/n
      
      no.params <- length(dataConfigs)*2
      
      #If there are continuous parents, regress current node on continuous parents
    } else {
      f <- as.formula(paste(node, "~", paste(continuousParents, collapse = " + ")))
      
      #If there are no discrete parents, do a single regression
      if(length(discreteParents) == 0){
        logl <- as.vector(logLik(lm(f, data = familyData)))/n
        no.params <- length(continuousParents) + 2
        
        #If there are discrete and continuous parents, do a regression for each
        #parental configuration    
      } else {
        dataConfigs <- split(familyData, familyData[,discreteParents])
        
        #Log-likelihood is sum of log-likelihoods for each observation. Hence it
        #is already proportional to the counts for each configuration
        logl <- sum(sapply(dataConfigs, function(dataConfig){
          return(as.vector(logLik(lm(f, data = dataConfig)))) }))/n
        
        no.params <- length(dataConfigs)*(length(continuousParents)+2)
      }
    }
  }
  
  return(list(logl = logl, df = no.params))
}


#' Compute penalized likelihood for a node in a CGN
#'
#' @inheritParams computeNAL.node
#' @param penalty string determining how fast the complexity penalty grows in n.
#'   Options are "bic" (log(n)/n), "aic" (1/n) or any number larger than 0, in
#'   which case the penalty grows as n^(-number)
#'
#' @return score for given node
#' @export
#'
#' @examples
computeScore.node <- function(node, dag, dat, penalty, no.nodes = NULL){
  
  NAL <- computeNAL.node(node, dag, dat)
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
  
  if(is.null(no.nodes)){no.nodes <- nnodes(dag)}
  n <- nrow(dat)
  
  if(bic){penaltyFactor <- 0.5*log(n)/n}
  else if (aic){penaltyFactor <- 1/n}
  else {penaltyFactor <- 1/no.nodes * n^(-alpha)}
  
  nodeScore <-  logl - penaltyFactor*no.params
  return(nodeScore)
}


#' Compute score for CGN
#'
#' @inheritParams computeScore.node
#' @param cl optional parallel cluster. Should have relevant functions, packages
#'   and data pre-loaded
#'
#' @return score for full network on given dataset
#' @export
#'
#' @examples
computeScore <- function(dag, dat, penalty, cl = NULL){
  
  node.names <- names(dat)
  score.nodes <- sapply(node.names, computeScore.node, dag = dag, dat = dat,
                        penalty = penalty)
  
  score <- mean(score.nodes, na.rm = TRUE)*nnodes(dag)
  
  return(score)
}


