library(bnlearn)

#' ** Functions for computing score **

#' Computes penalized Node-Average Likelihood (NAL) of a Conditional Gaussian
#' Bayesian Network (CGBN) on a given dataset
#'
#' Computes the NAL for discrete and continuous nodes in a given DAG and
#' substracts from it a complexity penalty. NAL handles incomplete data by using
#' locally complete data. That is, each term in the log-likelihood only depends
#' on data for one node and its parents. Each term can thus be computed using
#' observations for which this data is complete. This gives the NAL.
#'
#' Penalization takes the form of (1/nnodes(dag)) n^(-alpha) * nparams(dag) or
#' of BIC or AIC
#'
#' @param dag bn object specifying parental set of each node
#' @param dat dataframe with data for each node
#' @param penalty string determining penalization. Should be "bic", "aic" or a
#'   number larger than zero
#'
#' @return score of DAG
#'
#' @export
computeScore <- function(dag, dat, penalty, cl = NULL){
  
  node.names <- names(dat)
  score.nodes <- sapply(node.names, computeScore.node, dag = dag, dat = dat,
                        penalty = penalty)
  
  score <- mean(score.nodes, na.rm = TRUE)*nnodes(dag)
  
  return(score)
}

#' Compute penalized NAL for a discrete or continuous node in a CGN
#'
#' @param node name of node to be scored
#' @inheritParams computeScore
#'
#' @return Score for given node
#' 
#' @export
computeScore.node <- function(node, dag, dat, penalty, no.nodes = NULL){
  
  #Compute NAL and retrieve number of degrees of freedom
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
  
  #Compute and return score
  nodeScore <-  logl - penaltyFactor*no.params
  return(nodeScore)
}

#' ** Functions for computing NAL **

#' Computes total NAL of a CGBN
#'
#' @inheritParams computeScore
#'
#' @return NAL
#' 
#' @note NAL is averaged over observations, i.e. for complete data, 
#' computeNAL(dag, dat) = logLik(dag, dat)/nrow(dat) 
#' 
#' @export
computeNAL <- function(dag, dat){
  node.names <- names(dat)
  logl.nodes <- sapply(node.names, function(node)
  {computeNAL.node(node, dag, dat)$logl})
  logl <- sum(logl.nodes)
  return(logl)
}


#' Computes NAL for a discrete or continuous node in a CGBN
#'
#' @inheritParams computeScore.node
#'
#' @return NAL and number of fitted parameters for a given node
#'
#' @export
computeNAL.node <- function(node, dag, dat){
  if(!node %in% nodes(dag) || !node %in% names(dat))
  {stop("Invalid node name: does not appear in either DAG or data")}
  
  if(is.factor(dat[[node]])){
    return(computeNAL.discrete(node, dag, dat))
  }
  else{
    return(computeNAL.gaussian(node, dag, dat))
  }
}


#' * Functions for scoring discrete nodes *

#' Computes NAL for a discrete node in a CGBN
#' 
#' NAL for discrete nodes is based on fitting a categorical distribution for
#' each configuration of discrete parents
#'
#' @inheritParams computeScore.node
#'
#' @return NAL and number of fitted parameters for given node
#'
#' @note returns -Inf if the node has continuous parents, as this is not allowed
#'   in a CGBN
#'
#' @export
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


#' Computes contingency table for a node and its parents using locally complete
#' data
#'
#' @inheritParams computeScore.node
#'
#' @return Contingency table. First dimension corresponds to levels of node,
#'   further dimensions correspond to levels of parents
#'
#' @export
computeCount.node <- function(node, dag, dat){
  currentParents <- parents(dag, node)
  familydata <- na.omit(dat[,c(node, currentParents)])
  counts <- table(familydata)
  return(counts)
}


#' Determines which parents of a given node are continuous RV
#'
#' @inheritParams computeScore.node
#'
#' @return logical vector with an entry for each parent. TRUE if parent is
#'   numeric in dat, FALSE if parent is a factor
#'   
#' @export
parentsContinuous <- function(node, dag, dat){
  return(sapply(parents(dag, node), function(parent){is.numeric(dat[[parent]])}))
}


#' Computes NAL for continuous nodes in a CGBN
#'
#' NAL for continuous nodes is based on fitting a linear regression of with the
#' continuous parents as regressors, for each configuration of discrete parents
#'
#' @inheritParams computeScore.node
#'
#' @return NAL and number of fitted parameters for given node
#'
#' @export
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
      #number of observations for each configuration. Note that sum of densities
      #is already on the order of the number of observations
      configLogL <- function(dataConfig){
        x <- dataConfig[, node]
        if(nrow(dataConfig) > 10)
        {return(sum(dnorm(x, mean(x), sd(x), log = TRUE)))}
        else
        {return(0)}
      }
      logl <- sum(sapply(dataConfigs, configLogL))/n
      
      no.params <- length(dataConfigs)*2
      
      #If there are continuous parents, regress current node on continuous parents
    } else {
      f <- as.formula(paste(node, "~", paste(continuousParents, collapse = " + ")))
      
      #If there are no discrete parents, do a single regression
      #If there is no data available in this case (all observations for this 
      #parental set contain at least one missing value), set logL to negative
      #infinity. Parental set is deemed too complex in this case
      if(length(discreteParents) == 0){
        if(n > 5*length(continuousParents)){
          logl <- as.vector(logLik(lm(f, data = familyData)))/n
        } else {
          logl <- -Inf
        }
        no.params <- length(continuousParents) + 2
        
        #If there are discrete and continuous parents, do a regression for each
        #parental configuration    
      } else {
        dataConfigs <- split(familyData, familyData[,discreteParents])
        
        #Log-likelihood is sum of log-likelihoods for each observation. Hence it
        #is already proportional to the counts for each configuration
        configLogL <- function(dataConfig, f){
          x <- dataConfig[, node]
          if(nrow(dataConfig) > 5*length(continuousParents)){
            return(as.vector(logLik(lm(f, data = dataConfig))))
          } else {
            return(0)
          }
        }
        logl <- sum(sapply(dataConfigs, configLogL, f = f))/n
        
        no.params <- length(dataConfigs)*(length(continuousParents)+2)
      }
    }
  }
  
  return(list(logl = logl, df = no.params))
}

