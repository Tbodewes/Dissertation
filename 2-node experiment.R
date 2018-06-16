library(parallel)
source("Functions.R") #loads bnlearn

#Two node experiment from Balov:
twoNodeExperiment.instance <- function(runs, n, beta, penalty, debug = FALSE){
  singleRun <- function(n, beta, penalty, debug = FALSE){
    x1 <- factor(rbinom(n, 1, 0.6))
    x2 <- factor(rbinom(n, 1, 0.7))
    
    z <- rbinom(n, 1, 1-beta)
    x1[z == 1] <- NA
    
    df <- data.frame(X1 = x1, X2 = x2)
    
    g1 <- empty.graph(c("X1", "X2"))
    g2 <- set.arc(g1, from = "X1", to = "X2" )
    
    score1 <- computeScore(g1, df, penalty)
    score2 <- computeScore(g2, df, penalty)
    
    if(debug){cat("Score 1:", score1," Score 2:", score2, "\n")}
    return(score1 < score2)
  }
  
  cat("n:", n, "beta:", beta, "penalty:", penalty, "\n")
  outcomes <- replicate(runs, singleRun(n, beta, penalty, debug))
  return(mean(outcomes))
}

twoNodeExperiment <- function(runs, n.vec, beta.vec, penalty.vec, 
                              parallel = FALSE){
  input <- expand.grid(beta.vec, penalty.vec, n.vec, stringsAsFactors = FALSE)
  names(input) <- c("beta", "penalty", "n")
  input <- input[c("n", "penalty", "beta")]
  
  if(parallel){
    
    #Set up functions and input format to allow use of parLapply (as)
    twoNodeExperiment.par <- function(inputVec, runs){
      n <- as.numeric(inputVec[1])
      penalty <- inputVec[[2]]
      beta <- as.numeric(inputVec[3])
      
      return(twoNodeExperiment.instance(runs, n, beta, penalty))
    }
    
    input.list <- split(input, seq(nrow(input)))
    
    no.cores <- detectCores() - 1
    cl <- makeCluster(no.cores)
    on.exit(stopCluster(cl))
    
    invisible(clusterEvalQ(cl, library(bnlearn)))
    
    functionsToLoad <- c("twoNodeExperiment.instance", "computeScore",
                         "computeScore.node", "computeNAL.node",
                         "computeCounts")
    invisible(clusterExport(cl, functionsToLoad))
    
    output <- parSapplyLB(cl, input.list, twoNodeExperiment.par, runs = runs)
  }
  else{
    output <- mapply(twoNodeExperiment.instance, n = input$n, beta = input$beta,
                     penalty = input$penalty, runs = runs)
  }
  
  result <- cbind(input, output)
  
  return(result)
}

runs <- 1000
n.vec <- 10^(2:5)
beta.vec <- c(1, 0.99, 0.95, 0.9, 0.75)
penalty.vec <- c(seq(0.2, 0.8, by = 0.1), "bic", "aic")

result <- twoNodeExperiment(runs, n.vec, beta.vec, penalty.vec, parallel = TRUE)

save("result", file = "2-node experiment.Rdata")
