library(parallel)
library(data.table)

#' Obtains DAG from bootstrap samples and fits parameters using EM
#'
#' Uses model averaging to determine the strength of each arc from a list of
#' bootstrap samples and constructs a DAG which only includes arcs with strength
#' above a given threshold. Then fits the parameters of this DAG using
#' parametric EM.
#'
#' @param samples List of bootstrapped DAGs
#' @param threshold Minimum strength of an arc to be included in DAG
#' @inheritParams em.parametric
#'
#' @return bn.fit object containing averaged DAG with fitted parameters
#' 
#' @export
#'
#' @examples
fit.boot <- function(samples, threshold, dat, particles, obs.perParam,
                     max.iter = 1, parallel = FALSE){
  #Average bootstrapped DAGs
  strength <- custom.strength(samples, names(dat))
  dag <- averaged.network(strength, names(dat), threshold)
  
  #Direct undirected arcs using CGN assumptions 
  undir <- undirected.arcs(dag)
  if(nrow(undir) > 0){
    for(i in 1:nrow(undir)){
      from <- undir[i, 1]
      to <- undir[i, 2]
      if(is.factor(dat[[from]]) && is.numeric(dat[[to]])){
        dag <- set.arc(dag, from, to)
      }
    }
  }
  if(nrow(undirected.arcs(dag)) > 0){stop("Graph is not fully directed")}
  
  #Learn parameters in a single iteration of parametric EM
  dag.fit <- em.parametric(dag, dat, max.iter = max.iter, particles = particles,
                           obs.perParam = obs.perParam, parallel = parallel)$dag
  return(dag.fit)
}


#' Computes in- and out- of sample predictions for a given Bayesian network
#'
#' @param dat Dataframe containing data to use for prediction. Need not be
#'   complete, missing values are marginalized out
#' @param varName Column name of variable from dat to predict
#' @param dag.fit bn.fit object to use for prediction
#' @param train if FALSE (default), an out-of-sample prediction procedure is
#'   used. If TRUE (only available for predicting continuous variables), the
#'   values fitted during training are used
#' @inheritParams em.parametric
#'
#' @return Two column dataframe containing fitted and actual values
#' 
#' @export
#'
#' @note If there is evidence in the prediction dataset that was not in the
#'   training dataset, the likelihood weighted sampling procedure used for
#'   prediction fails. This particularly happens if large parts of the training
#'   dataset are missing or if nodes have many discrete parents (or both). In
#'   the former case, using parametric EM to fit the BN can resolve the issue.
#'   In the latter case one can set unidentifiable parameters using a uniform
#'   prior
#'
#' @examples
computePredictions <- function(varName, dat, dag.fit, train = FALSE, 
                               particles = 100, parallel = FALSE){
  colNum <- which(names(dat) == varName)
  cat("Predicting", varName, fill = TRUE)
  
  if(train){
    fits <- dag.fit[[colNum]]$fitted.values
    actuals <- dat[[colNum]]
    
  } else {
    #Only impute when an actual is available (otherwise it is a waste of time)
    complete.actuals <- complete.cases(dat[, colNum])
    dat.toImpute <- dat[complete.actuals,]
    actuals <- dat.toImpute[[colNum]]
    
    #Simultaneously impute missing input and predict output
    dat.toImpute[, colNum] <- as.numeric(NA)
    
    if(parallel){
      dat.imputed <- parImpute(dag.fit, dat.toImpute, particles = particles)
    } else {
      dat.imputed <- impute(dag.fit, dat.toImpute, method = "bayes-lw", n = particles)
    }
    fits <- dat.imputed[[colNum]]
  }
  
  return(data.frame(fits = fits, actuals = actuals))
}

#' Computes the normalized root mean square error for given predictions
#'
#' Normalizes the RMSE by dividing by the range of the observed values. Only
#' applicable to predictions of continuous variables.
#'
#' @param pred Dataframe with a column fits containing fitted values and a
#'   column actuals containing observed values. Input is assumed to be complete.
#'
#' @return Scalar normalized RMSE
#'
#' @export
#'
#' @examples
computeRMSE <- function(pred){
  fits <- pred$fits
  actuals <- pred$actuals
  
  RMSE <- sqrt(mean((fits - actuals)^2, na.rm = TRUE))
  normRMSE <- RMSE/(max(actuals, na.rm = TRUE) - min(actuals, na.rm = TRUE))
  return(normRMSE)
}

#' Computes average validation RMSE for given bootstrap DAGs and threshold
#'
#' Obtains averaged DAG, fits parameters using parametric EM and computes
#' out-of-sample predictions for a given validation set. Returns RMSE averaged
#' over nodes to predict.
#'
#' @inheritParams fit.boot
#' @param dat.train Dataframe to use for parameter fitting
#' @param dat.val Dataframe to use for prediction
#' @inheritParams em.parametric
#' @param penalty Optional argument specifying which penalty was used to learn
#'   DAGs
#' @param predictRange Indices of columns in dat.val that should be predicted
validate <- function(threshold, samples, dat.train, dat.val, particles, obs.perParam, 
                     penalty = "", predictRange = 8:24){
  case <- paste(penalty, "_", threshold, sep = "")
  cat("Penalty_threshold:", case, "\n")
  
  #Average bootstrapped DAGS for a given threshold and fit parameters using a single
  #iteration of parametric EM
  dag.fit <- fit.boot(samples, threshold, dat.train, particles = particles,
                      obs.perParam = obs.perParam, parallel = TRUE)
  
  #Compute predictions for each variable using fitted DAG
  pred.list <- lapply(names(dat.val)[predictRange], computePredictions, 
                      dat = dat.val, dag.fit = dag.fit,
                      particles = particles, parallel = TRUE)
  
  #Compute average RMSE for fitted DAG
  RMSE.av <- mean(sapply(pred.list, computeRMSE))
  print(RMSE.av)
  return(RMSE.av)
}