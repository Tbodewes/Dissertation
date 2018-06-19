library(magrittr)
library(parallel)
source("Functions.R")


# load the data.
data(alarm)

modelstring = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                    "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                    "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                    "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
alarm.dag <- model2network(modelstring)
alarm.fit <- bn.fit(alarm.dag, alarm)

generateData <- function(dag, n, no.missing){
  df.full <- rbn(x = dag, n)
  
  replaceMiss <- function(x){x[sample(1:length(x), no.missing)] <- NA; return(x)} 
  
  df.miss <- replaceMiss  %>%
  apply(X = df.full, MARGIN =  1, .) %>% 
  t(.) %>%
  as.data.frame()
  
  return(df.miss)
}




computeF <- function(dag.true, dag.fitted){
  node.names <- nodes(dag.true)
  
  computeConMatrix.node <- function(node, dag.true, dag.fitted){
    parents.true <- incoming.arcs(dag.true, node)[,1]
    parents.fitted <- incoming.arcs(dag.fitted, node)[,1]
    TP <- sum(parents.true %in% parents.fitted)
    FP <- sum(!parents.fitted %in% parents.true)
    FN <- sum(!parents.true %in% parents.fitted)
    
    return(c(TP, FP, FN))
  }
  
  conMatrix <- colSums(t(sapply(node.names, computeConMatrix.node,
                      dag.true = dag.true, dag.fitted = dag.fitted)))
  TP <- conMatrix[1]
  FP <- conMatrix[2]
  FN <- conMatrix[3]
  
  precision <- TP/(TP + FP)
  recall <- TP/(TP + FN)
  F1 <- 2/(1/precision + 1/recall)
  
  return(F1)
}

alarmExperiment.instance <- function(n, no.missing, penalty, alarm.fit){
  cat("n:", n, "Missing:", no.missing, "Penalty:", penalty, fill = TRUE)
  
  simulatedData <- generateData(alarm.fit, n, no.missing)
  
  nodes.ordered <- node.ordering(alarm.fit)
  
  dag.fitted <- fit.dag.ordered(node.names = nodes.ordered, 
                        max.parents = 3, 
                        dat = simulatedData, 
                        penalty = penalty)
  
  F1 <- computeF(alarm.fit, dag.fitted)
  
  return(F1)
}


alarmExperiment <- function(n.vec, no.missing.vec, penalty.vec, alarm.fit){

  input <- expand.grid(no.missing.vec, penalty.vec, n.vec, 
                       stringsAsFactors = FALSE)
  names(input) <- c("no.missing", "penalty", "n")
  input <- input[c("n", "penalty", "no.missing")]
  
  output <- mapply(alarmExperiment.instance, 
                   n = input$n, 
                   no.missing = input$no.missing,
                   penalty = input$penalty, 
                   alarm.fit = list(alarm.fit))
  
  result <- cbind(input, output)
  
  return(result)
}

n.vec <- c(500, 2500, 5000, 2.5E4, 5E4, 1E5, 2.5E5)
no.missing.vec <- c(0, 2, 4)
penalty.vec <- c(0.3, 0.4, 0.5, 0.75, "bic", "aic")


result.alarmExperiment <- alarmExperiment(n.vec, no.missing.vec,
                                          penalty.vec, alarm.fit)

debug(alarmExperiment.instance)
alarmExperiment.instance(500, 0, "0.25", alarm.fit)

save(result.alarmExperiment, file = "Alarm experiment.Rdata")
