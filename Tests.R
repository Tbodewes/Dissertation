source("Scoring.R")
source("Fit DAG.R")
library(Rgraphviz)
library(profvis)

#' **Test scoring function and local search for discrete data**
dat.discrete <- learning.test
dag.discrete <- model2network("[A][C][F][B|A][D|A:C][E|B:F]")

#' Test by comparing to inbuilt function in complete data case
computeNAL(dag.discrete, dat.discrete)
logLik(dag.discrete, dat.discrete)/nrow(dat.discrete)

#' Introduce missing values
n <- dim(dat.discrete)[1]
p <- dim(dat.discrete)[2]
prob.missing <- 0.2
miss <- as.data.frame(matrix(rbinom(n*p, 1, prob.missing), nrow = n, ncol = p))
discrete.miss <- dat.discrete
discrete.miss[miss == 1] <- NA

#' Verify that function works with missing values, while logLik does not
computeNAL(dag.discrete, discrete.miss)
try(logLik(dag.discrete, discrete.miss))
#' computeNAL works well and gives a good approximation to the average
#' log-likelihood found using the full dataset if prob.missing is not too large.
#' If prob.missing is close to 1, there will be nodes and parent configurations
#' such that there are no valid observations. In this case the log-likelihood is
#' NaN (a warning is given)

#' In full data case, computed scores agree with inbuilt functions for discrete
#' dataset
computeScore(dag.discrete, dat.discrete, penalty = "bic")
BIC(dag.discrete, dat.discrete)/nrow(dat.discrete)

computeScore(dag, dat.discrete, penalty = "aic")
AIC(dag.discrete, dat.discrete)/nrow(dat.discrete)

#' Local search gives correct DAG up to Markov equivalence in discrete case
discrete.fitted <- fit.dag(dat.discrete, "bic", parallel = FALSE)
graphviz.compare(cpdag(discrete.fitted), cpdag(dag.discrete))

#' Check that number of comparisons is similar to inbuilt function
discrete.fitted$learning$ntests
tabu(dat.discrete)$learning$ntests
#' Number of comparisons is 30% larger


#Load alarm data and true DAG
dat.alarm <- alarm
dag.alarm <- model2network(paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                                 "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                                 "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                                 "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                                 "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                                 "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                                 "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = ""))

#' Fit alarm DAG and check that parallel implementation works and is faster
#' Also run serial case twice to compare if this gives the same results
system.time(alarm.fitted.ser1 <- fit.dag(dat.alarm, "bic", parallel = FALSE))
system.time(alarm.fitted.ser2 <- fit.dag(dat.alarm, "bic", parallel = FALSE))
system.time(alarm.fitted.par <- fit.dag(dat.alarm, "bic", parallel = TRUE))
system.time(alarm.fitted.inbuilt <- tabu(dat.alarm))

#' Check that number of comparisons is similar to inbuilt function
alarm.fitted.ser1$learning$ntests
alarm.fitted.ser2$learning$ntests
alarm.fitted.par$learning$ntests
alarm.fitted.inbuilt$learning$ntests
graphviz.compare(alarm.fitted.ser1, alarm.fitted.par)
#' Make similar number of comparisons (3759 vs 4039), though slight difference
#' will already distort computational experiments. Number of queries stays
#' the same between runs. Parallel implementation also finds the same number
#' of queries. Finally, serial and parallel implementation learn same graph

shd(alarm.fitted.ser1, dag.alarm)
shd(alarm.fitted.inbuilt, dag.alarm)
#' Inbuilt method also gives higher SHD. Need to correct these two for valid 
#' comparison with structural EM, or write custom implementation of structural
#' EM

alarm.fitted <- fit.dag(dat.alarm, "bic", restarts = 3, parallel = FALSE)
graphviz.compare(cpdag(alarm.fitted), cpdag(dag.alarm))



shd(alarm.fitted, dag.alarm)
shd(alarm.fitted.inbuilt, dag.alarm)


#' ** Test scoring function and local search for Gaussian data **
dat.gauss <- gaussian.test
dag.gauss <- model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")

#' NAL and inbuilt function agree exactly in full-data case for Gaussian data
computeNAL(dag.gauss, dat.gauss)
logLik(dag.gauss, dat.gauss)/nrow(dat.gauss)

#' Introduce missing values
n <- dim(dat.gauss)[1]
p <- dim(dat.gauss)[2]
prob.missing <- 0.2
miss <- as.data.frame(matrix(rbinom(n*p, 1, prob.missing), nrow = n, ncol = p))
gauss.miss <- dat.gauss
gauss.miss[miss == 1] <- NA

#' Verify that function works with missing values, while logLik does not
computeNAL(dag.gauss, gauss.miss)
try(logLik(dag.gauss, gauss.miss))
#' NAL gives approximately same likelihood as before while logLik fails 

#' In full data case, computed scores agree with inbuilt functions 
computeScore(dag.gauss, dat.gauss, penalty = "bic")
BIC(dag.gauss, dat.gauss)/nrow(dat.gauss)

computeScore(dag, dat.gauss, penalty = "aic")
AIC(dag.gauss, dat.gauss)/nrow(dat.gauss)
#' Values disagree very slightly. Must be due to different way of counting
#' parameters
#' TODO: Check in bnlearn code how parameters for Gaussians are counted

#' Local search learns correct graph for Gaussian data
gauss.fitted <- fit.dag(dat.gauss, "bic", parallel = FALSE)
graphviz.compare(cpdag(gauss.fitted), cpdag(dag.gauss))

#' ** Test scoring function and local search for mixed data **
dat.mix <- clgaussian.test
dag.mix <- model2network("[A][B][C][H][D|A:H][F|B:C][E|B:D][G|A:D:E:F]")

#' Very minor disagreement, not sure why. Could be due to rounding somewhere
computeNAL(dag.mix, dat.mix)
logLik(dag.mix, dat.mix)/nrow(dat.mix)

#' Introduce missing values
n <- dim(dat.mix)[1]
p <- dim(dat.mix)[2]
prob.missing <- 0.2
miss <- as.data.frame(matrix(rbinom(n*p, 1, prob.missing), nrow = n, ncol = p))
mix.miss <- dat.mix
mix.miss[miss == 1] <- NA

#' Verify that function works with missing values, while logLik does not
computeNAL(dag.mix, mix.miss)
try(logLik(dag.mix, mix.miss))
#' NAL works and gives value close to full-data one, while logLik fails

#' In full data case, computed scores agree with inbuilt functions
computeScore(dag.mix, dat.mix, penalty = "bic")
BIC(dag.mix, dat.mix)/nrow(dat.mix)

computeScore(dag.mix, dat.mix, penalty = "aic")
AIC(dag.mix, dat.mix)/nrow(dat.mix)

#' Local search learns correct graph for mixed data
mix.fitted <- fit.dag(dat.mix, "bic", parallel = FALSE)
graphviz.compare(cpdag(mix.fitted), cpdag(dag.mix))

#' ** Tests dag fitting for known node order and restricted parental set size **

#' Ordered fitting algorithm finds precisely the right DAG in small case
nodes.ordered <- node.ordering(dag.mix)
order.fitted <- fit.dag.ordered(nodes.ordered, max.parents = 4, dat.mix, "bic")
graphviz.compare(cpdag(order.fitted), cpdag(dag.mix))


#' Determine maximum number of parents
max(sapply(nodes(dag.alarm), function(node, dag)
  {length(parents(dag, node))}, dag = dag.alarm))
#' Max parents in alarm is 4, but this is computationally infeasible. Use 3

nodes.ordered.alarm <- node.ordering(dag.alarm) 
alarm.order.fitted <- fit.dag.ordered(nodes.ordered.alarm, max.parents = 3,
                                      dat.alarm, "bic")
graphviz.compare(cpdag(alarm.order.fitted), cpdag(dag.alarm))
#' Gives almost correct graph (2 FN, 1 FP). Found configurations actually
#' score higher than the true configurations, should be resolved in larger
#' datasets.
computeScore.node("CCHL", dag.alarm, dat.alarm, "bic")
computeScore.node("CCHL", alarm.order.fitted, dat.alarm, "bic")

computeScore.node("HYP", dag.alarm, dat.alarm, "bic")
computeScore.node("HYP", alarm.order.fitted, dat.alarm, "bic")

#' ** Test white- and blacklisting **

#Graph obeys white- and blacklisting rules
blacklist <- data.frame(matrix(c("A", "B"), 1))
whitelist <- data.frame(matrix(c("A", "E"), 1))
discrete.fitted.restricted <- fit.dag(dat.discrete, "bic", parallel = FALSE,
                                      blacklist = blacklist, 
                                      whitelist = whitelist)
graphviz.compare(discrete.fitted.restricted, dag.discrete)


#' **Test internal functions hill climber and tabu search **

testGraph <- empty.graph(c("A", "B", "C"))
arcs(testGraph) <- matrix(c("A", "B", "B", "C"), byrow = T, ncol = 2)
node.from <- "A"

#Test findMoveType and makeMove
node.to <- "C"
moveType <- findMoveType(node.from, node.to, testGraph) #should be add
modelstring(makeMove(node.from, node.to, testGraph, moveType)) 
#should be [A][B|A][C|A:B]

node.to <- "B"
moveType <- findMoveType(node.from, node.to, testGraph) #should be delete
modelstring(makeMove(node.from, node.to, testGraph, moveType)) 
#should be [A][B][C|B]

node.from <- "C"
moveType <- findMoveType(node.from, node.to, testGraph) #should be reverse
modelstring(makeMove(node.from, node.to, testGraph, moveType)) 
#should be [A][C][B|A:C]


node.to <- "A"
findMoveType(node.from, node.to, testGraph) #should be add
try(modelstring(makeMove(node.from, node.to, testGraph, moveType)))
#should be error


#Test computeNoQueries
testGraph <- empty.graph(c("A","B", "C"))
moves <- initializeMoves(nodes(testGraph))
computeNoQueries(moves, testGraph) #Should be 6
testGraph <- set.arc(testGraph, "A", "B")
#Adding an edge leads to one possible move becoming a reversal. Hence we 
#would need one extra query
computeNoQueries(moves, testGraph) #Should be 7



#Test findBestAllowedMove
testGraph <- empty.graph(LETTERS[1:3])
testGraph2 <- set.arc(testGraph, "A", "B")
tabuList <- list(cpdag(testGraph2))

graphInList(testGraph, tabuList) #Should be FALSE (not in tabu list)

graphInList(testGraph2, tabuList) #Should be TRUE (graph in tabu list)

moves <- initializeMoves(LETTERS[1:3])
moves$score <- c(1, -1, 2, -2, -3, -4)

findBestAllowedMove(moves, testGraph, tabuList) #Should be C to A, as A to B
#and B to A are not allowed (give graph that is in tabu list)

#Debugging shows that graphs in tabu phase are not Markov equivalent
#Question: does longer tabu search and more restarts lead to a better graph?
#Also, what is the impact of tabu search and restarts on number of queries and
#computation time?
system.time(result.custom.none <- fit.dag(alarm, "bic", parallel = TRUE, 
                                       tabuSteps = 0, restarts = 0))
result.custom.none$queries
shd(result.custom.none$dag, alarm.actual)
computeScore(result.custom.none$dag, alarm, "bic")
#40 seconds, 3638 queries, shd = 30 and score = -11.04

system.time(result.custom.long <- fit.dag(alarm, "bic", parallel = TRUE, 
                                       tabuSteps = 25, restarts = 10))
result.custom.long$queries
shd(result.custom.long$dag, alarm.actual)
computeScore(result.custom.long$dag, alarm, "bic")
#1308 seconds, 77534 queries, shd = 29 and score = -11.03
#Running long tabu search and large number of restarts does not seem worth it

#Profile and identify potential bottlenecks
profvis(fit.dag(alarm, "bic", parallel = FALSE))
#Majority of time (20 of 35 seconds) is spent computing counts. Most of this
#time (16s) is spent on na.omit. This is already an optimized C function
#and can likely not be improved much. It explains some of the discrepancy
#with the inbuilt functions: these do not have to spend time on omitting 
#NA entries.
#followed by 7 of 35 seconds used for initially computing score for each
#move. Finding best allowed moves accounts for ~2 seconds.

