source("Scoring.R")
source("Fit DAG.R")
library(Rgraphviz)
library(magrittr)
library(profvis)

#' **Test scoring function for discrete data**
#' Learn an appropriate dag from discrete test dataset
dag <- hc(learning.test)

#' Test by comparing to inbuilt function in complete data case
computeNAL(dag, learning.test)
logLik(dag, learning.test)/nrow(learning.test)

#Introduce missing values
n <- dim(learning.test)[1]
p <- dim(learning.test)[2]
prob.missing <- 0.2
miss <- as.data.frame(matrix(rbinom(n*p, 1, prob.missing), nrow = n, ncol = p))
dat.miss <- learning.test
dat.miss[miss == 1] <- NA

#Verify that function works with missing values, while logLik does not
computeNAL(dag, dat.miss)
try(logLik(dag, dat.miss))

#computeNAL works well and gives a good approximation to the average log-likelihood
#found using the full dataset if prob.missing is not too large. If prob.missing
#is close to 1, there will be nodes and parent configurations such that there are
#no valid observations. In this case the log-likelihood is NaN (a warning is given)

#In full data case, computed scores agree with inbuilt functions for discrete
#dataset
computeScore(dag, learning.test, penalty = "bic")
BIC(dag, learning.test)/nrow(learning.test)

computeScore(dag, learning.test, penalty = "aic")
AIC(dag, learning.test)/nrow(learning.test)


#### Tests Balov alarm experiment ####

#Number of parent configurations to evaluate
choose(1:36, 0) + choose(1:36, 1) + choose(1:36, 2) + choose(1:36, 3)

nodes.ordered <- node.ordering(alarm.dag)

system.time(dag.fitted <- fit.dag(nodes.ordered, max.parents = 3,
                      df = alarm, penalty = "bic"))
dag.true <- subgraph(alarm.dag, nodes.ordered)


graphviz.compare(dag.true, dag.fitted)
#Gives approximately the correct graph, likely for statistical reasons

####Tests hill climber and tabu search ####

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

dag.custom <- fit.dag(learning.test, "bic", parallel = FALSE, 
                      tabuSearch = FALSE)$dag
modelstring(dag.custom)
dag.inbuilt <- hc(learning.test)
modelstring(dag.inbuilt)

graphviz.compare(dag.custom, dag.inbuilt)
graphviz.plot(cpdag(dag.custom))
graphviz.plot(cpdag(dag.inbuilt))
#Custom method finds Markov equivalent graph to inbuilt method

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

#' Test whether sequential and parallel versions execute without error and
#' give the same DAG
fit.dag(learning.test, "bic", parallel = FALSE)
fit.dag(learning.test, "bic", parallel = TRUE)


alarm.custom <- fit.dag(alarm, "bic", parallel = TRUE, 
                      tabuSearch = FALSE)$dag
alarm.inbuilt <- hc(alarm)

modelstring = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                    "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                    "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                    "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
alarm.actual = model2network(modelstring)

graphviz.compare(cpdag(alarm.custom), cpdag(alarm.inbuilt))
graphviz.compare(cpdag(alarm.actual), cpdag(alarm.custom))
graphviz.compare(cpdag(alarm.actual), cpdag(alarm.inbuilt))
#All give very different results


#Test tabu search and compare running times
data(alarm)
system.time(result.custom.tabu <- fit.dag(alarm, "bic", parallel = TRUE))
alarm.custom.tabu <- result.custom.tabu$dag
result.custom.tabu$queries
#About 35 seconds in parallel, 4034 queries (including 10 step tabu search and 
#3 random restarts)

system.time(alarm.inbuilt.tabu <- tabu(alarm))
alarm.inbuilt.tabu$learning$ntests
#About 2 seconds and 3759 queries (same order of magnitude, no random restarts)

graphviz.compare(alarm.inbuilt.tabu, alarm.custom.tabu)
#Methods give similar graphs, but often disagree on arc direction. Similar
#number of relative false positives as relative false negatives (though this
#is somewhat obscured by arcs with different directions also being marked)

graphviz.compare(alarm.actual, alarm.inbuilt.tabu)
graphviz.compare(alarm.actual, alarm.custom.tabu)
#Performance appears comparable

shd(alarm.actual, alarm.inbuilt.tabu)
shd(alarm.actual, alarm.custom.tabu)

computeScore(alarm.inbuilt.tabu, alarm, "bic")
computeScore(alarm.custom.tabu, alarm, "bic")
computeScore(alarm.actual, alarm, "bic")
#Custom method does not give graph much worse than true graph, neither in
#SHD nor in score. Custom method outperforms inbuilt method in SHD and score.
#Score for actual DAG and graph from inbuilt search is not properly defined 
#as some parent configurations have no observations, here we extend the 
#definition of score to only use nodes for which it is defined

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

debug(fit.dag)
fit.dag(learning.test, "bic", parallel = FALSE, restarts = 3)
