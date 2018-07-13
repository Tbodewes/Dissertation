invisible(lapply(c("Scoring.R", "Fit DAG.R", "Experiments.R"), source))

alarm.string <- paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                      "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                      "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                      "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                      "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                      "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                      "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
dag.alarm <- model2network(alarm.string)
alarm.fit <- bn.fit(dag.alarm, alarm)

load("ecoli70.rda")
ecoli.fit <- bn 

load("tuscania-bn.Rdata")
grapes.fit <- G.fit

dag.vec <- list(alarm.fit, ecoli.fit, grapes.fit)
dat.vec <- list(NULL, NULL, NULL)
dag.names <- list("Alarm", "Ecoli", "Grapes")
k.vec <- c(10, 50, 100, 500)
beta.vec <- c(0, 0.2, 0.4)
replications <- 20
penalties <- c("0.25", "0.45", "0.55", "0.75", "bic", "aic")

#Run unordered fitting experiment
result.df <- experiment.full(dag.vec = dag.vec, dag.names = dag.names,
                             dat.vec = dat.vec,
                             k.vec = k.vec,  beta.vec = beta.vec,
                             replications = replications,
                             penalties = penalties, str.em = FALSE,
                             parallel = TRUE, ordered = FALSE)

save(result.df, file = "Experiment_unordered")

#Run ordered fitting experiment
dag.vec.pruned <- lapply(dag.vec, pruneFit, max.parents = 3)
result.df.ordered <- experiment.full(dag.vec = dag.vec.pruned, 
                                     dag.names = dag.names, dat.vec = dat.vec,
                                     k.vec = k.vec,  beta.vec = beta.vec,
                                     replications = replications,
                                     penalties = penalties, str.em = FALSE,
                                     parallel = TRUE, ordered = TRUE)

save(result.df.ordered, file = "Experiment_ordered")

#Run EM experiment
penalties.em = c("0.25")
k.vec.em <- c(10, 50, 100)
result.em <- experiment.full(dag.vec = dag.vec, dag.names = dag.names,
                             k.vec = k.vec.em,  beta.vec = beta.vec,
                             replications = replications, dat.vec = dat.vec,
                             penalties = penalties.em, str.em = TRUE,
                             parallel = TRUE, ordered = FALSE)

save(result.em, file = "Experiment_EM")