#TODO: EM for different penalizations

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
sangiovese.fit <- G.fit

dag.vec <- list(alarm.fit, ecoli.fit, sangiovese.fit)
dat.vec <- list(NULL, NULL, NULL)
dag.names <- list("Alarm", "Ecoli", "sangiovese")
k.vec <- c(10, 50, 100, 500)
beta.vec <- c(0, 0.2, 0.4)
replications <- 5
penalties <- c("0.25", "0.45", "0.75", "bic", "aic")

#Run unordered fitting experiment
result.unordered.alarm <- experiment.full(dag.vec = dag.vec[1], 
                                          dag.names = dag.names[1],
                                          dat.vec = dat.vec[1],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = FALSE)

save(result.unordered.alarm, file = "Result_unordered_alarm")

result.unordered.ecoli <- experiment.full(dag.vec = dag.vec[2], 
                                          dag.names = dag.names[2],
                                          dat.vec = dat.vec[2],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = FALSE)

save(result.unordered.ecoli, file = "Result_unordered_ecoli")

result.unordered.sangiovese <- experiment.full(dag.vec = dag.vec[3], 
                                          dag.names = dag.names[3],
                                          dat.vec = dat.vec[3],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = FALSE)

save(result.unordered.sangiovese, file = "Result_unordered_sangiovese")

result.unordered <- rbind(result.unordered.alarm, result.unordered.ecoli,
                               result.unordered.sangiovese)
save(result.unordered, file = "Result_unordered")


#Run ordered fitting experiment
dag.vec.pruned <- lapply(dag.vec, pruneFit, max.parents = 3)
result.ordered.alarm <- experiment.full(dag.vec = dag.vec.pruned[1], 
                                          dag.names = dag.names[1],
                                          dat.vec = dat.vec[1],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = TRUE)

save(result.ordered.alarm, file = "Result_ordered_alarm")

result.ordered.ecoli <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                          dag.names = dag.names[2],
                                          dat.vec = dat.vec[2],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = TRUE)

save(result.ordered.ecoli, file = "Result_ordered_ecoli")

result.ordered.sangiovese <- experiment.full(dag.vec = dag.vec.pruned[3], 
                                               dag.names = dag.names[3],
                                               dat.vec = dat.vec[3],
                                               k.vec = k.vec,  beta.vec = beta.vec,
                                               replications = replications,
                                               penalties = penalties, str.em = FALSE,
                                               parallel = TRUE, ordered = TRUE)

save(result.ordered.sangiovese, file = "Result_ordered_sangiovese")

result.ordered <- rbind(result.ordered.alarm, result.ordered.ecoli,
                               result.ordered.sangiovese)
save(result.ordered, file = "Result_ordered")

#Run EM experiment
penalties.em <-  c("0.25")
k.vec.em <- c(10, 50, 100)
beta.vec.em <- c("0.2", "0.4")

result.em.alarm <- experiment.full(dag.vec = dag.vec[1], dag.names = dag.names[1],
                                   dat.vec = dat.vec[1],
                                   k.vec = k.vec.em,  beta.vec = beta.vec,
                                   replications = replications, 
                                   penalties = penalties.em, str.em = TRUE,
                                   parallel = TRUE, ordered = FALSE)
save(result.em.alarm, file = "Result_EM_alarm")

result.em.ecoli <- experiment.full(dag.vec = dag.vec[2], dag.names = dag.names[2],
                                   dat.vec = dat.vec[2],
                                   k.vec = k.vec.em,  beta.vec = beta.vec,
                                   replications = replications, 
                                   penalties = penalties.em, str.em = TRUE,
                                   parallel = TRUE, ordered = FALSE)
save(result.em.ecoli, file = "Result_EM_ecoli")

result.em.sangiovese <- experiment.full(dag.vec = dag.vec[3], dag.names = dag.names[3],
                                    dat.vec = dat.vec[3],
                                    k.vec = k.vec.em,  beta.vec = beta.vec,
                                    replications = replications, 
                                    penalties = penalties.em, str.em = TRUE,
                                    parallel = TRUE, ordered = FALSE)
save(result.em.sangiovese, file = "Result_EM_sangiovese")

result.em <- rbind(result.em.alarm, result.em.ecoli,
                               result.em.sangiovese)
save(result.em, file = "Result_em")
