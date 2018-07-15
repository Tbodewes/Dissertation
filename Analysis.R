#TODO: EM for different penalizations

invisible(lapply(c("Scoring.R", "Fit DAG.R", "Experiments.R"), source))
library(dplyr)
library(ggplot2)


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
dag.names <- list("Alarm", "Ecoli", "Sangiovese")
k.vec <- c(10, 50, 100, 250)
beta.vec <- c(0, 0.2, 0.4)
replications <- 5
penalties <- c("0.25", "0.45", "0.75", "bic")

penalties.em <-  c("0.35")
k.vec.em <- c(10, 50, 100)
beta.vec.em <- c(0.1, 0.2, 0.4)

dag.vec.pruned <- lapply(dag.vec, pruneFit, max.parents = 3)
ecoli.pruned <- pruneFit(ecoli.fit, 2)
alarm.pruned <- pruneFit(alarm.fit, 2)

#Run unordered fitting experiment
#7,900 seconds
system.time(result.unordered.alarm <- experiment.full(dag.vec = dag.vec[1], 
                                          dag.names = dag.names[1],
                                          dat.vec = dat.vec[1],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = FALSE))

saveRDS(result.unordered.alarm, file = "Result_unordered_alarm.RDS")

#19,000 seconds
system.time(result.unordered.ecoli <- experiment.full(dag.vec = dag.vec[2], 
                                          dag.names = dag.names[2],
                                          dat.vec = dat.vec[2],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = FALSE))

saveRDS(result.unordered.ecoli, file = "Result_unordered_ecoli.RDS")

# 3100 seconds
system.time(result.unordered.sangiovese <- experiment.full(dag.vec = dag.vec[3], 
                                          dag.names = dag.names[3],
                                          dat.vec = dat.vec[3],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = FALSE))

saveRDS(result.unordered.sangiovese, file = "Result_unordered_sangiovese.RDS")

result.unordered.alarm <- readRDS("Result_unordered_alarm.RDS")
result.unordered.ecoli <- readRDS("Result_unordered_ecoli.RDS")
result.unordered.sangiovese <- readRDS("Result_unordered_sangiovese.RDS")

result.unordered <- rbind(result.unordered.alarm, result.unordered.ecoli,
                               result.unordered.sangiovese)
saveRDS(result.unordered, file = "Result_unordered.RDS")

head(result.unordered)

pdf(file = "graph_unordered.pdf", height = 9, width = 9)
result.unordered %>% group_by(dag, k, beta, penalty) %>% 
  summarize(SHD.av = mean(SHD.NAL), T.av = mean(T.NAL)) %>%
  ggplot(aes(k, SHD.av)) + facet_grid(dag~beta, scales = "free_y") +
  geom_line(aes(colour = penalty))
dev.off()

#Run ordered fitting experiment
#6,500 seconds
system.time(result.ordered.alarm <- experiment.full(dag.vec = list(alarm.pruned), 
                                          dag.names = dag.names[1],
                                          dat.vec = dat.vec[1],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = TRUE))

saveRDS(result.ordered.alarm, file = "Result_ordered_alarm.RDS")

#Run ecoli experiment with max parental set size of 2 for computational reasons
#25,000 seconds
system.time(result.ordered.ecoli <- experiment.full(dag.vec = list(ecoli.pruned), 
                                          dag.names = dag.names[2],
                                          dat.vec = dat.vec[2],
                                          k.vec = k.vec,  beta.vec = beta.vec,
                                          replications = replications,
                                          penalties = penalties, str.em = FALSE,
                                          parallel = TRUE, ordered = TRUE))

saveRDS(result.ordered.ecoli, file = "Result_ordered_ecoli.RDS")

#3,300 seconds
system.time(result.ordered.sangiovese <- experiment.full(dag.vec = dag.vec.pruned[3], 
                                               dag.names = dag.names[3],
                                               dat.vec = dat.vec[3],
                                               k.vec = k.vec,  beta.vec = beta.vec,
                                               replications = replications,
                                               penalties = penalties, str.em = FALSE,
                                               parallel = TRUE, ordered = TRUE))

saveRDS(result.ordered.sangiovese, file = "Result_ordered_sangiovese.RDS")

result.ordered.alarm <- readRDS("Result_ordered_alarm.RDS")
result.ordered.ecoli <- readRDS("Result_ordered_ecoli.RDS")
result.ordered.sangiovese <- readRDS("Result_ordered_sangiovese.RDS")

result.ordered <- rbind(result.ordered.alarm, result.ordered.ecoli,
                               result.ordered.sangiovese)
saveRDS(result.ordered, file = "Result_ordered.RDS")

pdf("graph_ordered.pdf", height = 9, width = 9)
result.ordered %>% group_by(dag, k, beta, penalty) %>% 
  summarize(SHD.av = mean(SHD.NAL), T.av = mean(T.NAL)) %>%
  ggplot(aes(k, SHD.av)) + facet_grid(dag~beta) +
  geom_line(aes(colour = penalty))
dev.off()

#Run EM experiment
#32,000 seconds
system.time(result.em.alarm <- experiment.full(dag.vec = dag.vec[1], dag.names = dag.names[1],
                                   dat.vec = dat.vec[1],
                                   k.vec = k.vec.em,  beta.vec = beta.vec.em,
                                   replications = replications, 
                                   penalties = penalties.em, str.em = TRUE,
                                   parallel = FALSE, ordered = FALSE))
saveRDS(result.em.alarm, file = "Result_EM_alarm.RDS")

system.time(result.em.ecoli <- experiment.full(dag.vec = dag.vec[2], dag.names = dag.names[2],
                                   dat.vec = dat.vec[2],
                                   k.vec = k.vec.em,  beta.vec = beta.vec.em,
                                   replications = replications, 
                                   penalties = penalties.em, str.em = TRUE,
                                   parallel = TRUE, ordered = FALSE))
saveRDS(result.em.ecoli, file = "Result_EM_ecoli.RDS")

#5,000 seconds
system.time(result.em.sangiovese <- experiment.full(dag.vec = dag.vec[3], 
                                                    dag.names = dag.names[3],
                                    dat.vec = dat.vec[3],
                                    k.vec = k.vec.em,  beta.vec = beta.vec.em,
                                    replications = replications, 
                                    penalties = penalties.em, str.em = TRUE,
                                    parallel = TRUE, ordered = FALSE))
saveRDS(result.em.sangiovese, file = "Result_EM_sangiovese.RDS")

result.em.alarm <- readRDS("Result_EM_alarm.RDS")
result.em.ecoli <- readRDS("Result_EM_ecoli.RDS")
result.em.sangiovese <- readRDS("Result_EM_sangiovese.RDS")

result.em <- rbind(result.em.alarm, result.em.ecoli,
                               result.em.sangiovese)
saveRDS(result.em, file = "Result_em.RDS")

p.em <- result.em %>% 
  mutate(SHD.rel = SHD.NAL/SHD.EM, Q.rel = Q.NAL/Q.EM, T.rel = T.NAL/T.EM) %>%
  group_by(dag, k, beta, penalty) %>% 
  summarize(SHD.av = mean(SHD.rel, na.rm = TRUE), Q.av = mean(Q.rel, na.rm = TRUE),
            T.av = mean(T.rel, na.rm = TRUE)) %>%
  ggplot(aes(k)) + facet_grid(dag~beta, scales = "free_y")

pdf("Graph_EM_SHD.pdf", height = 9, width = 9)  
p.em + geom_line(aes(y = SHD.av))
dev.off()

pdf("Graph_EM_Q.pdf", height = 9, width = 9)  
p.em + geom_line(aes(y = Q.av))
dev.off()

pdf("Graph_EM_T.pdf", height = 9, width = 9)  
p.em + geom_line(aes(y = T.av))
dev.off()