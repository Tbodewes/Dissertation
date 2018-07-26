#TODO: EM for different penalizations

library(tidyverse)

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
dag.names <- list("Alarm", "Ecoli", "Sangiovese")
k.vec <- c(10, 50, 100, 250, 500, 1000)
beta.vec <- c(0, 0.05, 0.1, 0.2)
replications <- 10
penalties <- c("0.10", "0.25", "0.40", "0.60", "bic", "aic")

penalties.em <-  c("0.25")
k.vec.em <- c(10, 50, 100, 250, 500, 1000)
beta.vec.em <- c(0.05, 0.1, 0.2)

dag.vec.pruned <- mapply(pruneFit, dag = dag.vec, max.parents = c(2, 2, 3))


# ORDERED FITTING
system.time(result.ordered.alarm <- experiment.full(dag.vec = dag.vec.pruned[1], 
                                                    dag.names = dag.names[1],
                                                    dat.vec = dat.vec[1],
                                                    k.vec = k.vec,  beta.vec = beta.vec,
                                                    replications = replications,
                                                    penalties = penalties, str.em = FALSE,
                                                    parallel = TRUE, ordered = TRUE))

saveRDS(result.ordered.alarm, file = "Result_ordered_alarm.RDS")


#25,000 seconds
# system.time(result.ordered.ecoli <- experiment.full(dag.vec = dag.vec.pruned[2], 
#                                           dag.names = dag.names[2],
#                                           dat.vec = dat.vec[2],
#                                           k.vec = k.vec,  beta.vec = beta.vec,
#                                           replications = replications,
#                                           penalties = penalties, str.em = FALSE,
#                                           parallel = TRUE, ordered = TRUE))
# 
# saveRDS(result.ordered.ecoli, file = "Result_ordered_ecoli.RDS")


system.time(result.ordered.ecoli0 <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                                     dag.names = dag.names[2],
                                                     dat.vec = dat.vec[2],
                                                     k.vec = k.vec,  beta.vec = beta.vec[1],
                                                     replications = replications,
                                                     penalties = penalties, str.em = FALSE,
                                                     parallel = TRUE, ordered = TRUE))
saveRDS(result.ordered.ecoli0, file = "Result_ordered_ecoli0.RDS")

system.time(result.ordered.ecoli0.05 <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                                        dag.names = dag.names[2],
                                                        dat.vec = dat.vec[2],
                                                        k.vec = k.vec,  beta.vec = beta.vec[2],
                                                        replications = replications,
                                                        penalties = penalties, str.em = FALSE,
                                                        parallel = TRUE, ordered = TRUE))
saveRDS(result.ordered.ecoli0.05, file = "Result_ordered_ecoli05.RDS")

# system.time(result.ordered.ecoli0.1 <- experiment.full(dag.vec = dag.vec.pruned[2], 
#                                                        dag.names = dag.names[2],
#                                                        dat.vec = dat.vec[2],
#                                                        k.vec = k.vec,  beta.vec = beta.vec[3],
#                                                        replications = replications,
#                                                        penalties = penalties, str.em = FALSE,
#                                                        parallel = TRUE, ordered = TRUE))
# saveRDS(result.ordered.ecoli0.1, file = "Result_ordered_ecoli1.RDS")

# system.time(result.ordered.ecoli0.2 <- experiment.full(dag.vec = dag.vec.pruned[2], 
#                                                        dag.names = dag.names[2],
#                                                        dat.vec = dat.vec[2],
#                                                        k.vec = k.vec,  beta.vec = beta.vec[4],
#                                                        replications = replications,
#                                                        penalties = penalties, str.em = FALSE,
#                                                        parallel = TRUE, ordered = TRUE))
# saveRDS(result.ordered.ecoli0.2, file = "Result_ordered_ecoli2.RDS")


#3,300 seconds
system.time(result.ordered.sangiovese <- experiment.full(dag.vec = dag.vec.pruned[3], 
                                               dag.names = dag.names[3],
                                               dat.vec = dat.vec[3],
                                               k.vec = k.vec,  beta.vec = beta.vec,
                                               replications = replications,
                                               penalties = penalties, str.em = FALSE,
                                               parallel = TRUE, ordered = TRUE))

saveRDS(result.ordered.sangiovese, file = "Result_ordered_sangiovese.RDS")



# UNORDERED FITTING


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
# system.time(result.unordered.ecoli <- experiment.full(dag.vec = dag.vec[2], 
#                                                       dag.names = dag.names[2],
#                                                       dat.vec = dat.vec[2],
#                                                       k.vec = k.vec,  beta.vec = beta.vec,
#                                                       replications = replications,
#                                                       penalties = penalties, str.em = FALSE,
#                                                       parallel = TRUE, ordered = FALSE))
# 
# saveRDS(result.unordered.ecoli, file = "Result_unordered_ecoli.RDS")


system.time(result.unordered.ecoli0 <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                                     dag.names = dag.names[2],
                                                     dat.vec = dat.vec[2],
                                                     k.vec = k.vec,  beta.vec = beta.vec[1],
                                                     replications = replications,
                                                     penalties = penalties, str.em = FALSE,
                                                     parallel = TRUE, ordered = FALSE))
saveRDS(result.unordered.ecoli0, file = "Result_unordered_ecoli0.RDS")

system.time(result.unordered.ecoli0.05 <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                                        dag.names = dag.names[2],
                                                        dat.vec = dat.vec[2],
                                                        k.vec = k.vec,  beta.vec = beta.vec[2],
                                                        replications = replications,
                                                        penalties = penalties, str.em = FALSE,
                                                        parallel = TRUE, ordered = FALSE))
saveRDS(result.unordered.ecoli0.05, file = "Result_unordered_ecoli05.RDS")

system.time(result.unordered.ecoli0.1 <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                                       dag.names = dag.names[2],
                                                       dat.vec = dat.vec[2],
                                                       k.vec = k.vec,  beta.vec = beta.vec[3],
                                                       replications = replications,
                                                       penalties = penalties, str.em = FALSE,
                                                       parallel = TRUE, ordered = FALSE))
saveRDS(result.unordered.ecoli0.1, file = "Result_unordered_ecoli1.RDS")

system.time(result.unordered.ecoli0.2 <- experiment.full(dag.vec = dag.vec.pruned[2], 
                                                       dag.names = dag.names[2],
                                                       dat.vec = dat.vec[2],
                                                       k.vec = k.vec,  beta.vec = beta.vec[4],
                                                       replications = replications,
                                                       penalties = penalties, str.em = FALSE,
                                                       parallel = TRUE, ordered = FALSE))
saveRDS(result.unordered.ecoli0.2, file = "Result_unordered_ecoli2.RDS")






# 3100 seconds
system.time(result.unordered.sangiovese <- experiment.full(dag.vec = dag.vec[3], 
                                                           dag.names = dag.names[3],
                                                           dat.vec = dat.vec[3],
                                                           k.vec = k.vec,  beta.vec = beta.vec,
                                                           replications = replications,
                                                           penalties = penalties, str.em = FALSE,
                                                           parallel = TRUE, ordered = FALSE))

saveRDS(result.unordered.sangiovese, file = "Result_unordered_sangiovese.RDS")





# COMPARISON WITH EM


#32,000 seconds
# system.time(result.em.alarm <- experiment.full(dag.vec = dag.vec[1], dag.names = dag.names[1],
#                                    dat.vec = dat.vec[1],
#                                    k.vec = k.vec.em,  beta.vec = beta.vec.em,
#                                    replications = replications, 
#                                    penalties = penalties.em, str.em = TRUE,
#                                    parallel = FALSE, ordered = FALSE))
# saveRDS(result.em.alarm, file = "Result_em_alarm.RDS")

system.time(result.em.alarm0.05 <- experiment.full(dag.vec = dag.vec[1], 
                                                  dag.names = dag.names[1],
                                                  dat.vec = dat.vec[1],
                                                  k.vec = k.vec,  beta.vec = beta.vec[2],
                                                  replications = replications,
                                                  penalties = penalties.em, str.em = TRUE,
                                                  parallel = TRUE, ordered = FALSE))
saveRDS(result.em.alarm0.05, file = "Result_em_alarm05_second.RDS")

system.time(result.em.alarm0.1 <- experiment.full(dag.vec = dag.vec[1], 
                                                  dag.names = dag.names[1],
                                                  dat.vec = dat.vec[1],
                                                  k.vec = k.vec,  beta.vec = beta.vec[3],
                                                  replications = replications,
                                                  penalties = penalties.em, str.em = TRUE,
                                                  parallel = TRUE, ordered = FALSE))
saveRDS(result.em.alarm0.1, file = "Result_em_alarm1_second.RDS")

system.time(result.em.alarm0.2 <- experiment.full(dag.vec = dag.vec[1], 
                                                  dag.names = dag.names[1],
                                                  dat.vec = dat.vec[1],
                                                  k.vec = k.vec,  beta.vec = beta.vec[4],
                                                  replications = replications,
                                                  penalties = penalties.em, str.em = TRUE,
                                                  parallel = TRUE, ordered = FALSE))
saveRDS(result.em.alarm0.2, file = "Result_em_alarm2_second.RDS")




# system.time(result.em.ecoli <- experiment.full(dag.vec = dag.vec[2], dag.names = dag.names[2],
#                                    dat.vec = dat.vec[2],
#                                    k.vec = k.vec.em,  beta.vec = beta.vec.em,
#                                    replications = replications, 
#                                    penalties = penalties.em, str.em = TRUE,
#                                    parallel = TRUE, ordered = FALSE))
saveRDS(result.em.ecoli, file = "Result_em_ecoli.RDS")

system.time(result.em.ecoli0.05 <- experiment.full(dag.vec = dag.vec[2], 
                                                          dag.names = dag.names[2],
                                                          dat.vec = dat.vec[2],
                                                          k.vec = k.vec,  beta.vec = beta.vec[2],
                                                          replications = replications,
                                                          penalties = penalties.em, str.em = TRUE,
                                                          parallel = TRUE, ordered = FALSE))
saveRDS(result.em.ecoli0.05, file = "Result_em_ecoli05.RDS")

system.time(result.em.ecoli0.1 <- experiment.full(dag.vec = dag.vec[2], 
                                                         dag.names = dag.names[2],
                                                         dat.vec = dat.vec[2],
                                                         k.vec = k.vec,  beta.vec = beta.vec[3],
                                                         replications = replications,
                                                         penalties = penalties.em, str.em = TRUE,
                                                         parallel = TRUE, ordered = FALSE))
saveRDS(result.em.ecoli0.1, file = "Result_em_ecoli1.RDS")

system.time(result.em.ecoli0.2 <- experiment.full(dag.vec = dag.vec[2], 
                                                         dag.names = dag.names[2],
                                                         dat.vec = dat.vec[2],
                                                         k.vec = k.vec,  beta.vec = beta.vec[4],
                                                         replications = replications,
                                                         penalties = penalties.em, str.em = TRUE,
                                                         parallel = TRUE, ordered = FALSE))
saveRDS(result.em.ecoli0.2, file = "Result_em_ecoli2_second.RDS")





#5,000 seconds
system.time(result.em.sangiovese <- experiment.full(dag.vec = dag.vec[3], 
                                                    dag.names = dag.names[3],
                                    dat.vec = dat.vec[3],
                                    k.vec = k.vec.em,  beta.vec = beta.vec.em,
                                    replications = replications, 
                                    penalties = penalties.em, str.em = TRUE,
                                    parallel = TRUE, ordered = FALSE))
saveRDS(result.em.sangiovese, file = "Result_em_sangiovese.RDS")


#' Aggregate results

result.ordered.ecoli0 <- readRDS("Result_ordered_ecoli0.RDS")
result.ordered.ecoli0.05 <- readRDS("Result_ordered_ecoli05.RDS")
result.ordered.ecoli0.1 <- readRDS("Result_ordered_ecoli1.RDS")
result.ordered.ecoli0.2 <- readRDS("Result_ordered_ecoli2.RDS")

result.ordered.ecoli <- rbind(result.ordered.ecoli0, result.ordered.ecoli0.05,
                              result.ordered.ecoli0.1, result.ordered.ecoli0.2)

saveRDS(result.ordered.ecoli, file = "Result_ordered_ecoli.RDS")

result.unordered.ecoli0 <- readRDS("Result_unordered_ecoli0.RDS")
result.unordered.ecoli0.05 <- readRDS("Result_unordered_ecoli05.RDS")
result.unordered.ecoli0.1 <- readRDS("Result_unordered_ecoli1.RDS")
result.unordered.ecoli0.2 <- readRDS("Result_unordered_ecoli2.RDS")

result.unordered.ecoli <- rbind(result.unordered.ecoli0, result.unordered.ecoli0.05,
                                result.unordered.ecoli0.1, result.unordered.ecoli0.2)

saveRDS(result.unordered.ecoli, file = "Result_unordered_ecoli_first.RDS")


result.em.ecoli0.05 <- readRDS("Result_em_ecoli05.RDS")
result.em.ecoli0.1 <- readRDS("Result_em_ecoli1.RDS")
result.em.ecoli0.2 <- readRDS("Result_em_ecoli2.RDS")

result.em.ecoli <- rbind(result.em.ecoli0.05, result.em.ecoli0.1, result.em.ecoli0.2)

saveRDS(result.em.ecoli, file = "Result_em_ecoli.RDS")


result.ordered.alarm <- readRDS("Result_ordered_alarm.RDS")
result.ordered.ecoli <- readRDS("Result_ordered_ecoli.RDS")
result.ordered.sangiovese <- readRDS("Result_ordered_sangiovese.RDS")

result.ordered <- rbind(result.ordered, result.ordered.alarm, result.ordered.ecoli,
                        result.ordered.sangiovese)
saveRDS(result.ordered, file = "Result_ordered.RDS")


result.unordered.alarm <- readRDS("Result_unordered_alarm.RDS")
result.unordered.ecoli <- readRDS("Result_unordered_ecoli.RDS")
result.unordered.sangiovese <- readRDS("Result_unordered_sangiovese.RDS")

result.unordered <- rbind(result.unordered, result.unordered.alarm, result.unordered.ecoli,
                          result.unordered.sangiovese)
saveRDS(result.unordered, file = "Result_unordered.RDS")


result.em.alarm <- readRDS("Result_em_alarm.RDS")
result.em.ecoli <- readRDS("Result_em_ecoli.RDS")
result.em.sangiovese <- readRDS("Result_em_sangiovese.RDS")

result.em <- rbind(result.em.alarm, result.em.ecoli,
                               result.em.sangiovese)
saveRDS(result.em, file = "Result_em.RDS")


result.ordered <- readRDS("Result_ordered_first.RDS")


levels(result.ordered$penalty)[1] <- "0.10"
levels(result.unordered$penalty)[1] <- "0.10"

#' Plot results


pdf("Graph_ordered.pdf", height = 6, width = 8)
result.ordered %>% filter(beta != 0.4)  %>% group_by(dag, k, beta, penalty) %>% 
  summarize(SHD.av = mean(SHD.NAL), se = sd(SHD.NAL)/n(), 
            lb = quantile(SHD.NAL, 0.25), ub = quantile(SHD.NAL, 0.75)) %>%
  ggplot(aes(k, SHD.av)) +
  geom_line(aes(colour = penalty, linetype = penalty)) +
  geom_point(aes(colour = penalty), size = 1) +
  #geom_errorbar(aes(ymin = SHD.av - 2*se, ymax = SHD.av + 2*se, 
  #                  colour = penalty, width = .1)) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, 
  #                  colour = penalty, width = .1)) +
  facet_grid(dag~beta, labeller = label_both) +
  labs(x = "Relative sample size (k)", y = "Scaled SHD") +
  ggtitle("Results for known node order") +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4)
dev.off()

pdf(file = "Graph_unordered.pdf", height = 6, width = 8)
result.unordered %>% filter(beta != 0.4) %>% group_by(dag, k, beta, penalty) %>% 
  summarize(SHD.av = mean(SHD.NAL)) %>%
  ggplot(aes(k, SHD.av)) + 
  geom_line(aes(colour = penalty, linetype = penalty)) +
  geom_point(aes(colour = penalty), size = 1) +
  facet_grid(dag~beta, scales = "free_y", labeller = label_both) +
  labs(x = "Relative sample size (k)", y = "Scaled SHD") +
  ggtitle("Results for unknown node order") +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4)
dev.off()




result.em <- readRDS("Result_em.RDS")
names(result.em) <- c("dag", "k", "beta", "penalty", "SHD.NAL", "QUE.NAL", "TIM.NAL", "SHD.SEM",
                      "QUE.SEM", "TIM.SEM")

sum(is.na(result.em$SHD.SEM))

result.em.cast <- result.em %>% 
  mutate(id = row_number()) %>%
  melt(id.vars = c("dag", "k", "beta", "penalty", "id")) %>%
  mutate(variable = as.character(variable),
         metric = as.factor(substr(variable, 1, 3)), 
         method = as.factor(substr(variable, 5, 7))) %>%
  select(-variable)

result.em.tidy <- na.omit(spread(result.em.cast, key = method, value = value)) %>%
  mutate(ratio = NAL/SEM) %>% group_by(dag, k, beta, penalty, metric) %>%
  summarise(ratio.av = mean(ratio), se = sd(ratio)/n(), lb = ratio.av - 2*se,
            ub = ratio.av + 2*se)

emplot.height <- 4.5
emplot.width <- 7.5

pdf("Graph_EM_SHD.pdf", height = emplot.height, width = emplot.width) 
filter(result.em.tidy, metric == "SHD", beta != 0.4) %>% 
  ggplot(aes(k, ratio.av)) +
  geom_line(size = 0.05) + 
  facet_grid(dag ~ beta, scale = "free_y", labeller = label_both) +
  geom_point() +
  geom_errorbar(aes(ymin = lb, ymax = ub,
                    width = 4)) +
  labs(x = "Relative sample size (k)", y = "SHD of NAL / SHD of S-EM") +
  ggtitle("SHD comparison of NAL and Structural EM") + 
  geom_hline(yintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4)
dev.off()


pdf("Graph_EM_Q.pdf", height = emplot.height, width = emplot.width) 
filter(result.em.tidy, metric == "QUE", beta != 0.4) %>% 
  ggplot(aes(k, ratio.av)) +
  geom_line() + facet_grid(dag ~ beta, scale = "free_y",labeller = label_both) +
  geom_point() +
  geom_errorbar(aes(ymin = lb, ymax = ub,
                    width = 4)) +
  labs(x = "Relative sample size (k)", y = "Queries of NAL / queries of S-EM") +
  ggtitle("Computational cost comparison of NAL and Structural EM") + 
  geom_hline(yintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4)
dev.off()

pdf("Graph_EM_T.pdf", height = emplot.height, width = emplot.width) 
filter(result.em.tidy, metric == "TIM", beta != 0.4) %>% 
  ggplot(aes(k, ratio.av)) +
  geom_line() + facet_grid(dag ~ beta, scale = "free_y", labeller = label_both) +
  
  geom_point() +
  geom_errorbar(aes(ymin = lb, ymax = ub,
                    width = 4)) +
  labs(x = "Relative sample size (k)", y = "Time of NAL / time of S-EM") + 
  geom_hline(yintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4)
dev.off()



computePenalizations <- function(n.vec, penalty.vec, no.nodes = 30){
  
  singlePenalty <- function(n, penalty, no.nodes){

    bic = F
    aic = F
    alpha = NULL
    
    if(penalty == "bic"){bic = TRUE}
    else if(penalty == "aic"){aic = TRUE}
    else if(!is.na(is.numeric(penalty))){alpha = as.numeric(penalty)}
    else {stop("Penalty must be 'bic', 'aic' or a number")}
    
    if(bic){penaltyFactor <- 0.5*log(n)/n}
    else if (aic){penaltyFactor <- 1/n}
    else {penaltyFactor <- (1/no.nodes)* n^(-alpha)}
    
    return(penaltyFactor)
  }
  
  input <- expand.grid(n.vec, penalty.vec, no.nodes, stringsAsFactors = FALSE)
  names(input) <- c("n", "penalty", "nodes")
  output <- mapply(singlePenalty, n = input$n, penalty = input$penalty, 
                   no.nodes = input$nodes)
  result <- cbind(input, output)
  
  return(result)
}

penalizations <- computePenalizations(seq(1000, 1000000, 1000), penalties, 30)

pdf("penalizationRates.pdf", height = 3, width = 8)
ggplot(penalizations, (aes(x = n, y = output))) +
  geom_line(aes(colour = penalty, linetype = penalty)) + 
  scale_y_log10() +
  labs(title = paste("Penalization rate behaviour"),
       x = "Sample size (x1000)",
       y = expression(lambda[n])) +
  scale_x_continuous(labels = function(x){x/1000})
dev.off()

head(penalizations)


# Phenotype graph

phen.names <- c("NoAHR", "baseline", "106a", "126#", "1290", "146b",
                "17", "185", "199a", "19b", "25", "30a", "328")
dag.phen <- empty.graph(phen.names)
phen.arcs <- data.frame(from = c(phen.names[1], phen.names[1], phen.names[7], phen.names[13],
                                 phen.names[1], phen.names[6], phen.names[1], phen.names[9],
                                 phen.names[10], phen.names[10], phen.names[1], phen.names[4],
                                 phen.names[11], phen.names[7]),
                        to = c(phen.names[2], phen.names[3], phen.names[3], phen.names[3],
                               phen.names[5], phen.names[5], phen.names[8], phen.names[8],
                               phen.names[8], phen.names[11], phen.names[12], phen.names[12],
                               phen.names[12], phen.names[13]), stringsAsFactors = FALSE)
                                 
arcs(dag.phen) <-  phen.arcs
graphviz.plot(dag.phen)
