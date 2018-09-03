library(tidyverse)
library(scales)
library(extrafont)

loadfonts(device = "win")
par(family = "serif")

invisible(lapply(c("Scoring.R", "Fit DAG.R", "Experiments.R"), source))

#' ** CONFIGURE EXPERIMENTS **
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
replications <- 20
penalties <- c("0.10", "0.25", "0.40", "0.60", "bic", "aic")

penalties.em <-  c("0.10", "0.25", "0.40")
k.vec.em <- c(10, 50, 100, 250, 500)
beta.vec.em <- c(0.05, 0.1, 0.2)

dag.vec.pruned <- mapply(pruneFit, dag = dag.vec, max.parents = c(2, 2, 3))

#' ** RUN EXPERIMENTS **

#' * Known node ordering *
system.time(result.ordered<- experiment.full(dag.vec = dag.vec.pruned, 
                                                    dag.names = dag.names,
                                                    dat.vec = dat.vec,
                                                    k.vec = k.vec,  beta.vec = beta.vec,
                                                    replications = replications,
                                                    penalties = penalties, str.em = FALSE,
                                                    parallel = TRUE, ordered = TRUE))

saveRDS(result.ordered, file = "Result_ordered.RDS")

#' * Unknown node ordering *
system.time(result.unordered<- experiment.full(dag.vec = dag.vec, 
                                             dag.names = dag.names,
                                             dat.vec = dat.vec,
                                             k.vec = k.vec,  beta.vec = beta.vec,
                                             replications = replications,
                                             penalties = penalties, str.em = FALSE,
                                             parallel = TRUE, ordered = FALSE))

saveRDS(result.unordered, file = "Result_unordered.RDS")

#' * Comparison of NAL and S-EM *
system.time(result.em <- experiment.full(dag.vec = dag.vec,
                                         dag.names = dag.names,
                                         dat.vec = dat.vec,
                                         k.vec = k.vec.em,  beta.vec = beta.vec.em,
                                         replications = replications,
                                         penalties = penalties.em, str.em = TRUE,
                                         parallel = FALSE, ordered = FALSE))
saveRDS(result.em, file = "Result_em.RDS")


#' ** PROCESS RESULTS **
result.ordered <- readRDS("Result_ordered.RDS")
result.unordered <- readRDS("Result_unordered.RDS")
result.em <- readRDS("Result_em.RDS")

names(result.ordered)[c(1, 4)] <- c("DAG", "Penalty")
result.ordered$DAG <- toupper(result.ordered$DAG)
result.ordered$Penalty <- toupper(result.ordered$Penalty)
result.ordered$beta <- paste("beta:~", result.ordered$beta, sep = "")

names(result.unordered)[c(1, 4)] <- c("DAG", "Penalty")
result.unordered$DAG <- toupper(result.unordered$DAG)
result.unordered$Penalty <- toupper(result.unordered$Penalty)
result.unordered$beta <- paste("beta:~", result.unordered$beta, sep = "")

names(result.em)[c(1, 4)] <- c("DAG", "Penalty")
result.em$DAG <- toupper(result.em$dag)
result.em$Penalty <- toupper(result.em$Penalty)
result.em$beta <- paste("beta:~", result.em$beta, sep = "")

result.em.cast <- result.em %>% 
  mutate(id = row_number()) %>%
  melt(id.vars = c("dag", "k", "beta", "Penalty", "id")) %>%
  mutate(variable = as.character(variable),
         metric = as.factor(substr(variable, 1, 3)), 
         method = as.factor(substr(variable, 5, 7))) %>%
  select(-variable)

result.em.tidy <- na.omit(spread(result.em.cast, key = method, value = value)) %>%
  mutate(ratio = NAL/SEM) %>% group_by(dag, k, beta, Penalty, metric) %>%
  summarise(ratio.av = mean(ratio), se = sd(ratio)/n(), lb = ratio.av - 2*se,
            ub = ratio.av + 2*se)

#' ** PRODUCE FIGURES **
titleSize <- 12.5
textSize <- 11.5
theme.png <- theme(plot.title = element_text(size = titleSize, face = "bold", family = "serif"), 
                   plot.subtitle = element_text(size = textSize, face = "italic", family = "serif"),
                   axis.title = element_text(size = textSize, family = "serif"),
                   legend.title = element_text(size = textSize, family = "serif"),
                   legend.text = element_text(size = textSize, family = "serif"),
                   strip.text = element_text(size = textSize, family = "serif"))

png("Graph_ordered.png", height = 6, width = 8, units = "in", 
    res = 1000)
result.ordered  %>% group_by(DAG, k, beta, Penalty) %>% 
  summarize(SHD.av = mean(SHD.NAL), se = sd(SHD.NAL)/n(), 
            lb = SHD.av - 2*se, ub = SHD.av + 2*se) %>%
  ggplot(aes(k, SHD.av)) +
  geom_line(aes(colour = Penalty, linetype = Penalty)) +
  geom_point(aes(colour = Penalty), size = 1) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, colour = Penalty), size = 0.1) +
  facet_grid(DAG~beta, labeller = labeller(beta = label_parsed, DAG = label_value)) +
  labs(x = "Relative sample size (k)", y = "Scaled SHD") +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4) +
  ggtitle("Structure learning performance of NAL for known node ordering",
          subtitle = paste("NAL often gives better results for incomplete data if",
                           "a penalization factor under 0.5 is used")) +
  theme.png
dev.off()

png(file = "Graph_unordered.png", height = 6, width = 8, units = "in", 
    res = 1000)
result.unordered  %>% group_by(DAG, k, beta, Penalty) %>% 
  summarize(SHD.av = mean(SHD.NAL), se = sd(SHD.NAL)/n(), 
            lb = SHD.av - 2*se, ub = SHD.av + 2*se) %>%
  ggplot(aes(k, SHD.av)) +
  geom_line(aes(colour = Penalty, linetype = Penalty)) +
  geom_point(aes(colour = Penalty), size = 1) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, colour = Penalty), size = 0.1) +
  facet_grid(DAG~beta, scales = "free_y", 
             labeller = labeller(beta = label_parsed, DAG = label_value)) +
  labs(x = "Relative sample size (k)", y = "Scaled SHD") +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4) +
  ggtitle("Structure learning performance of NAL for unknown node ordering",
          subtitle = paste("NAL tends to perform best for unknown node ordering", 
                           "with penalization rates under 0.25")) +
  theme.png
dev.off()


png(file = "Que_unordered.png", height = 6, width = 8, units = "in", 
    res = 1000)
result.unordered  %>% group_by(DAG, k, beta, Penalty) %>% 
  summarize(QUE.av = mean(QUE.NAL)) %>%
  ggplot(aes(k, QUE.av)) + 
  geom_line(aes(colour = Penalty, linetype = Penalty)) +
  geom_point(aes(colour = Penalty), size = 1) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, colour = Penalty), size = 0.1) +
  facet_grid(DAG~beta, scales = "free_y", 
             labeller = labeller(beta = label_parsed, DAG = label_value)) +
  labs(x = "Relative sample size (k)", y = "Calls to scoring function") +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4) +
  ggtitle("Computational cost of NAL for unknown node ordering",
          subtitle = "Tabu search with NAL scoring converges faster for harsher penalization") +
  theme.png
dev.off()



emplot.height <- 4.5

png("Graph_EM_SHD.png", height = emplot.height, width = 8, units = "in", 
    res = 1000) 
filter(result.em.tidy, metric == "SHD") %>% 
  ggplot(aes(k, ratio.av)) +
  geom_line(aes(colour = Penalty, linetype = Penalty)) +
  geom_point(aes(colour = Penalty), size = 1) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, colour = Penalty), size = 0.1) +
  facet_grid(dag~beta, scales = "free_y", 
             labeller = labeller(beta = label_parsed, dag = label_value)) +
  labs(x = "Relative sample size (k)", y = "SHD ratio (NAL / S-EM)")+ 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4) +
  ggtitle("Comparison of NAL and S-EM on structure learning performance",
          subtitle = paste("NAL with appropriate penalization", 
                           "learns graphs of similar quality as Structural EM")) +
  theme.png +
  scale_color_manual(values = hue_pal()(6)[1:3])
dev.off()


png("Graph_EM_Q.png", height = emplot.height, width = 8, units = "in", 
    res = 1000) 
filter(result.em.tidy, metric == "QUE") %>% 
  ggplot(aes(k, ratio.av)) +
  geom_line(aes(colour = Penalty, linetype = Penalty)) +
  geom_point(aes(colour = Penalty), size = 1) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, colour = Penalty), size = 0.1) +
  facet_grid(dag~beta, scales = "free_y", 
             labeller = labeller(beta = label_parsed, dag = label_value)) +
  labs(x = "Relative sample size (k)", y = "Calls to scoring function ratio (NAL / S-EM)")+ 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4) +
  ggtitle("Comparison of NAL and S-EM on number of calls to scoring function",
          subtitle = paste("NAL with appropriate penalization", 
                           "requires 50-75% fewer calls to the scoring function than S-EM")) +
  theme.png +
  scale_color_manual(values = hue_pal()(6)[1:3])
dev.off()

png("Graph_EM_T.png", height = emplot.height, width = 8, units = "in", 
    res = 1000) 
filter(result.em.tidy, metric == "TIM") %>% 
  ggplot(aes(k, ratio.av)) +
  geom_line(aes(colour = Penalty, linetype = Penalty)) +
  geom_point(aes(colour = Penalty), size = 1) +
  #geom_errorbar(aes(ymin = lb, ymax = ub, colour = Penalty), size = 0.1) +
  facet_grid(dag~beta, scales = "free_y", 
             labeller = labeller(beta = label_parsed, dag = label_value)) +
  labs(x = "Relative sample size (k)", y = "Computation time ratio (NAL / S-EM)")+ 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 5, alpha = 1/4) +
  ggtitle("Comparison of NAL and S-EM on computation time",
          subtitle = paste("NAL with appropriate penalization converges four to ten times", 
                           "faster than S-EM")) +
  theme.png +
  scale_color_manual(values = hue_pal()(6)[1:3])
dev.off()