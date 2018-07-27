library(data.table)
library(tidyverse)



#' ** Structure learning and hyperparameter tuning **
mehra.daily <- readRDS("mehra_daily.RDS")

training <- filter(mehra.daily, Year <= 2004)
validation <- filter(mehra.daily, between(Year, 2005, 2006))
testing <- filter(mehra.daily, Year > 2006)

invisible(lapply(c("Scoring.R", "Fit DAG.R", "Experiments.R", 
                   "MEHRA functions.R"), source))

penalties <- c("0.10", "0.25", "0.40")
thresholds <- c(0.50, 0.65, 0.80, 0.95, 0.99)
no.boot <- 200
particles <- 500
obs.perParam <- 100

bl <-  {data.frame(
  "from" = c(rep("Region", 9), rep("Zone", 9), rep("Type", 9),
             rep("Year", 9), rep("Season", 9), rep("Month", 9),
             rep("Latitude", 9), rep("Longitude", 9), rep("Altitude", 9),
             rep("Day", 9), rep("t2m", 10), rep("ws", 10),
             rep("wd", 10), rep("tp", 10), rep("blh", 10),
             rep("ssr", 10), rep("no2", 10), rep("so2", 10), 
             rep("co", 10), rep("o3", 10), rep("pm10", 10),
             rep("pm2.5", 10), rep("CVD60", 22), rep("CVD.ly", 23)),
  
  "to" = c("Zone", "Type", "Year", "Season", "Month", 
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Type", "Year", "Season", "Month",  
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Year", "Season", "Month", 
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Month",
            "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", 
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season",
           "Month", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season",
           "Month",  "Latitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", 
           "Month", "Latitude", "Longitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", 
           "Month", "Latitude", "Longitude", "Altitude",
           
           "Region", "Zone", "Type", "Year", "Season", "Month", 
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", "Month", 
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", "Month", 
           "Latitude", "Longitude", "Altitude", "Day", 
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season", "Month",
           "Latitude", "Longitude", "Altitude", "Day",
           
           "Region", "Zone", "Type", "Year", "Season",
           "Month", "Latitude", "Longitude", "Altitude",
           "t2m", "ws", "wd", "tp", "blh", "ssr", "no2", 
           "o3", "so2", "co", "pm10", "pm2.5", "Day",
           
           "Region", "Zone", "Type", "Year", "Season",
           "Month", "Latitude", "Longitude", "Altitude",
           "t2m", "ws", "wd", "tp", "blh", "ssr", "no2", 
           "o3", "so2", "co", "pm10", "pm2.5", "Day", "CVD60"),
           stringsAsFactors = FALSE)
}


# Tuning phase
bootSamples10 <- bootstrap.dag(samples = no.boot, dat = training, penalty = "0.10",  
                               blacklist = bl, parallel = TRUE)
saveRDS(bootSamples10, "bootSamples10.rds")

bootSamples25 <- bootstrap.dag(samples = no.boot, dat = training, penalty = "0.25",  
                               blacklist = bl, parallel = TRUE)
saveRDS(bootSamples25, "bootSamples25.rds")

bootSamples40 <- bootstrap.dag(samples = no.boot, dat = training, penalty = "0.40",  
                               blacklist = bl, parallel = TRUE)
saveRDS(bootSamples40, "bootSamples40.rds")

bootSamples10 <- readRDS("bootSamples10.rds")
bootSamples25 <- readRDS("bootSamples25.rds")
bootSamples40 <- readRDS("bootSamples40.rds")

RMSE10 <- lapply(thresholds, validate, samples = bootSamples10, penalty = "0.10", 
                 dat.train = training, dat.val = validation, particles = particles,
                 obs.perParam = obs.perParam)
saveRDS(RMSE10, "RMSE10.rds")

RMSE25 <- lapply(thresholds, validate, samples = bootSamples25, penalty = "0.25", 
                 dat.train = training, dat.val = validation, particles = particles,
                 obs.perParam = obs.perParam)
saveRDS(RMSE25, "RMSE25.rds")

RMSE40 <- lapply(thresholds, validate, samples = bootSamples40, penalty = "0.40", 
                 dat.train = training, dat.val = validation, particles = particles,
                 obs.perParam = obs.perParam)
saveRDS(RMSE40, "RMSE40.rds")

#Testing phase
penalty.best <- TBD
threshold.best <- TBD 

training.full <- rbind(training.validation)

bootSamples.best <- bootstrap.dag(samples = no.boot, dat = training.full,   
                                  penalty = penalty.best, blacklist = bl,
                                  parallel = TRUE)

strength.best <- custom.strength(bootSamples.best, names(mehra.daily))

mehra.result <- validate(threshold.best, bootSamples.best, penalty = penalty.best, 
                         dat.train = training.full, dat.val = testing, 
                         particles = particles, obs.perParam = obs.perParam, 
                         returnFull = TRUE)



#' ** Data extraction and transformation **
testing <- readRDS("C:/Users/tjbodewe/Desktop/Data/testing.rds")
training <- readRDS("C:/Users/tjbodewe/Desktop/Data/training.rds")
mehra <- rbindlist(list(training, testing))
rm(testing, training)
gc()

mehra.daily <- mehra %>% select(-Hour) %>%
  group_by(Region, Zone, Type, Year, Month, Day, Latitude, Longitude, Altitude, Season) %>%
  summarise_all(funs(mean), na.rm = TRUE)

mehra.daily <- mehra.daily %>% 
  group_by(Region, Zone, Type, Latitude, Longitude, Altitude) %>%
  mutate(Year = as.numeric(Year) + 1980) %>%
  mutate(CVD.ly = lead(CVD60, 365)) %>%
  select(Region, Zone, Type, Year, Month, Day, Season, everything()) %>%
  as.data.frame()

saveRDS(mehra.daily, file = "MEHRA_daily.RDS")
rm(mehra) 
gc()



mehra.daily <- readRDS("mehra_daily.RDS")

training <- filter(mehra.daily, Year <= 2006)
validation <- filter(mehra.daily, between(Year, 2007, 2010))
testing <- filter(mehra.daily, Year > 2010)

head(mehra.daily)
lapply(training, function(x) mean(is.na(x))) #Over 99% missing for pm2.5!
lapply(validation, function(x) mean(is.na(x)))
lapply(testing, function(x) mean(is.na(x)))

mehra.bootsamples <- bootstrap.dag(training, penalty = "0.25", blacklist = bl, samples = 200) 
saveRDS(mehra.bootsamples, "MEHRA_bootsamples.RDS")

mehra.bootsamples <- readRDS("MEHRA_bootsamples.RDS")

strength <- custom.strength(mehra.bootsamples, c(names(mehra.daily), "CVD.lm", "CVD.l1000"))

dag.mehra.99 <- averaged.network(strength, c(names(mehra.daily), "CVD.lm", "CVD.l1000"), 0.99)

dag.mehra.50 <- averaged.network(strength, c(names(mehra.daily), "CVD.lm", "CVD.l1000"), 0.50)

dag.mehra.95 <- averaged.network(strength, names(mehra.daily), 0.95)

#CGN structure forces this arc to be directed
dag.mehra.99 <- set.arc(dag.mehra.99, "Month", "pm2.5")
dag.mehra.95 <- set.arc(dag.mehra.95, "Month", "pm2.5")

pdf("MEHRA_dag99.pdf", width = 10, height = 10)
graphviz.plot(subgraph(dag.mehra.99, names(mehra.daily)))
dev.off()

pdf("MEHRA_dag95.pdf", width = 15, height = 15)
graphviz.plot(subgraph(dag.mehra.95, names(mehra.daily)))
dev.off()

mehra.fit <- bn.fit(dag.mehra.99, training)
mehra.em <- em.parametric(dag.mehra.99, sample_n(training, 10000), max.iter = 1, debug = TRUE)$dag
mehra.fit[[15]]

dag.test <- fit.dag(mehra.daily, "0.25")


RMSE.train <- sapply(names(mehra.daily)[8:24], computeRMSE, 
                     dat = training, dag.fit = mehra.fit, train = TRUE)


cl <- makeCluster(detectCores() - 1)
invisible(clusterEvalQ(cl, library(bnlearn)))

RMSE.testing <- parSapply(cl, names(mehra.daily)[8:24], computeRMSE, 
                       dat = sample_n(testing, 1000), dag.fit = mehra.em, train = FALSE,
                       impute = TRUE)
stopCluster(cl)
RMSE.testing

BN <- mehra.fit
col2remove <- c()
for (nCol in names(which(sapply(testing, class) == 'factor'))){
  stringCol <- paste("BN$", nCol, "$prob", sep = "")
  levelsCol <- names(which( eval(parse(text = stringCol)) != 0 ))
  if (!all(unique(testing[,nCol]) %in% levelsCol)) {
    col2remove <- c(col2remove, nCol)
  }
}
col2remove <- which(names(testing) %in% col2remove)

#Core of problem: combinations of factors and missing values that don't 
#occur in the training set but do occur in the testing set
with(training, sum(!is.na(pm10) & !is.na(pm2.5) & Region == "East Midlands"))
with(testing, mean(!is.na(pm10) & !is.na(pm2.5) &  Region == "East Midlands"))

#' ** EDA **
#'
#' Hypotheses:
#' - Can neglect day 
#' - 

head(mehra)
mehra[1:50, -c(1,2)]

mehra.sample <- sample_n(mehra, 1E5)
summary(mehra.sample)
attach(mehra.sample)
table(Day)
table(Month)
table(Hour)
hist(Latitude)
table(Year)
hist(tp)
mehra.sample <- mutate(mehra.sample, Year.num = as.numeric(Year))

#' * Temporal variation *
#' Based on this analysis, year can be safely turned into a continuous variable. We can also
#' remove day and 

#' Seems to be an about linearly decreasing trend in mortality over the years
ggplot(mehra.sample, aes(x = Year, y = CVD60)) + geom_boxplot(varwidth = TRUE) + coord_flip()


#' No clear yearly trend in any of the pollution factors. Certainly no obvious non-linear relationship
#' Found relationship between year and many of the variables in the MEHRA paper is somewhat striking
ggplot(mehra.sample, aes(x = Year, y = o3)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Year, y = pm2.5)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Year, y = pm10)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Year, y = so2)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Year, y = no2)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Year, y = co)) + geom_boxplot(varwidth = TRUE) + coord_flip()

#' Clear monthly differences in mortality rates, especially winter vs summer. Difference within 
#' seasons seems large enough to include month instead of season 
ggplot(mehra.daily, aes(x = Month, y = CVD60)) + geom_boxplot(varwidth = TRUE)

#' Some month-related effects seem to be present in air pollution (especially o3)
ggplot(mehra.sample, aes(x = Month, y = o3)) + geom_boxplot(varwidth = TRUE)
ggplot(mehra.daily, aes(x = Month, y = pm2.5)) + geom_boxplot(varwidth = TRUE)
ggplot(mehra.sample, aes(x = Month, y = pm10)) + geom_boxplot(varwidth = TRUE)
ggplot(mehra.sample, aes(x = Month, y = so2)) + geom_boxplot(varwidth = TRUE)
ggplot(mehra.sample, aes(x = Month, y = no2)) + geom_boxplot(varwidth = TRUE)
ggplot(mehra.sample, aes(x = Month, y = co)) + geom_boxplot(varwidth = TRUE)

ggplot(mehra.daily, aes(x = Month, y = t2m)) + geom_boxplot(varwidth = TRUE)

#' Day does not seem to have any effect on outcomes of interest and can safely be dropped
ggplot(mehra.sample, aes(x = Day, y = CVD60)) + geom_boxplot(varwidth = TRUE) + coord_flip()

ggplot(mehra.sample, aes(x = Day, y = o3)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = pm2.5)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = pm10)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = so2)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = no2)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = co)) + geom_boxplot(varwidth = TRUE) + coord_flip()

ggplot(mehra.sample, aes(x = Day, y = wd)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = ws)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = t2m)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = blh)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = tp)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Day, y = ssr)) + geom_boxplot(varwidth = TRUE) + coord_flip()


# * Relations between weather and air pollution *
mehra.weather <- select(mehra.sample, wd, ws, t2m, blh, tp, ssr)
mehra.air <- select(mehra.sample, o3, pm2.5, pm10, so2, no2, co)
mehra.envir <- cbind(mehra.weather, mehra.air)
ggpairs(mehra.air)
ggpairs(mehra.weather)



#' *Geographic variation*

#' There is some variation in latitude and longitude in each region and zone. Modelling as a gaussian seems
#' not completely unreasonable (though skewness is apparent)
ggplot(mehra.sample, aes(x = Region, y = Latitude)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Region, y = Longitude)) + geom_boxplot(varwidth = TRUE) + coord_flip() 
ggplot(mehra.sample, aes(x = Zone, y = Latitude)) + geom_boxplot(varwidth = TRUE) + coord_flip()
ggplot(mehra.sample, aes(x = Zone, y = Longitude)) + geom_boxplot(varwidth = TRUE) + coord_flip()


mean(is.na(pm2.5))

dat.yorkshire <- testing %>% filter(Region == "Yorkshire and The Humber", Zone == "Yorkshire & Humberside")
dat.yorkshire <- mutate(dat.yorkshire, Hour.num = as.numeric(Hour))
dat.yorkshire <- select(dat.yorkshire, -"Hour.num <- as.numeric(Hour)")
head(dat.yorkshire)
dat.yorkshire %>% filter(Year == "2006", Month == "01") %>%
ggplot(., aes(x = as.numeric(rownames(.)), y = Hour.num)) + geom_line()
head(Day)
head(as.numeric(Day))


#' ** Transform data **
head(mehra)




lapply(mehra.daily, function(x) mean(is.na(x)))

mehra.sub <- filter(mehra.daily, Region == "East Midlands")
p1 <- filter(mehra.sub, between(Latitude, 52.55, 52.56)) %>%
  ggplot(aes(x = Date))
p1 + geom_line(aes(y = t2m))
p1 + geom_line(aes(y = blh))






mehra.sample <- sample_n(mehra.daily, 5E5)
system.time(dag.mehra <- fit.dag(select(mehra.sample, - Date, - Day), penalty = "0.25", blacklist = bl))
graphviz.plot(cpdag(dag.mehra))