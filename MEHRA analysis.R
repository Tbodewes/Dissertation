library(data.table)
library(tidyverse)

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

#' ** Exploratory data analysis **

head(mehra.daily)
mehra.eda <- mehra.daily %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))



library(raster)
UK.regions <- getData("GADM", country = "GBR", level = 2)
detach(package:raster)

levels(mehra.daily$Zone) %in% UK.regions$NAME_2
#' Names do not match

UK.regions.f <- fortify(UK.regions, region = "NAME_2")

geoMortality <- mehra.eda %>% select(Longitude, Latitude, CVD60) %>% 
  group_by(Longitude, Latitude) %>% 
  summarize_all(mean)

png("CVD60Map.png", height = 6, width = 5.8, units = "in", 
    res = 500)
ggplot(data = UK.regions.f, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "white", color = "grey") +
  coord_map() +
  lims(x = c(-6, 2), y = c(50, 55.5)) +
  geom_point(data = geoMortality,
             aes(x = Longitude, y = Latitude, color = CVD60), 
             size = 1.5, inherit.aes = FALSE) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  labs(x = "Longitude", y = "Longitude",
       title = "Mortality tends to be higher in the north and west of England") + 
  theme(plot.title = element_text(size=12, face = "italic")) +
  ggtitle("English air quality measurement stations and mortality rates",
          subtitle = "Mortality tends to be higher in the north and west of England") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()

png("Data_availability.png", height = 6, width = 8, units = "in", 
    res = 800)
mehra.eda %>% select(no2, o3, so2, pm10, co, pm2.5, Date, Region) %>%
  gather(key = "Type", value = "Value", -Date, -Region) %>% 
  filter(is.finite(Value)) %>%
  group_by(Date, Region, Type) %>%
  summarize(Count = n()) %>% 
  ggplot(aes(Date, Region)) + facet_wrap(~Type) +
  geom_point(aes(color = Count)) + 
  scale_color_gradient(name = "Measurements\nper day", low = "dodgerblue", high = "black") +
  ggtitle("Number of air quality measurements per region per day",
          subtitle = "Availability of air quality measurements depends on region and year") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()



lapply(mehra.daily[18:23], function(x) mean(is.finite(x)))
lapply(mehra.daily[18:23], function(x) sum(is.finite(x)))


pollutants <- names(training)[18:23]
weather <- names(training)[11:16]

png("CVD60byYear.png", height = 4, width = 8, units = "in", 
    res = 500)
mehra.eda %>% select(Date, CVD60) %>%
  sample_frac(1) %>% 
  ggplot(aes(Date, CVD60)) +
  geom_point(alpha = 1/10, shape = ".") +
  ylim(c(0,0.4)) +
  ggtitle("Distribution of mortality over time",
          subtitle = "Mortality displays a decreasing trend and seasonality") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()

pdf("CVD60byType.pdf", height = 4, width = 4)
mehra.eda %>% select(Type, CVD60) %>% 
  ggplot(aes(Type, CVD60)) +
  geom_violin() +
  coord_flip() +
  ggtitle("Dist. of mortality by type",
          subtitle = "Dist. of mortality is similar across types") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic")) 
dev.off()

pdf("CVD60byMonth.pdf", height = 4, width = 4)
mehra.eda %>% select(Season, Month, CVD60) %>% 
  ggplot(aes(Month, CVD60)) +
  geom_violin() +
  coord_flip() +
  ggtitle("Distribution of mortality by month",
          subtitle = "Mortality tends to be higher in winter") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))  
dev.off()


png("pollutantsOverTime.png", height = 4, width = 8, units = "in", 
    res = 500)
temp <- mehra.eda %>% select(one_of(pollutants), Date) %>%
  gather(key = "Pollutant", value = "Value", -Date) %>%
  na.omit() %>% 
  filter((Pollutant == "co" & Value < 7.5) | (Pollutant == "no2" & Value < 200) |
           (Pollutant == "pm10" & Value < 150) | (Pollutant == "o3" & Value < 200) | 
           (Pollutant == "pm2.5" & Value < 90) | (Pollutant == "so2" & Value < 150) )

ggplot(temp, aes(Date, Value)) +
  facet_wrap(~Pollutant, scales = "free_y") +
  geom_point(alpha = 1/10, shape = ".") +
  ggtitle("Pollutant concentrations over time",
          subtitle = "CO and SO2 show decreasing trends and all pollutants show some seasonality") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic")) 
dev.off()


png("pollutantsByType.png", height = 5.3, width = 8, units = "in", 
    res = 500)
mehra.eda %>% select(one_of(pollutants), Type) %>%
  gather(key = "Pollutant", value = "Value", -Type) %>% 
  group_by(Pollutant, Type) %>% 
  ggplot(aes(Type, Value)) +
  facet_wrap(~Pollutant, scales = "free_x") +
  geom_violin() +
  coord_flip() +
  labs(y = expression(paste("Concentration (", mu, "g /", m^3, ")", sep = ""))) +
  ggtitle("Pollutant concentrations by environment type",
          subtitle = "Distributions of pollutant concentrations vary substantially across environment types") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic")) 
dev.off()

png("weatherByMonth.png", height = 4, width = 8, units = "in",
    res = 500)
mehra.eda %>% select(one_of(weather), Month) %>% 
  gather(key = "Weather", value = "Value", -Month) %>%
  filter(Weather %in% c("blh", "ssr", "t2m")) %>% 
  ggplot(aes(Month, Value)) +
  facet_wrap(~Weather, scales = "free_x") +
  geom_violin() +
  coord_flip() +
  labs(x = "Month", y = "Meter, Watt per square meter per second and degrees Kelvin respectively") +
  ggtitle("Distribution of weather variables by month",
          subtitle = "SSR and T2M tend to be substantially higher around the summer") +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic")) 
dev.off()

png("Pairs_airquality.png", height = 6, width = 8, units = "in", 
    res = 500)
mehra.eda %>% sample_n(5e5) %>% 
  select(CVD60, one_of(pollutants)) %>%
  mutate(CVD60 = ifelse(CVD60 > 0.4, NA, CVD60),
         co = ifelse(co > 8, NA, co),
         no2 = ifelse(no2 > 200, NA, no2),
         pm10 = ifelse(pm10 > 150, NA, pm10),
         o3 = ifelse(o3 > 150, NA, o3),
         pm2.5 = ifelse(pm2.5 > 90, NA, pm2.5),
         so2 = ifelse(so2 > 150, NA, so2)) %>%
  ggpairs(lower = list(continuous = function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) +
      geom_point(..., alpha = 1/10, shape = ".")
  })) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marginal distributions and pairwise dependence of pollutant concentrations",
          subtitle = paste("Pollutant concentrations are substantially correlated with each other and",
                           "only weakly correlated with mortality")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic")) 
dev.off()

png("Pairs_weather.png", height = 6, width = 8, units = "in", 
    res = 500)
mehra.eda %>% sample_n(1e5) %>% 
  select(CVD60, one_of(weather)) %>%
  ggpairs(lower = list(continuous = function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) +
      geom_point(..., alpha = 1/20, shape = ".")
  })) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marginal distributions and pairwise dependence of weather variables",
          subtitle = paste("Weather variables and mortality seem to correlate",
                           "in 2 groups: (CVD60 + SSR + T2M) and (BLH + WS + TP + SSR)")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()


pointsWithCor <- function(data, mapping, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor <- cor(x, y, use = "complete.obs")
  
  ggplot(data = data, mapping = mapping) +
    geom_point(..., alpha = 1/20, shape = ".") +
    geom_text(data = data.frame(
      x = max(x, na.rm = TRUE),
      y = max(y, na.rm = TRUE),
      lab = round(cor, digits = 3)),
      mapping = aes(x = x, y = y, label = lab),
      hjust = 1, vjust = 1,
      size = 3,
      inherit.aes = FALSE)
}


png("Pairs_pol_weather.png", height = 6, width = 8, units = "in", 
    res = 500)
mehra.eda %>% sample_n(5e5) %>% 
  mutate(co = ifelse(co > 8, NA, co),
         no2 = ifelse(no2 > 200, NA, no2),
         pm10 = ifelse(pm10 > 150, NA, pm10),
         o3 = ifelse(o3 > 150, NA, o3),
         pm2.5 = ifelse(pm2.5 > 90, NA, pm2.5),
         so2 = ifelse(so2 > 150, NA, so2)) %>%
  ggduo(columnsX = 18:23, columnsY = 11:16, 
        types =list(continuous = pointsWithCor)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Pairwise dependencies between pollutants and weather variables",
          subtitle = paste("Between-group correlation is clearest between",
                           "pollutants O3, CO and PM2.5/10 and weather variables SSR and BLH")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()







temp <- mehra.eda %>% select(Month, ssr, t2m) %>%
  group_by(Month) %>% 
  sample_n(1E4)

ggplot(temp, aes(ssr, t2m)) +
  facet_wrap(~Month) +
  geom_point(alpha = 1/20) +
  stat_cor()

  
mehra.eda %>% select(Month, o3) %>% 
  ggplot(aes(Month, o3)) +
  geom_violin()


temp <- mehra.eda %>% select(one_of(weather), Date) %>%
  gather(key = "Weather", value = "Value", -Date) %>% 
  group_by(Date, Weather) %>% 
  sample_frac(0.05)
  # filter((Pollutant == "co" & Value < 7.5) | (Pollutant == "no2" & Value < 200) |
  #          (Pollutant == "pm10" & Value < 150) | (Pollutant == "o3" & Value < 200) | 
  #          (Pollutant == "pm2.5" & Value < 90) | (Pollutant == "so2" & Value < 150) )

ggplot(temp, aes(Date, Value)) +
  facet_wrap(~Weather, scales = "free_y") +
  geom_point(alpha = 1/10, shape = ".") +
  ggtitle("Weather variables over time")
  

#' ** Structure learning and hyperparameter tuning **
mehra.daily <- readRDS("mehra_daily.RDS")

n <- nrow(mehra.daily)
test.indices <- sample(1:n, round(0.2*n))
val.indices <- sample((1:n)[-test.indices], round(0.2*n))
train.indices <- (1:n)[-c(test.indices, val.indices)]

save(train.indices, val.indices, test.indices, file = "MEHRA_split.rdata")

load("MEHRA_split.rdata")
if(!exists("mehra.daily")) {mehra.daily <- readRDS("mehra_daily.RDS")}

training <- mehra.daily[train.indices,]
validation <- mehra.daily[val.indices,]
testing <- mehra.daily[test.indices,]

round(cbind(sapply(training, function(x) mean(!is.na(x))),
            sapply(validation, function(x) mean(!is.na(x))),
            sapply(testing, function(x) mean(!is.na(x))))[18:23, ], 2)

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

RMSE10 <- sapply(thresholds, validate, samples = bootSamples10, penalty = "0.10", 
                 dat.train = training, dat.val = validation, particles = particles,
                 obs.perParam = obs.perParam)
saveRDS(RMSE10, "RMSE10.rds")


bootSamples25 <- bootstrap.dag(samples = no.boot, dat = training, penalty = "0.25",  
                               blacklist = bl, parallel = TRUE)
saveRDS(bootSamples25, "bootSamples25.rds")

RMSE25 <- sapply(thresholds, validate, samples = bootSamples25, penalty = "0.25", 
                 dat.train = training, dat.val = validation, particles = particles,
                 obs.perParam = obs.perParam)
saveRDS(RMSE25, "RMSE25.rds")



bootSamples40 <- bootstrap.dag(samples = no.boot, dat = training, penalty = "0.40",  
                               blacklist = bl, parallel = TRUE)
saveRDS(bootSamples40, "bootSamples40.rds")

RMSE40 <- sapply(thresholds, validate, samples = bootSamples40, penalty = "0.40", 
                 dat.train = training, dat.val = validation, particles = particles,
                 obs.perParam = obs.perParam)
saveRDS(RMSE40, "RMSE40.rds")


RMSE10 <- readRDS("RMSE10.rds")
RMSE25 <- readRDS("RMSE25.rds")
RMSE40 <- readRDS("RMSE40.rds")

RMSE <- c(RMSE10, RMSE25, RMSE40)

input <- expand.grid(thresholds, penalties)

mehra.tuning <- cbind(input, RMSE)
names(mehra.tuning) <- c("Threshold", "Penalty", "RMSE")

pdf("RMSE_plot.pdf", height = 2.2, width = 7.5)
ggplot(mehra.tuning, aes(Threshold, RMSE)) + 
  geom_line(aes(color = Penalty, linetype = Penalty)) +
  geom_point(aes(color = Penalty)) +
  labs(y = "Average normalized RMSE", x = "Arc inclusion threshold") +
  ggtitle("Average normalized RMSE for different hyperparameter configurations",
          subtitle = paste("Penalties of 0.10 and 0.25 outperform 0.40 and",
                           "give models with comparable performance")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()


bootSamples10 <- readRDS("bootSamples10.rds")
bootSamples25 <- readRDS("bootSamples25.rds")
bootSamples40 <- readRDS("bootSamples40.rds")

strength10 <- custom.strength(bootSamples10, names(training)) %>% 
  anti_join(bl, by = c("from", "to"))

strength25 <- custom.strength(bootSamples25, names(training)) %>% 
  anti_join(bl, by = c("from", "to"))

ECDF.df <- rbind(strength10,strength25) %>% 
  mutate(penalty = rep(c("0.10", "0.25"), each = n()/2)) %>% 
  select(penalty, strength) %>%
  group_by(penalty) %>% 
  arrange(strength) %>% 
  mutate(index = seq(0, 1, length.out = n()))

findYforThreshold <- function(df, threshold, bl){
  df <- df %>% 
    arrange(strength) %>% 
    mutate(index = seq(0, 1, length.out = n()))
  return(with(df, min(index[strength >= threshold])))
}

thresholds <- c(attributes(strength10)$threshold, 
                attributes(strength25)$threshold,
                0.80,
                0.99)

y.thresholds <- mapply(findYforThreshold, 
                       df = rep(list(strength10, strength25), 2),
                       threshold = thresholds)

threshold.df <- data.frame(penalty = rep(c("0.10", "0.25"), 2),
                           type = rep(c("auto", "custom"), each = 2),
                           threshold = thresholds,
                           y = y.thresholds)

pdf("ECDF.pdf", height = 3.5, width = 7.5)
ggplot(threshold.df) +
  facet_grid(~penalty, labeller = label_both) +
  geom_point(data = ECDF.df, aes(strength, index)) + 
  geom_segment(aes(x = threshold, y = y,
                   xend = 0, yend = y, linetype = type)) +
  geom_segment(aes(x = threshold, y = y, 
                   xend = threshold, yend = 0, linetype = type))  +
  geom_text(data = filter(threshold.df, type == "auto"), 
            aes(x = 0, y = y - 0.04, label = round(y,2)),
            hjust = 0, size = 3) +
  geom_text(data = filter(threshold.df, type == "custom"), 
            aes(x = 0, y = y + 0.04, label = round(y,2)),
            hjust = 0, size = 3) +
  geom_text(aes(x = threshold - 0.02, y = 0, label = round(threshold, 2)),
            hjust = 1, vjust = 0, size = 3) +
  scale_linetype_manual(name = "Threshold", values = c(3,2)) +
  ylim(c(0,1)) +
  labs(x = "Arc strength", 
       y = "Prop. of arcs with strength below threshold") + 
  ggtitle("Empirical CDF of arc strengths for 0.10 and 0.25 penalties",
          subtitle = paste("0.10 penalty with optimal threshold 0.80 gives",
                           "a simpler network than 0.25 with threshold 0.99")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()




#Testing phase
penalty.best <- "0.10"
threshold.best <- 0.80

training.full <- rbind(training, validation)

bootSamples.best <- bootstrap.dag(samples = no.boot, dat = training.full,   
                                  penalty = penalty.best, blacklist = bl,
                                  parallel = TRUE)
saveRDS(bootSamples.best, paste("Bootsamples_full", penalty.best, ".rds", sep = ""))

mehra.fit <- fit.boot(bootSamples.best, threshold.best, training.full, 
                      particles = particles, obs.perParam = NULL,
                      parallel = TRUE)
saveRDS(mehra.fit, paste("Mehra_fit_full_", penalty.best,"_", threshold.best,
                         ".rds", sep = ""))


#' ** Model validation **

bootSamples.best <- readRDS(paste("Bootsamples_full", penalty.best, ".rds", sep = ""))
mehra.fit <- readRDS(paste("Mehra_fit_full_", penalty.best,"_", threshold.best,
                           ".rds", sep = ""))

mehra.dag <- subgraph(bn.net(mehra.fit), names(training)[1:23])

graphviz.plot(cpdag(mehra.dag))
#Logical arc orientations are wind direction -> pm2.5 and CVD60 -> CVD.ly
mehra.dag <- set.arc(mehra.dag, "wd", "pm2.5")
mehra.dag <- set.arc(mehra.dag, "CVD60", "CVD.ly")

nodes(mehra.dag)[8:24] <- c("LAT", "LON", "ALT", "T2M", "WS", "WD", "TP", 
                            "BLH", "SSR", "CVD60", "NO2", "O3", "SO2", "PM10",
                            "CO", "PM2.5", "CVD.ly")

svg("Mehra_dag.svg", height = 6, width = 8)
graphviz.plot(mehra.dag)
dev.off()



narcs(mehra.dag)
nparams(bn.fit(mehra.dag, training[,-24]))

prednames <- names(training)[8:24]
mehra.pred <- lapply(prednames, computePredictions, 
                     dat = testing, dag.fit = mehra.fit,
                     particles = particles, parallel = TRUE)
saveRDS(mehra.pred, paste("Mehra_pred_full_", penalty.best,"_", threshold.best,
                          ".rds", sep = ""))

mehra.trainPred <- lapply(prednames, computePredictions, 
                          dat = rbind(training, validation), dag.fit = mehra.fit,
                          particles = particles, parallel = TRUE,
                          train = TRUE)

mehra.pred <- readRDS(paste("Mehra_pred_full_", penalty.best,"_", threshold.best,
                            ".rds", sep = ""))

RMSE.train.best <- sapply(mehra.trainPred, computeRMSE)
RMSE.best <- sapply(mehra.pred, computeRMSE)

RMSE.train.vit <- c(0.21, 0.25, 0.10, 0.05, 0.04, 0.14, 0.03, 0.07, 0.11, 0.04, 0.01,
0.05, 0.01, 0.00, 0.03, 0.02, NA)
RMSE.test.vit <- c(0.49,	0.40,	0.94,	0.09,	0.19,	0.34,	0.06,	0.10,	0.13,	0.27,	0.08,
                   0.12,	0.03,	0.03,	0.06,	0.03, NA)
RMSE.df <- rbind(data.frame(Analysis = "Current", Name = prednames,
                 test = RMSE.best, train = RMSE.train.best),
      data.frame(Analysis = "Vitolo et al.", Name = prednames,
                 test = RMSE.test.vit, train = RMSE.train.vit)) %>%
  mutate(Diff = ifelse(test - train > 0, test-train, 0),
         train = ifelse(train > test, test, train)) %>% 
  gather(key = "Type", value = "Value", -Analysis, -Name) %>% 
  filter(Name != "CVD.ly")

saveRDS(RMSE.df, "RMSE_df.rds")

pdf("RMSEcomparison.pdf", height = 6, width = 8)
ggplot(RMSE.df, aes(Name, Value)) +
  facet_grid(~Analysis) +
  geom_col(data = filter(RMSE.df, Type != "test"), aes(fill = Type)) +
  coord_flip() +
  geom_text(data = filter(RMSE.df, Type == "test"), 
            aes(x = Name, y = -0.05, label = round(Value, 2)), hjust = 1,
            size = 3) +
  ylim(c(-0.15, 1.1)) +
  theme(axis.title.y = element_blank()) +
  labs(y = "RMSE") +
  scale_fill_discrete(labels = c("Out-of-sample", "In-sample")) + 
  ggtitle("RMSE of in- and out-of-sample predictions, compared to results from Vitolo et al.",
          subtitle = paste("Our model has similar performance as the model from",
                           "Vitolo et al. and seems to overfit less")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic")) +
  guides(fill = guide_legend(reverse = TRUE))
dev.off()

mehra.pred.named <- rbindlist(mapply(function(df, name){df$name = name; return(df)},
                           df = mehra.pred, name = prednames, SIMPLIFY = FALSE))

mehra.trainPred.named <- mapply(function(df, name){df$name = name; return(df)},
                             df = mehra.trainPred, name = prednames, SIMPLIFY = FALSE) %>% 
  rbindlist




png("qqplots.png", height = 5, width = 7, units = "in",
    res = 500)
mehra.trainPred.named %>% 
  filter(name != "CVD.ly") %>% 
  group_by(name) %>% 
  sample_frac(0.10) %>% 
  mutate(residuals = actuals - fits) %>%
  ggplot(aes(sample = residuals)) +
  facet_wrap(~name, scales = "free") +
  geom_qq() +
  geom_qq_line() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "Normal quantiles", y = "Residual quantiles") + 
  ggtitle("Normal QQ-plots for training residuals of continuous variables",
          subtitle = paste("Residuals show strong signs of non-normality")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()

mehra.dag <- bn.net(mehra.fit)

#Compare feature correlation and residual correlation 



mehra.residuals <- mehra.trainPred.named %>% 
  mutate(residuals = actuals - fits) %>% 
  select(name, residuals) %>%
  group_by(name) %>% 
  mutate(index = 1:n()) %>% 
  spread(key = 1, value = 2) %>% 
  select(-index)

temp <- mehra.residuals %>% 
  select("CVD60", one_of(weather), one_of(pollutants)) %>% 
  cor(.,use = "complete.obs")

diag(temp) <- c(0)

resid.cor <- temp %>% as.data.frame  %>%
  mutate(name1 = rownames(.)) %>% 
  gather(key = "name2", value = "corr", -name1) %>%
  filter(as.vector(upper.tri(temp, diag = TRUE))) %>% 
  mutate(type = "feature")

feature.cor <- mehra.daily %>% 
  select("CVD60", one_of(weather), one_of(pollutants)) %>% 
  cor(.,use = "complete.obs") %>% as.data.frame  %>%
  mutate(name1 = rownames(.)) %>% 
  gather(key = "name2", value = "corr", -name1) %>%
  filter(as.vector(lower.tri(temp, diag = FALSE))) %>% 
  mutate(type = "resid")

cor.df <- rbind(resid.cor, feature.cor) 

png("corrComparison.png", height = 4, width = 7, units = "in",
    res = 600)  
ggplot(cor.df, aes(name2, name1)) +
  geom_raster(aes(fill = abs(corr))) +
  scale_x_discrete(limits = c(rownames(temp))) +
  scale_y_discrete(limits = c(rownames(temp))) +
  scale_fill_gradient(name = "Absolute\ncorrelation", 
                      low = "white", high = "#132B43") +
  geom_segment(aes(x = "CVD60", y = "CVD60", xend = "pm2.5", yend = "pm2.5"),
               size = 1, linetype = "dashed") +
  geom_text(data = data.frame( x = "CVD60", y = "pm2.5", label = "Original correlation"),
            aes(y = y, x = x, label = label),
            size = 4, hjust = 0, fontface = "bold")+
  geom_text(data = data.frame(x = "pm2.5", y = "CVD60", label = "Residual correlation"),
            aes(y = y, x = x, label = label),
            size = 4, hjust = 1, fontface = "bold") +
  labs(x = NULL, y = NULL) + 
  ggtitle("Comparison of correlation of original features with correlation of residuals",
          subtitle = paste("The model captures most of the linear dependence between variables")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()





  ggplot(aes(name1, name2)) +
  geom_raster(aes(fill = abs(corr))) +
  scale_fill_gradient(low = "lightgrey", high = "#132B43")

#' ** Use network for inference **
pollutants <- names(training)[18:23]
weather <- names(training)[11:16]

sapply(pollutants, dsep, bn = mehra.dag, y = "Type")
sapply(pollutants, dsep, bn = mehra.dag, y = "Year")
sapply(pollutants, dsep, bn = mehra.dag, y = "Month")

#' According to the network, mortality is fully independent of:
names(training)[sapply(names(training), dsep, bn = mehra.dag, y = "CVD60")]
names(training)[sapply(names(training), dsep, bn = mehra.dag, y = "Year")]
names(training)[sapply(names(training), dsep, bn = mehra.dag, y = "Month")]

#' No pollutant is d-separated from any other pollutant
input <- expand.grid(pollutants, pollutants, stringsAsFactors = FALSE)
output <- mapply(dsep, x = input[[1]], y = input[[2]], MoreArgs = list(bn = mehra.net))
cbind(input, output)

#' SSR is d-separated from WD and BLH
input <- expand.grid(weather, weather, stringsAsFactors = FALSE)
output <- mapply(dsep, x = input[[1]], y = input[[2]], MoreArgs = list(bn = mehra.net))
cbind(input, output)[output,]

dsep("ssr", "t2m", "Month", bn = mehra.dag)

#' SSR is d-separated from PM10 and PM2.5
input <- expand.grid(weather, pollutants, stringsAsFactors = FALSE)
output <- mapply(dsep, x = input[[1]], y = input[[2]], MoreArgs = list(bn = mehra.net))
cbind(input, output)[output,]

cl <- makeCluster(detectCores() - 1)

samples.type <- cpdist(mehra.fit, nodes = c("Type", "no2", "so2"),
               evidence = list(Type = levels(training$Type)),
               method = "lw", n = 100000)

pdf("pollutantsTypeSim.pdf", height = 4, width = 8)
samples.type %>% 
  gather(key = "Pollutant", value = "Concentration", -Type) %>% 
  ggplot(aes(Type, Concentration)) +
  facet_wrap(~Pollutant, scales = "free_x") +
  geom_violin() +
  coord_flip() +
  labs(y = expression(paste("Concentration (", mu, "g /", m^3, ")", sep = "")), 
       title = "Simulated distribution of pollutants conditional on environmental type")
dev.off()


baseline <- samples.type %>% 
  select(-Type) %>% 
  summarize_all(funs(weighted.mean(., weights)), - weights)

new <- samples.type %>% 
  group_by(Type) %>% 
  summarize_all(funs(weighted.mean(., weights)), - weights)
  
cl <- makeCluster(detectCores()-1)

samples.year <- cpdist(mehra.fit, nodes = c("Year", "Month", "CVD60", pollutants),
                       evidence = list(Year = c(1981, 2014)),
                       method = "lw", cl = cl, n = 1000000)

temp <- samples.year %>% 
  mutate(Year = round(Year, 0),
         Date = as.Date(paste(Year, Month, "30", sep = "-"))) %>%
  select(-Year, -Month) %>% 
  gather(key = "Key", value = "Value", -Date) %>% 
  group_by(Date, Key) %>% 
  summarize(mean = mean(Value), sd = sd(Value))

png("mortalitySim.png", height = 3.8, width = 8, units = "in",
    res = 500)
ggplot(filter(temp, Key == "CVD60"), aes(Date, mean)) +
  geom_point(data = mehra.eda %>% sample_frac(1) %>% 
               select(Date, CVD60),
             aes(x = Date, y = CVD60), alpha = 1/50, shape = ".",
             col = "darkblue") +
  geom_line(size = 0.7) +
  geom_line(aes(y = mean + 2*sd), linetype = 2, size = 0.7) +
  geom_line(aes(y = mean - 2*sd), linetype = 2, size = 0.7) +
  ylim(c(0, 0.4)) +
  labs(title = "The model captures the overall pattern in mortality well",
       y = "Mortality rate") + 
  ggtitle("Simulated mean and two-sigma bands for mortality compared with true mortality",
          subtitle = paste("The model captures the temporal pattern in mortality well")) +
  theme(plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"))
dev.off()
  

samples.year$weights <- head(attributes(samples.year)$weights)

ggplot(samples.year) + geom_line(aes(Year, CVD60))
