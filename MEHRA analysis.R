library(data.table)
library(tidyverse)

#' ** DATA PRE-PROCESSING **
#Data originally obtained from repository associated to MEHRA paper
mehra <- readRDS("C:/Users/tjbodewe/Desktop/Data/mehra.rds") 

names(mehra)[12:24] <- toupper(names(mehra))[12:24]

mehra.daily <- mehra %>% select(-Hour) %>%
  group_by(Region, Zone, Type, Year, Month, Day, Latitude, Longitude, Altitude, Season) %>%
  summarise_all(funs(mean), na.rm = TRUE)

mehra.daily <- mehra.daily %>% 
  group_by(Region, Zone, Type, Latitude, Longitude, Altitude) %>%
  mutate(Year = as.numeric(Year) + 1980) %>%
  select(Region, Zone, Type, Year, Month, Day, Season, everything()) %>%
  as.data.frame()

saveRDS(mehra.daily, file = "MEHRA_daily.RDS")
rm(mehra) 
gc()

mehra.daily <- readRDS("MEHRA_daily.RDS")

#' ** EXPLORATORY DATA ANALYSIS **

mehra.eda <- mehra.daily %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))

titleSize <- 12.5
textSize <- 11.5
theme.png <- theme(plot.title = element_text(size = titleSize, face = "bold", family = "serif"), 
                     plot.subtitle = element_text(size = textSize, face = "italic", family = "serif"),
                     axis.title = element_text(size = textSize, family = "serif"),
                     legend.title = element_text(size = textSize, family = "serif"),
                     legend.text = element_text(size = textSize, family = "serif"),
                     strip.text = element_text(size = textSize, family = "serif"))


#' * Autocorrelation *
pollutants <- names(mehra.daily)[18:23]
weather <- names(mehra.daily)[11:16]

#' Find a station that has observations for all pollutants
detectNA <- mehra.daily %>% group_by(Longitude, Latitude) %>% summarize_all(function(x){mean(!is.na(x))})
chosenLat <- 50.37167
chosenLon <- 4.142361

#' * Generate autocorrelation plot *
autocor.hourly <- mehra %>% filter(abs(Latitude - chosenLat) < 1e-5, abs(Longitude + chosenLon) < 1e-4) %>% 
  select(one_of(pollutants, weather)) %>% 
  lapply(.,function(x){acf(x, lag.max = 40, type = "correlation", plot = FALSE, na.action = na.pass)$acf}) %>%
  data.frame() %>% 
  mutate(lag = 0:40) %>% 
  gather(key = "key", value = "value", -lag) %>% 
  mutate(type = "Hourly data")

autocor.daily <- mehra.daily %>% filter(abs(Latitude - 50.37167) < 1e-5, abs(Longitude + 4.142361) < 1e-4) %>% 
  select(one_of(pollutants, weather)) %>% 
  lapply(.,function(x){acf(x, lag.max = 40, type = "correlation", plot = FALSE, na.action = na.pass)$acf}) %>%
  data.frame() %>% 
  mutate(lag = 0:40) %>%  
  gather(key = "key", value = "value", -lag) %>% 
  mutate(type = "Daily data")

autocor <- rbind(autocor.hourly, autocor.daily)
autocor$type <- factor(autocor$type, levels = c("Hourly data", "Daily data"))
autocor$key <- factor(autocor$key, levels = c(pollutants, weather))

png("Autocorrelation_EDA.png", height = 5, width = 8, units = "in", 
    res = 500)
autocor %>% ggplot(aes(x = lag, y = key)) +
  facet_wrap(~type) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(name = "Auto-\ncorrelation", low = "white", high = "black") +
  labs(x = "Lag", y = NULL)+
  ggtitle("Autocorrelation for hourly and daily data",
          subtitle = "Hourly data has high autocorrelation, averaging by day reduces this for most variables") +
  theme.png
dev.off()

#' * Generate map with overlaid mortality *
library(raster)
UK.regions <- getData("GADM", country = "GBR", level = 2)
detach(package:raster)

UK.regions.f <- fortify(UK.regions, region = "NAME_2")

geoMortality <- mehra.eda %>% select(Longitude, Latitude, CVD60) %>% 
  group_by(Longitude, Latitude) %>% 
  summarize_all(mean)  

png("CVD60Map.png", height = 7, width = 6, units = "in", 
    res = 500)
ggplot(data = UK.regions.f, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "white", color = "grey") +
  coord_map() +
  lims(x = c(-6, 2), y = c(50, 55.5)) +
  geom_point(data = geoMortality,
             aes(x = Longitude, y = Latitude, color = CVD60), 
             size = 1.5, inherit.aes = FALSE) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Air quality measurement stations with overlaid mortality rates",
          subtitle = "Mortality tends to be higher in the north and west of England") +
  theme.png +
  theme(plot.margin = grid::unit(c(0,0,-1.5,0), "cm"))
dev.off()

#' * Generate data availability plot *
png("Data_availability.png", height = 6, width = 8, units = "in", 
    res = 700)
mehra.eda %>% select(NO2, O3, SO2, PM10, CO, PM2.5, Date, Region) %>%
  gather(key = "Type", value = "Value", -Date, -Region) %>% 
  filter(is.finite(Value)) %>%
  group_by(Date, Region, Type) %>%
  summarize(Count = n()) %>% 
  ggplot(aes(Date, Region)) + facet_wrap(~Type) +
  geom_point(aes(color = Count)) + 
  scale_color_gradient(name = "Measurements\nper day", low = "dodgerblue", high = "black") +
  ggtitle("Number of air quality measurements per region per day",
          subtitle = "Availability of air quality measurements depends on Region and Year") +
  theme.png
dev.off()


# Count number of available observations
lapply(mehra.daily[18:23], function(x) sum(is.finite(x)))
lapply(mehra.daily[18:23], function(x) mean(is.finite(x)))

#' * Generate plots for distribution of mortality *
png("CVD60byYear.png", height = 4, width = 8, units = "in", 
    res = 500)
mehra.eda %>% select(Date, CVD60) %>%
  sample_frac(1) %>% 
  ggplot(aes(Date, CVD60)) +
  geom_point(alpha = 1/10, shape = ".") +
  ylim(c(0,0.4)) +
  ggtitle("Mortality over time",
          subtitle = "Mortality displays a decreasing trend and seasonality") +
  theme.png
dev.off()

png("CVD60byType.png", height = 4, width = 4, units = "in", 
    res = 1000)
mehra.eda %>% select(Type, CVD60) %>% 
  ggplot(aes(Type, CVD60)) +
  geom_violin() +
  coord_flip() +
  ggtitle("Mortality by type",
          subtitle = "Mortality is similar across types") +
  theme.png
dev.off()

png("CVD60byMonth.png", height = 4, width = 4, units = "in", 
    res = 1000)
mehra.eda %>% select(Season, Month, CVD60) %>% 
  ggplot(aes(reverse.levels(Month), CVD60)) +
  geom_violin() +
  coord_flip() +
  xlab("Month") +
  ggtitle("Mortality by month",
          subtitle = "Mortality tends to be higher in winter") +
  theme.png
dev.off()

#' * Generate plots for distribution of pollutants *
png("pollutantsOverTime.png", height = 4, width = 8, units = "in", 
    res = 500)
temp <- mehra.eda %>% select(one_of(pollutants), Date) %>%
  gather(key = "Pollutant", value = "Value", -Date) %>%
  na.omit() %>% 
  filter((Pollutant == "CO" & Value < 7.5) | (Pollutant == "NO2" & Value < 200) |
           (Pollutant == "PM10" & Value < 150) | (Pollutant == "O3" & Value < 200) | 
           (Pollutant == "PM2.5" & Value < 90) | (Pollutant == "SO2" & Value < 150) )

ggplot(temp, aes(Date, Value)) +
  facet_wrap(~Pollutant, scales = "free_y") +
  geom_point(alpha = 1/10, shape = ".") +
  ggtitle("Pollutant concentrations over time",
          subtitle = "CO and SO2 show decreasing trends and all pollutants show some seasonality") +
  ylab(expression(paste("Concentration (", mu, "g /", m^3, ")", sep = ""))) +
  theme.png
dev.off()


png("pollutantsByType.png", height = 5.3, width = 8, units = "in", 
    res = 1000)
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
  theme.png
dev.off()


#' * Generate scatterplot matrices *
require(GGally)

#' Use custom function to round correlation to two significant digits
mycor <- function(data, mapping, alignPercent = 0.6, method = "pearson", 
                  use = "complete.obs", corAlignPercent = NULL, corMethod = NULL, 
                  corUse = NULL, sgnf=2, ...) {
  if (!is.null(corAlignPercent)) {
    stop("'corAlignPercent' is deprecated.  Please use argument 'alignPercent'")
  }
  if (!is.null(corMethod)) {
    stop("'corMethod' is deprecated.  Please use argument 'method'")
  }
  if (!is.null(corUse)) {
    stop("'corUse' is deprecated.  Please use argument 'use'")
  }
  useOptions <- c("all.obs", "complete.obs", "pairwise.complete.obs", 
                  "everything", "na.or.complete")
  use <- pmatch(use, useOptions)
  if (is.na(use)) {
    warning("correlation 'use' not found.  Using default value of 'all.obs'")
    use <- useOptions[1]
  } else {
    use <- useOptions[use]
  }
  cor_fn <- function(x, y) {
    round(cor(x, y, method = method, use = use),2)
  }
  xCol <- substring(deparse(mapping$x),2)
  yCol <- substring(deparse(mapping$y),2)
  if (GGally:::is_date(data[[xCol]]) || GGally:::is_date(data[[yCol]])) {
    if (!identical(class(data), "data.frame")) {
      data <- fix_data(data)
    }
    for (col in c(xCol, yCol)) {
      if (GGally:::is_date(data[[col]])) {
        data[[col]] <- as.numeric(data[[col]])
      }
    }
  }
  if (is.numeric(GGally:::eval_data_col(data, mapping$colour))) {
    stop("ggally_cor: mapping color column must be categorical, not numeric")
  }
  colorCol <- deparse(mapping$colour)
  singleColorCol <- ifelse(is.null(colorCol), NULL, paste(colorCol, 
                                                          collapse = ""))
  if (use %in% c("complete.obs", "pairwise.complete.obs", "na.or.complete")) {
    if (length(colorCol) > 0) {
      if (singleColorCol %in% colnames(data)) {
        rows <- complete.cases(data[,c(xCol, yCol, colorCol)])
      } else {
        rows <- complete.cases(data[,c(xCol, yCol)])
      }
    } else {
      rows <- complete.cases(data[,c(xCol, yCol)])
    }
    if (any(!rows)) {
      total <- sum(!rows)
      if (total > 1) {
        warning("Removed ", total, " rows containing missing values")
      } else if (total == 1) {
        warning("Removing 1 row that contained a missing value")
      }
    }
    data <- data[rows, ]
  }
  xVal <- data[[xCol]]
  yVal <- data[[yCol]]
  if (length(names(mapping)) > 0) {
    for (i in length(names(mapping)):1) {
      tmp_map_val <- deparse(mapping[names(mapping)[i]][[1]])
      if (tmp_map_val[length(tmp_map_val)] %in% colnames(data)) 
        mapping[[names(mapping)[i]]] <- NULL
      if (length(names(mapping)) < 1) {
        mapping <- NULL
        break
      }
    }
  }
  if (length(colorCol) < 1) {
    colorCol <- "ggally_NO_EXIST"
  }
  if ((singleColorCol != "ggally_NO_EXIST") && (singleColorCol %in% 
                                                colnames(data))) {
    cord <- plyr::ddply(data, c(colorCol), function(x) {
      cor_fn(x[[xCol]], x[[yCol]])
    })
    colnames(cord)[2] <- "ggally_cor"
    cord$ggally_cor <- signif(as.numeric(cord$ggally_cor), 
                              sgnf)
    lev <- levels(data[[colorCol]])
    ord <- rep(-1, nrow(cord))
    for (i in 1:nrow(cord)) {
      for (j in seq_along(lev)) {
        if (identical(as.character(cord[i, colorCol]), 
                      as.character(lev[j]))) {
          ord[i] <- j
        }
      }
    }
    cord <- cord[order(ord[ord >= 0]), ]
    cord$label <- GGally:::str_c(cord[[colorCol]], ": ", cord$ggally_cor)
    xmin <- min(xVal, na.rm = TRUE)
    xmax <- max(xVal, na.rm = TRUE)
    xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
                  (xmax - xmin))
    ymin <- min(yVal, na.rm = TRUE)
    ymax <- max(yVal, na.rm = TRUE)
    yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * 
                  (ymax - ymin))
    p <- ggally_text(label = GGally:::str_c("Corr: ", signif(cor_fn(xVal, 
                                                                    yVal), sgnf)), mapping = mapping, xP = 0.5, yP = 0.9, 
                     xrange = xrange, yrange = yrange, color = "black", 
                     ...) + theme(legend.position = "none")
    xPos <- rep(alignPercent, nrow(cord)) * diff(xrange) + 
      min(xrange, na.rm = TRUE)
    yPos <- seq(from = 0.9, to = 0.2, length.out = nrow(cord) + 
                  1)
    yPos <- yPos * diff(yrange) + min(yrange, na.rm = TRUE)
    yPos <- yPos[-1]
    cordf <- data.frame(xPos = xPos, yPos = yPos, labelp = cord$label)
    cordf$labelp <- factor(cordf$labelp, levels = cordf$labelp)
    p <- p + geom_text(data = cordf, aes(x = xPos, y = yPos, 
                                         label = labelp, color = labelp), hjust = 1, ...)
    p
  }  else {
    xmin <- min(xVal, na.rm = TRUE)
    xmax <- max(xVal, na.rm = TRUE)
    xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
                  (xmax - xmin))
    ymin <- min(yVal, na.rm = TRUE)
    ymax <- max(yVal, na.rm = TRUE)
    yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * 
                  (ymax - ymin))
    p <- ggally_text(label = paste("Corr:\n", signif(cor_fn(xVal, 
                                                            yVal), sgnf), sep = "", collapse = ""), mapping, xP = 0.5, 
                     yP = 0.5, xrange = xrange, yrange = yrange, ...) + 
      theme(legend.position = "none")
    p
  }
}

png("Pairs_airquality.png", height = 6, width = 8, units = "in", 
    res = 500)
mehra.eda %>% sample_n(5e5) %>% 
  select(CVD60, one_of(pollutants)) %>%
  mutate(CVD60 = ifelse(CVD60 > 0.4, NA, CVD60), #Exclude the most extreme observations to
         CO= ifelse(CO > 8, NA, CO),             #improve readability of figure
         NO2 = ifelse(NO2 > 200, NA, NO2),
         PM10 = ifelse(PM10 > 150, NA, PM10),
         O3 = ifelse(O3 > 150, NA, O3),
         PM2.5 = ifelse(PM2.5 > 90, NA, PM2.5),
         SO2 = ifelse(SO2 > 150, NA, SO2)) %>%
  ggpairs(upper = list(continuous = mycor),
    lower = list(continuous = function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) +
      geom_point(..., alpha = 1/10, shape = ".")
  })) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marginal distributions and pairwise dependence of pollutant concentrations",
          subtitle = paste("Pollutant concentrations are substantially correlated with each other and",
                           "only weakly correlated with mortality")) +
  theme.png 
dev.off()

png("Pairs_weather.png", height = 6, width = 8, units = "in", 
    res = 500)
mehra.eda %>% sample_n(1e5) %>% 
  select(CVD60, one_of(weather)) %>%
  ggpairs(upper = list(continuous = mycor),
          lower = list(continuous = function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) +
      geom_point(..., alpha = 1/20, shape = ".")
  })) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marginal distributions and pairwise dependence of weather variables",
          subtitle = paste("Weather variables and CVD60 correlate",
                           "in two groups: (CVD60 + SSR + T2M) and (BLH + WS + TP + SSR)")) +
  theme.png 
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
      lab = round(cor, digits = 2)),
      mapping = aes(x = x, y = y, label = lab),
      hjust = 1, vjust = 1,
      size = 3,
      inherit.aes = FALSE)
}


png("Pairs_pol_weather.png", height = 6, width = 8, units = "in", 
    res = 500)
mehra.eda %>% sample_n(5e5) %>% 
  mutate(CO = ifelse(CO > 8, NA, CO),
         NO2 = ifelse(NO2 > 200, NA, NO2),
         PM10 = ifelse(PM10 > 150, NA, PM10),
         O3 = ifelse(O3 > 150, NA, O3),
         PM2.5 = ifelse(PM2.5 > 90, NA, PM2.5),
         SO2 = ifelse(SO2 > 150, NA, SO2)) %>%
  ggduo(columnsX = 18:23, columnsY = 11:16, 
        types =list(continuous = pointsWithCor)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Pairwise dependencies between pollutants and weather variables",
          subtitle = paste("Between-group dependency seems to be mediated by",
                           "pollutants O3 and PMx and weather variables SSR and BLH")) +
  theme.png 
dev.off()


  

#' ** STRUCTURE LEARNING AND HYPERPARAMETER TUNING **
invisible(lapply(c("Scoring.R", "Fit DAG.R", "Experiments.R", 
                   "MEHRA functions.R"), source))

if(!exists("mehra.daily")) {mehra.daily <- readRDS("mehra_daily.RDS")}

#' Generate and store training/validation/testing split
n <- nrow(mehra.daily)
test.indices <- sample(1:n, round(0.2*n))
val.indices <- sample((1:n)[-test.indices], round(0.2*n))
train.indices <- (1:n)[-c(test.indices, val.indices)]

save(train.indices, val.indices, test.indices, file = "MEHRA_split.rdata")

load("MEHRA_split.rdata")


training <- mehra.daily[train.indices,]
validation <- mehra.daily[val.indices,]
testing <- mehra.daily[test.indices,]

#' Configure hyperparameter tuning procedure
penalties <- c("0.10", "0.25", "0.40")
thresholds <- c(0.50, 0.65, 0.80, 0.95, 0.99)
no.boot <- 200
particles <- 500
obs.perParam <- 100

bl <-  {data.frame(
  "from" = c(rep("Region", 9), rep("Zone", 9), rep("Type", 9),
             rep("Year", 9), rep("Season", 9), rep("Month", 9),
             rep("Latitude", 9), rep("Longitude", 9), rep("Altitude", 9),
             rep("Day", 9), rep("T2M", 10), rep("WS", 10),
             rep("WD", 10), rep("TP", 10), rep("BLH", 10),
             rep("SSR", 10), rep("NO2", 10), rep("SO2", 10), 
             rep("CO", 10), rep("O3", 10), rep("PM10", 10),
             rep("PM2.5", 10), rep("CVD60", 22), rep("CVD.ly", 23)),
  
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
           "T2M", "WS", "WD", "TP", "BLH", "SSR", "NO2", 
           "O3", "SO2", "CO", "PM10", "PM2.5", "Day",
           
           "Region", "Zone", "Type", "Year", "Season",
           "Month", "Latitude", "Longitude", "Altitude",
           "T2M", "WS", "WD", "TP", "BLH", "SSR", "NO2", 
           "O3", "SO2", "CO", "PM10", "PM2.5", "Day", "CVD60"),
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


#' * Generate plot comparing RMSE for different hyperparameter configurations *
RMSE <- c(RMSE10, RMSE25, RMSE40)

input <- expand.grid(thresholds, penalties)

mehra.tuning <- cbind(input, RMSE)
names(mehra.tuning) <- c("Threshold", "Penalty", "RMSE")

png("RMSE_plot.png", height = 2, width = 8, units = "in", 
    res = 500)
ggplot(mehra.tuning, aes(Threshold, RMSE)) + 
  geom_line(aes(color = Penalty, linetype = Penalty)) +
  geom_point(aes(color = Penalty)) +
  labs(y = "Av. norm. RMSE", x = expression(Arc~inclusion~threshold~(gamma))) +
  ggtitle("Average normalized RMSE for different hyperparameter configurations",
          subtitle = expression(paste("Penalties ", alpha == 0.10, " and ", alpha == 0.25,
                      " outperform ", alpha == 0.40, 
                      " and give models with comparable performance"))) +
  scale_color_discrete(name = expression(Penalty~(alpha))) +
  scale_linetype_discrete(name = expression(Penalty~(alpha))) +
  theme.png
dev.off()


#' * Generate ECDF figure *
bootSamples10 <- readRDS("bootSamples10.rds")
bootSamples25 <- readRDS("bootSamples25.rds")

strength10 <- custom.strength(bootSamples10, names(mehra.daily)) %>% 
  anti_join(bl, by = c("from", "to"))

strength25 <- custom.strength(bootSamples25, names(mehra.daily)) %>% 
  anti_join(bl, by = c("from", "to"))

ECDF.df <- rbind(strength10,strength25) %>% 
  mutate(penalty = rep(c("alpha:~0.10", "alpha:~0.25"), each = n()/2)) %>% 
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

threshold.df <- data.frame(penalty = rep(c("alpha:~0.10", "alpha:~0.25"), 2),
                           type = rep(c("auto", "custom"), each = 2),
                           threshold = thresholds,
                           y = y.thresholds)

png("ECDF.png", height = 4, width = 8, units = "in", 
    res = 500)
ggplot(threshold.df) +
  facet_grid(~penalty, labeller = label_parsed) +
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
  scale_linetype_manual(name = "Threshold\ntype", values = c(3,2)) +
  ylim(c(0,1)) +
  labs(x = expression(Threshold~(gamma)), 
       y = "Proportion of arcs excluded") + 
  ggtitle("Empirical CDF of arc strengths",
          subtitle = expression(italic(paste("Configuration (", alpha,",",gamma,") = (0.10, 0.80) ",
                                      "leads to a simpler network than (",
                                      alpha,",",gamma,") = (0.25, 0.99)", sep = "")))) +
  theme.png
dev.off()


#' * Train model with optimal hyperparameter configuration *
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


#' ** RESULTS AND VALIDATION**

#' Generate figure of final DAG
bootSamples.best <- readRDS(paste("Bootsamples_full", penalty.best, ".rds", sep = ""))
mehra.fit <- readRDS(paste("Mehra_fit_full_0.10_0.8.rds", sep = ""))

mehra.dag <- bn.net(mehra.fit)

svg("Mehra_dag.svg", height = 6, width = 8)
graphviz.plot(mehra.dag)
dev.off()


#' Find number of arcs and number of parameters
narcs(mehra.fit)
nparams(mehra.fit)

#' Make in- and out-of-sample predictions
prednames <- names(mehra.daily)[8:23]
mehra.pred <- lapply(prednames, computePredictions, 
                     dat = testing, dag.fit = mehra.fit,
                     particles = particles, parallel = TRUE)
saveRDS(mehra.pred, paste("Mehra_pred_full_", penalty.best,"_", threshold.best,
                          ".rds", sep = ""))

mehra.trainPred <- lapply(prednames, computePredictions, 
                          dat = rbind(training, validation), dag.fit = mehra.fit,
                          particles = particles, parallel = TRUE,
                          train = TRUE)

mehra.pred <- readRDS(paste("Mehra_pred_full_0.10_0.8.rds", sep = ""))


#' * Generate RMSE figure *
RMSE.train.best <- sapply(mehra.trainPred, computeRMSE)
RMSE.best <- sapply(mehra.pred, computeRMSE)[-17]

RMSE.train.vit <- c(0.21, 0.25, 0.10, 0.05, 0.04, 0.14, 0.03, 0.07, 0.11, 0.04, 0.01,
0.05, 0.01, 0.00, 0.03, 0.02)
RMSE.test.vit <- c(0.49,	0.40,	0.94,	0.09,	0.19,	0.34,	0.06,	0.10,	0.13,	0.27,	0.08,
                   0.12,	0.03,	0.03,	0.06,	0.03)

prednames.caps <- c(prednames[1:3],toupper(prednames)[4:16])

RMSE.df <- rbind(data.frame(Analysis = "Current model", Name = prednames.caps,
                 test = RMSE.best, train = RMSE.train.best),
      data.frame(Analysis = "Model of Vitolo et al.", Name = prednames.caps,
                 test = RMSE.test.vit, train = RMSE.train.vit)) %>%
  mutate(Diff = ifelse(test - train > 0, test-train, 0),
         train = ifelse(train > test, test, train)) %>% 
  gather(key = "Type", value = "Value", -Analysis, -Name) %>% 
  filter(Name != "CVD.ly")

RMSE.df$Analysis <- factor(RMSE.df$Analysis, levels = c("Model of Vitolo et al.", "Current model"))

saveRDS(RMSE.df, "RMSE_df.rds")

png("RMSEcomparison.png", height = 6, width = 8, units = "in",
    res = 500)
ggplot(RMSE.df, aes(reorder(Name, Value), Value)) +
  facet_grid(~Analysis) +
  geom_col(data = filter(RMSE.df, Type != "test"), aes(fill = Type)) +
  coord_flip() +
  geom_text(data = filter(RMSE.df, Type == "test"), 
            aes(x = Name, y = -0.05, label = round(Value, 2)), hjust = 1,
            size = 3) +
  ylim(c(-0.15, 1.1)) +
  theme(axis.title.y = element_blank()) +
  labs(y = "Normalized RMSE") +
  scale_fill_discrete(name = "Predictions", labels = c("Out-of-sample", "In-sample")) + 
  ggtitle("Comparison of predictive performance of current model with original model",
          subtitle = paste("Our model appears to overfit less than the model of Vitolo et al.")) +
  theme.png +
  guides(fill = guide_legend(reverse = TRUE))
dev.off()


#' * Produce QQ-plots *
mehra.pred.named <- rbindlist(mapply(function(df, name){df$name = name; return(df)},
                           df = mehra.pred[1:16], name = prednames.caps, SIMPLIFY = FALSE))
                                                                                                                                                                                                          

mehra.trainPred.named <- rbindlist(mapply(function(df, name){df$name = name; return(df)},
                             df = mehra.trainPred, name = prednames.caps, SIMPLIFY = FALSE))

mehra.trainPred.named <- mehra.trainPred.named %>% 
  mutate(Longitude = rep(training.full$Longitude, 16),
         Latitude = rep(training.full$Latitude, 16),
         resid = actuals - fits)


png("qqplots.png", height = 5, width = 8, units = "in",
    res = 500)
mehra.trainPred.named %>% 
  filter(name != "CVD.ly") %>% 
  group_by(name) %>% 
  sample_frac(0.05) %>% 
  mutate(residuals = actuals - fits) %>%
  ggplot(aes(sample = residuals)) +
  facet_wrap(~name, scales = "free") +
  geom_qq() +
  geom_qq_line() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "Normal quantiles", y = "Residual quantiles") + 
  ggtitle("Normal QQ-plots for residuals",
          subtitle = paste("Residuals show clear signs of non-normality")) +
  theme.png
dev.off()

mehra.dag <- bn.net(mehra.fit)



#' * Generate correlation plots *
mehra.residuals <- mehra.trainPred.named %>% 
  mutate(residuals = actuals - fits) %>% 
  select(name, residuals) %>%
  group_by(name) %>% 
  mutate(index = 1:n()) %>% 
  spread(key = 1, value = 2) %>% 
  select(-index)

resid.cor <- mehra.residuals %>% 
  select("CVD60", one_of(weather), one_of(pollutants)) %>% 
  cor(.,use = "complete.obs") %>% 
  as.data.frame %>% 
  mutate(name1 = rownames(.)) %>% 
  gather(key = "name2", value = "corr", -name1) %>%
  mutate(type = "Residuals")

feature.cor <- mehra.daily %>% 
  select("CVD60", one_of(weather), one_of(pollutants)) %>%
  cor(.,use = "complete.obs") %>% 
  as.data.frame %>% 
  mutate(name1 = rownames(.)) %>% 
  gather(key = "name2", value = "corr", -name1) %>%
  mutate(type = "Original features")

cor.df <- rbind(resid.cor, feature.cor) %>% 
  mutate(corr = ifelse(name1 == name2, 0, corr))
cor.df$name1 <- factor(cor.df$name1, levels = c("CVD60", pollutants, weather))
cor.df$name2 <- factor(cor.df$name2, levels = c("CVD60", pollutants, weather))

png("corrComparison.png", height = 4, width = 8, units = "in",
    res = 600)
ggplot(cor.df, aes(name1, name2)) +
  facet_wrap(~type) +
  geom_raster(aes(fill = abs(corr))) +
    scale_fill_gradient(name = "Absolute\ncorrelation",
                        low = "white", high = "black", limits = c(0,1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL) + 
  ggtitle("Correlation of original features and of residuals",
          subtitle = paste("The model captures most of the linear dependence between variables,",
                           "especially for pollutants")) +
  theme.png
dev.off()

#' * Generate autocorrelation plots *
autocor.resid <- mehra.trainPred.named %>% filter(abs(Latitude - chosenLat) < 1e-5,
                                                  abs(Longitude + chosenLon) < 1e-4) %>%
  select(-fits, - actuals, - Latitude, - Longitude) %>% 
  group_by(name) %>% 
  mutate(id = 1:n()) %>% 
  spread(key = name, value = resid) %>%
  select(one_of(pollutants, weather)) %>% 
  lapply(.,function(x){acf(x, lag.max = 40, type = "correlation", plot = FALSE, na.action = na.pass)$acf}) %>%
  data.frame() %>% 
  mutate(lag = 0:40) %>%
  gather(key = "key", value = "value", -lag) %>% 
  mutate(type = "Residuals")


autocor.daily$type <- "Original features"
autocor2 <- rbind(autocor.daily, autocor.resid)
autocor2$type <- factor(autocor2$type, levels = c("Original features", "Residuals"))
autocor2$key <- factor(autocor2$key, levels = c(pollutants, weather))

png("Autocorrelation_resid.png", height = 5, width = 8, units = "in", 
    res = 500)
autocor2 %>% ggplot(aes(x = lag, y = key)) +
  facet_wrap(~type) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(name = "Auto-\ncorrelation", low = "white", high = "black") +
  labs(x = "Lag", y = NULL)+
  ggtitle("Autocorrelation of original features and of residuals",
          subtitle = "Not all autocorrelation is captured by the model") +
  theme.png
dev.off()

#' ** INFERENCE **

#' Check d-separations
sapply(pollutants, dsep, bn = mehra.dag, y = "Type")
sapply(pollutants, dsep, bn = mehra.dag, y = "Year")
sapply(pollutants, dsep, bn = mehra.dag, y = "Month")

#' According to the network, mortality is fully independent of:
names(training)[sapply(names(training), dsep, bn = mehra.dag, y = "CVD60")]
names(training)[sapply(names(training), dsep, bn = mehra.dag, y = "Year")]
names(training)[sapply(names(training), dsep, bn = mehra.dag, y = "Month")]

#' No pollutant is d-separated from any other pollutant
input <- expand.grid(pollutants, pollutants, stringsAsFactors = FALSE)
ouTPut <- mapply(dsep, x = input[[1]], y = input[[2]], MoreArgs = list(bn = mehra.net))
cbind(input, ouTPut)

#' SSR is d-separated from WD and BLH
input <- expand.grid(weather, weather, stringsAsFactors = FALSE)
ouTPut <- mapply(dsep, x = input[[1]], y = input[[2]], MoreArgs = list(bn = mehra.net))
cbind(input, ouTPut)[ouTPut,]

dsep("SSR", "T2M", "Month", bn = mehra.dag)

#' SSR is d-separated from PM10 and PM2.5
input <- expand.grid(weather, pollutants, stringsAsFactors = FALSE)
ouTPut <- mapply(dsep, x = input[[1]], y = input[[2]], MoreArgs = list(bn = mehra.net))
cbind(input, ouTPut)[ouTPut,]



#' * Generate figure with simulated mortality *
cl <- makeCluster(detectCores()-1)

samples.year <- cpdist(mehra.fit, nodes = c("Year", "Month", "CVD60", tolower(pollutants)),
                       evidence = list(Year = c(1981, 2014)),
                       method = "lw", cl = cl, n = 1000000)
stopCluster(cl)

#' No need to use weights as evidence is at a root of the network
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
             aes(x = Date, y = CVD60), alpha = 1, shape = ".",
             col = "darkgrey") +
  geom_line(size = 0.7) +
  geom_line(aes(y = mean + 2*sd), linetype = 2, size = 0.7) +
  geom_line(aes(y = mean - 2*sd), linetype = 2, size = 0.7) +
  ylim(c(0, 0.4)) +
  labs(y = "CVD60") + 
  ggtitle("Simulated mortality compared with actual observations",
          subtitle = "The model captures the temporal pattern in mortality well") +
  theme.png
dev.off()