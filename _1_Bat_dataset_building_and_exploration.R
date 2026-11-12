
# 2025 bat data analysis

# Libraries: ----
{library(lme4) # for glmer
library(lattice)
library(ggplot2)
library(DHARMa)
library(dplyr)
library(rgl)
library(mgcv)
library(performance) # r2()
library(emmeans) 
library(multcomp) # emmeans letters
library(multcompView)
library(ggeffects)
setwd("C:/Users/valavale/Documents/PhD/ANALYSIS/CODE/HighlandStats/GAM, intro to R etc/AllData")
source("HighstatLibV14.R")
library(mgcViz)
library(glmmTMB)
library(gamm4)
library(gratia)  # nice plotting for mgcv gam objects
library(lubridate)
library(tibble)}

# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#binomial-data
#https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html

#* 2.1 Bring in data ----


#setwd("D:/Batlab Assistant/LIGHT POLLUTION/R LAMPBATS/saving tables 0811")
#setwd("//ad.helsinki.fi/home/v/valavale/Desktop/backup bats and light")
setwd("C:/Users/valavale/Documents/backup bats and light")
Data_10min<-read.csv("Data_10min.csv")
Data10_with_E5L<-read.csv("Data10_with_E5L.csv")

(head(Data_10min, 2))

Dat<-Data_10min
Dat_all<-Data10_with_E5L

table(Dat_all$Device, Dat_all$Year)

nrow(Dat_all[which(Dat_all[Data_10min$Year=="2022",]$Device=="E5L"),])

Dat_all_2022<-Dat_all[which(Dat_all$Year=="2022"),]
Dat_all_2023<-Dat_all[which(Dat_all$Year=="2023"),]
# unique(Dat_all_2022$Year)
table2022<-table(Dat_all_2022$Device, Dat_all_2022$Ept)
table2023<-table(Dat_all_2023$Device, Dat_all_2023$Ept)
table2223<-cbind(table2022, table2023[-c(13:16),])

#* Section 3: Prepare the data----
#* Section 3.1: Create variables----

#' Separate different loations of H3 of differnet years
Dat_all[which(Dat_all$Device == "H3L" & Dat_all$Year == "2022"),]$Device<-"H3L_2022"
Dat_all[which(Dat_all$Device == "H3D" & Dat_all$Year == "2022"),]$Device <-"H3D_2022"
Dat_all[which(Dat_all$Device == "H3L" & Dat_all$Year == "2023"),]$Device <-"H3L_2023"
Dat_all[which(Dat_all$Device == "H3D" & Dat_all$Year == "2023"),]$Device <-"H3D_2023"


#'Define categorical covariates.
Dat_all$fYear <- factor(Dat_all$Year)
Dat_all$fTreatment <- factor(Dat_all$Treatment)
Dat_all$fMonth <- factor(Dat_all$Month)
Dat_all$fLamp <- factor(Dat_all$Lamp)
Dat_all$fCity <- factor(Dat_all$City)
Dat_all$fDevice <- factor(Dat_all$Device)
Dat_all$fHour <- factor(Dat_all$tunti)
Dat_all$fHalfHour <- factor(Dat_all$HalfHour)

Dat_all$Humid <- as.numeric(Dat_all$Humid)
Dat_all$Clouds <- as.numeric(Dat_all$Clouds)
#Dat_all$Temp
#Dat_all$Lux

#unique(Dat_all$fDevice) # all good
#str(Dat_all)
#Dat<-Dat_all[!(Dat_all$Month == 'Jun'& Dat_all$Device=="E5L"& Dat_all$Year=="2023"),]



#* CREATE Time since sunset
# Only hours with bat activity were considered for downstream analyses (22:00-04:00 in June and 21:00-05:50 in August). 

Dat_all <- Dat_all %>%
  filter(
    (fMonth == "Jun" & (Time >= "22:00:00" | Time <= "04:00:00")) |
      (fMonth == "Aug" & (Time >= "21:00:00" | Time < "06:00:00"))
  )


# Time since sunset
# Here time 1 is the first 10 minute interval considered in the analysis, even if the sun would set a bit later
sunsets <- tibble(
  Month  = c("Jun", "Aug"),
  Sunset = hms(c("22:00:00", "21:00:00"))) 

Dat_all<-Dat_all %>%
  left_join(sunsets, by = "Month") %>%
  mutate(
    t = hms(Time),
    # unwrap 00:00–06:00 (inclusive) to the next day
    t_unwrapped = if_else(t <= hms("06:00:00"), t + hours(24), t),
    
    # Convert Period to Duration to avoid class mismatch
    minutes_since_sunset = as.numeric(as.duration(t_unwrapped - Sunset) / dminutes(1)),
    TimeSunset = floor(minutes_since_sunset / 10) + 1   
  )
# checking times to data
# str(Dat_all)
# (unique(Dat_all[which(Dat_all$fMonth=="Jun"&Dat_all$fYear=="2022"),]$Time))
# (unique(Dat_all[which(Dat_all$fMonth=="Jun"&Dat_all$fYear=="2023"),]$Time))
# (unique(Dat_all[which(Dat_all$fMonth=="Aug"&Dat_all$fYear=="2022"),]$Date))
# (unique(Dat_all[which(Dat_all$fMonth=="Aug"&Dat_all$fYear=="2023"),]$Date))
# 


head(Dat_all,2)
#* CREATE Temperature_Z

# Mean temperature
Dat_all$DateTime<-paste(Dat_all$Date, Dat_all$Time) 
Dat_all$DateTime<- as.POSIXlt(Dat_all$DateTime, format = "%Y-%m-%d %H:%M:%S")


temp_mean_10min <- Dat_all %>%
  #group_by(fCity) %>%
  group_by(fCity, DateTime) %>%
  summarise(
    Temp_mean = mean(Temp, na.rm = TRUE),
    n_devices = n()
  )

Dat_all <- Dat_all %>%
  left_join(temp_mean_10min, by =c("fCity", "DateTime"))


# SD check
# x<-Dat_all %>%
#   group_by(DateTime) %>%
#   summarise(
#     Sd_temp = sd(Temp, na.rm = TRUE)
#   )
#sort(x$Sd_temp, decreasing = T )

# Z-score standardised temperature
Dat_all$Temp_z <- as.numeric(scale(Dat_all$Temp))
Dat_all$Temp_mean_z <- as.numeric(scale(Dat_all$Temp_mean))

# Z-score standardised temperature
Dat_all$Lux_z <- as.numeric(scale(Dat_all$Lux))
#Dat_all$Lux_mean_z <- as.numeric(scale(Dat_all$Lux_mean))


#* Section 3.2: Create Monthly data ----
Dat_all_jun<-Dat_all[which(Dat_all$Month=="Jun"),]
Dat_all_aug<-Dat_all[which(Dat_all$Month=="Aug"),]




#  Section 3.3: Check for OUTLIERS? PROBLEMS WITH DEVICES E1L AND H3D - not anymore ----

# testausmielessä otetaan pois ne rivit, joissa on liian suuria arvoja:
#Dat_all_2<-(Dat_all[which(Dat_all$Temp<30),])
# 
# # Drops rows where Temp is non-NA and > 30; keeps NA
# Dat_all_2 <- Dat_all[Dat_all$Temp <= 30 | is.na(Dat_all$Temp), ]
# Dat_all<-Dat_all_2
# 
# Dat_noNa<-(Dat_all[which(Dat_all$fYear=="2022"),])
# Dat_noNa<-Dat_noNa[!is.na(Dat_noNa$Lux),]
# max(Dat_noNa$Temp)
# max(Dat_noNa$Lux)
# 
# dat_over30<-(Dat_all[which(Dat_all$Temp>30),])
# mean(Dat_noNa[which(Dat_noNa$fDevice=="E1L"),]$Lux)
# 
# boxplot(Dat_noNa$Temp~Dat_noNa$fDevice)
# boxplot(Dat_noNa[which(Dat_noNa$fDevice=="H3D"),]$Lux~Dat_noNa[which(Dat_noNa$fDevice=="H3D"),]$fHour)
# 
# max(Dat_noNa[which(Dat_noNa$fDevice=="E1L"),]$Lux)
# max(Dat_noNa[which(Dat_noNa$fDevice=="E1D"),]$Lux)
# 
# max(Dat_noNa[which(Dat_noNa$fDevice=="E1L"),]$Temp)
# max(Dat_noNa[which(Dat_noNa$fDevice=="E1D"),]$Temp)
# 
# # klo 20 puuttuu 2022
# unique(Dat_noNa$fHour)

# temperature
# ALL GOOD
# Dat_noNa<-(Dat_all[which(Dat_all$fYear=="2023" & Dat_all$fMonth=="Jun"),])
# plot(Dat_noNa$Temp)
# boxplot(Dat_noNa[which(Dat_noNa$fDevice=="H3D"),]$Temp~Dat_noNa[which(Dat_noNa$fDevice=="H3D"),]$TimeSunset)
# boxplot(Dat_noNa[which(Dat_noNa$fDevice=="E1D"),]$Temp~Dat_noNa[which(Dat_noNa$fDevice=="E1D"),]$TimeSunset)






# Section 4: Data exploration----

#' How many zeros and ones do we have in the response variable?
table(Dat_all$Myo)
table(Dat_all$Ept)
#' Much more zeros than ones.
#' 56081 versus 319 for Myo = 0.5%
#' 54714 versus 1686 for Ept = 3%
#' If your question is: I have
#' 900 zeros and 100 ones...is that good or bad...then our
#' answer is: Don't know. Do a simulation study to assess the
#' effect of this. -- See later: no zero inflation in models




#* Subsection 4.1: Outliers----

#' Use Cleveland dotplots to see if there are any outliers in the
#' non-categorical covariates.
Myvar <- c("Lux_z", "Temp_z", "Clouds", "Rain_mm", "Humid", "Wind")
Mydotplot(Dat_all[,Myvar])
#' Some of them have a some weird clusters of larger values (temp, lux for example)


#* Subsection 4.2: Missing values----

#' Missing values?
#str(Dat_all)
colSums(is.na(Dat_all)) 
#' Quite a bit in the logger data: Temp, TempAvg, Lux, LuxAvg 
#' And in Rain_mm for some reason, wind also a bit
#' There are na also in temp_mean for those dates (like 2023-06-28) where 
#'   the recording was still on but the HOBO had stopped from logging temp/lux


#* Subsection 4.3: Categorical variables----

# Bat presence versus Treatment / Lamp:
table(Dat_all$Ept, Dat_all$Treatment)
# dark   HPS   LED
# 0 27652 10699 16363
# 1   548   521   617
table(Dat_all$Ept, Dat_all$Treatment, Dat_all$fMonth)
# Dark Light
# 0 27652 27062
# 1   548  1138

table(Dat_all$Ept, Dat_all$fLamp, Dat_all$fMonth)
# = Aug
#  dark   HPS   LED
#   0 14969  6107  8993
#   1   151    13     7
# 
# , ,  = Jun
#      dark   HPS   LED
#   0 12961  5089  7962
#   1   119    11    18

table(Dat_all$Ept, Dat_all$fTreatment)
#      Dark Light
# 0 27652 27062
# 1   548  1138
table(Dat_all$Myo, Dat_all$fTreatment)
# Dark Light
# 0 27930 28151
# 1   270    49

#' How many years?
table(Dat_all$fYear)
#' 2022  2023 
#' 21840 34560 

#' How many Treatment?
table(Dat_all$fTreatment)
#'  Dark Light 
#'  28200 28200 

#' How many Months?
table(Dat_all$fMonth)
#'  Aug   Jun 
#'  30240 26160  

#' How many Lamp?
table(Dat_all$fLamp)
#'   dark   HPS   LED 
#'   28200 11220 16980

#' How many City?
table(Dat_all$fCity)
#'      E     H 
#'      36144 20256 

#' How many Device?
table(Dat_all$fDevice)
#'  All 24 levels have 2532 or 1440 observations

# How many fHour?   
table(Dat_all$fHour)
#     0    1    2    3    4    5   20   21   22   23 
# 6336 6336 6336 6336 5136 4896 3456 4896 6336 6336 

#* Subsection 4.4: Collinearity covariates----

Dat_all_jun$logLux<-log1p(Dat_all_jun$Lux)
Dat_all_aug$logLux<-log1p(Dat_all_aug$Lux)

#' Make a pairplot.
Myvar <- c("logLux", "Temp_mean_z", "Clouds", "Rain_mm", "Humid", "Wind", "TimeSunset")
Mypairs(Dat_all_jun[,Myvar]) # Temp and timesunset correlated in June 0.57. 
Mypairs(Dat_all_aug[,Myvar])  # humid and wind correlated in aug
Mypairs(Dat_all[,Myvar])  


#' Is there  collinearity with factor Year and the continuous 
#' covariates? To answer this question, make conditional boxplots.
par(mfrow = c(2,3))
boxplot(Lux_z ~ factor(Treatment), data = Dat_all_jun)  #' 
boxplot(Temp_z ~ factor(Treatment), data = Dat_all_jun) #' 
boxplot(Clouds ~ factor(Treatment), data = Dat_all_jun)  #' 
boxplot(Rain_mm ~ factor(Treatment), data = Dat_all_jun)  #' 
boxplot(Humid ~ factor(Treatment), data = Dat_all_jun)    #' 
boxplot(Wind ~ factor(Treatment), data = Dat_all_jun)    #' 

# What about fHour? 
boxplot(LuxAvg ~ fHour, data = Dat_all_aug)  #' 
boxplot(TempAvg ~ fHour, data = Dat_all_aug)  #' 
boxplot(Clouds ~ fHour, data = Dat_all_aug)  #' 
boxplot(Rain_mm ~ fHour, data = Dat_all_aug)  #' 

#' Did this for all factors: City, Treatment, Lamp, Month, Year, Hour (see wd)
#' No problems with Treatment or Lamp, but some possible problems
#' with Year and Month (especially humid, Wind and Clouds),
#' and City (problem with Humid?)
#' Some problems with fHour and LuxAvg, and TempAvg (maybe less?). If a model 
#' is with fHour, maybe do at least without LuxAvg (think about temp)



#* Subsection 4.5: Relationships----

#' Plot RESPONSE of Ept versus X
ggplot(data = Dat_all,
       aes(x = fHour, 
           y = Ept)) +
  geom_point(alpha = 0.5) +
  # labs(x = "x", y = "Ept") +
  theme(text = element_text(size = 15)) +
  geom_smooth(se = FALSE, span = 1) #' Forcing the line to be smoother
#' Something going on somewhere, see plots in wd


#' Plot RESPONSE of Myo versus X
ggplot(data = Dat_all,
       aes(x = fHour, 
           y = Myo)) +
  geom_point(alpha = 0.5) +
  #labs(x = "x", y = "y") +
  theme(text = element_text(size = 15)) +
  #  geom_smooth(se = FALSE, span = 1) #' Forcing the line to be smoother
  geom_smooth(method = "gam", method.args = list(family = "binomial"), se = FALSE)

#' It would be good to tell the geom_smooth function that it should
#' use a Bernoulli smoother so that the line stays between 0 and 1. 
#' You can do that via:
#' geom_smooth(method = "gam", method.args = list(family = "binomial"), se = FALSE)
#' --> did it for Myo, see Myo plot: different geom_smooth()'s



#* Subsection 4.6: Conclusions data exploration----

#' 1. The response variable is binary.
#' 2. We should drop some covariates due to collinearity. Maybe Clouds and Humidity.
#'    We could also drop Rain_mm because there is barely any variation there.
#' 3. With a Bernoulli GLM, scatterplots should show
#'    sigmoid-type patterns. But some covariates maybe show a 
#'    non-linear relationship. A GAM may be a better approach
#'    than a GLM. Model validation will tell.
#' 4. A lot of zeros compared to ones - maybe can try clog-log, or simulation study?
#' 5. Many explanatory variables are telling about the same thing: Lamp / Treatment / LuxAvg,
#'    However, LuxAvg does not seem to be collinear with Lamp or Treatment. 
#'    






# SECTION 5. Feeding buzzes ----

# Read data
setwd("C:/Users/valavale/Documents/backup bats and light")
ept22<-read.csv("C:/Users/valavale/Documents/backup bats and light/Ept2022.csv")
ept23<-read.csv("C:/Users/valavale/Documents/backup bats and light/Ept2023.csv")


# Data from 2022
# There are no feeding buzzes counted for dark / unlit
HPS  <- c("E1L", "E2L", "E3L", "E4L", "E6L")
LED  <- c("E5L", "H1L", "H2L", "H3L", "H4L")
dark <- c("E1D", "E2D", "E3D", "E4D", "E5D", "E6D", 
          "H1D", "H2D", "H3D", "H4D")
ept22$treatment <- ifelse(ept22$Device %in% LED,  "LED",
                          ifelse(ept22$Device %in% HPS,  "HPS",
                                 ifelse(ept22$Device %in% dark, "Dark", NA)))
with(ept22, table(Device, treatment))


# Data from 2023
# There are all feeding buzzes calculated for all treatments
HPS  <- c("E3L", "E4L", "E6L", "E8L")
LED  <- c("E1L", "E2L", "E5L", "E7L", "H1L", "H2L", "H3L", "H4L")
dark <- c("E1D", "E2D", "E3D", "E4D", "E5D", "E6D", "E7D", "E8D",
          "H1D", "H2D", "H3D", "H4D")
ept23$treatment <- ifelse(ept23$Device %in% LED,  "LED",
                          ifelse(ept23$Device %in% HPS,  "HPS",
                                 ifelse(ept23$Device %in% dark, "Dark", NA)))
with(ept23, table(Device, treatment))

# Making a dataset of both years but without dark treatment
ept_nodark_bothyears<-rbind(ept22, ept23)
ept_nodark_bothyears <- droplevels(
  subset(ept_nodark_bothyears, treatment != "Dark")
)

# No effect of treatments HPS and LED when modeling both years.
mod <- glmer(FeedingBuzz ~ treatment + (1|Device),
             data = ept_nodark_bothyears,
             family = binomial)
summary(mod)



# JUST 2023 and also dark sites
# Here modeling all treatments for year 2023
mod <- glmer(FeedingBuzz ~ treatment + (1|Device),
             data = ept23,
             family = binomial)

summary(mod) # there is no effect of treatment to feeding buzzes. 


