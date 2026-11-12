
#' Inspired by 10.7717/peerj.6876  
#'        Hierarchical generalized additive models in ecology: an introduction with mgcv, Pedersen et al 2019
#' Tried first to run a GS model from the paper but this was not possible with my data (too long run), therefore i have this more simple gam
#' WHY WE USE bam()?
#'      - This is because the size of the dataset is very large (<48k rows)
#' See  https://stackoverflow.com/questions/71926851/different-estimates-between-bam-and-gam-model-mgcv-and-interaction-term-estima
#' See separate code for myotis bats
#' Note: dark = unlit

# libraries: 
library(dplyr)
library(ggplot2)
library(DHARMa)
library(gratia)
library(emmeans)
library(mgcv)




# A) LAMP MODEL EPT -----
Dat_all$fNight<-factor(interaction(Dat_all$Device, Dat_all$Date_12, drop = TRUE))
#table(Dat_all$Night)[1:5]
#str(Dat_all)
#bat_test <- Dat_all[Dat_all$Date_12 %in% sample(unique(Dat_all$Date_12), 20), ]

Dat_all <- Dat_all[order(Dat_all$fNight, Dat_all$TimeSunset), ]
Dat_all$AR.start <- c(TRUE, Dat_all$fNight[-1] != Dat_all$fNight[-nrow(Dat_all)])

#run model for rho
bat_mod <- bam( Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_all, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res <- residuals(bat_mod, type = "pearson")
rho_est <- acf(res, plot = FALSE)$acf[2]
rho_est

# Correlation check with temperature
cor(Dat_all$Temp_mean_z, Dat_all$TimeSunset, use = "complete.obs") # -0.3085
cor(Dat_all_jun$Temp_mean_z, Dat_all_jun$TimeSunset, use = "complete.obs") #  -0.5693471


# Model - GAM ----
set.seed(123)
bat_mod_AR1 <- bam(
  Ept ~ 
    fLamp * fMonth + fYear + fCity +
    s(TimeSunset, by = fLamp, k = 10) + # treatment specific time-relationship
    s(Temp_mean_z) + # temperature, collinearity with time -0.31
    s(fNight, bs="re") + # per night variability
    s(fDevice, bs = "re"),
  family = binomial,
  data = Dat_all,
  method = "fREML",
  rho = rho_est,
  AR.start = Dat_all$AR.start,
  discrete = TRUE
)

summary(bat_mod_AR1)
gam.check(bat_mod_AR1) # this is ok, k indexes are not low or significant


# Autocorrelation
acf(residuals(bat_mod_AR1, type = "pearson")) ; acf(residuals(bat_mod_AR1, type = "pearson"), plot = FALSE)$acf[2] # 0.06 - 0.11 is ok for the second bar!
res2 <- residuals(bat_mod_AR1, type = "pearson")

plot(Dat_all$TimeSunset, res2, pch = 16, cex = 0.3) # doesn't work if temp is part of model
abline(h = 0, col = "red")
acf(res2)

# DHARMa
sim_res <- simulateResiduals(bat_mod_AR1, n = 1000)
plot(sim_res)
#testDispersion(sim_res)
testZeroInflation(sim_res) #ok

plot(
  fitted(bat_mod_AR1),
  Dat_all$Ept,
  pch = 16,
  cex = 0.3
)

# 'multicollinearity' 
concurvity(bat_mod_AR1)


# effect
draw(bat_mod_AR1)
smooth_estimates(bat_mod_AR1)

# COMPARISON BETWEEN SMOOTHERS----
#' Here we check four combination of different smoothers / variables to add to the model
#' Including night-to-night variation improved model by a lot,
#' Including separate TimeSunset smooth per treatment improved model a bit (but enough)
#' Therefore, both of these additions are included in the final model (mod3 chosen as the best model based on AIC)

# ***************************************************************
# MODELS:

mod3 <- bam(
  Ept ~ 
    fLamp * fMonth + fYear + fCity +  # in all
    s(TimeSunset, by = fLamp, k = 10) + # in models 2 & 3
   # s(TimeSunset, k = 10) + # in models 1 & 4
     s(Temp_mean_z) + # in all
    s(fNight, bs="re") + # in models 3 & 4 
    s(fDevice, bs = "re"), # in all, removing this worsens aic
  family = binomial,
  data = Dat_all,
  method = "fREML",
  rho = rho_est,
  AR.start = Dat_all$AR.start,
  discrete = TRUE
)
# > mod1
# Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, k = 10) + 
#   s(Temp_mean_z) + s(fDevice, bs = "re")
# fREML score: 4562.792     
# 
# > mod2
# Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, by = fLamp, 
#        k = 10) + s(Temp_mean_z) + s(fDevice, bs = "re")
# fREML score: 4570.01   
# 
# > mod3
# Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, by = fLamp, 
#         k = 10) + s(Temp_mean_z) + s(fNight, bs = "re") + s(fDevice,                                                                                        bs = "re")
# fREML score: 4401.592     
# 
# > mod4
# Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, k = 10) + 
#   s(Temp_mean_z) + s(fNight, bs = "re") + s(fDevice, bs = "re")
# fREML score: 4394.699     

# ***************************************************************
# Different smoothers
AIC(mod1, mod2,mod3, mod4)
# df      AIC
# mod1  46.81363 8996.431 
# mod2  61.78839 8991.684
# mod3 291.62509 8586.042 
# mod4 277.62308 8591.733 
# 
# Difference between mod3/4 and mod1/2: 
#   mod3 and mod4 compared to m1/m2: adding night-to-night variation improved model a lot
# Difference between mod 3 and mod 4: 
#   mod3: separate TimeSunset smooths for each lamp, mod4: one TimeSunset smooth. AIC difference 5.7.
# That is modest but real support for treatment-specific nightly activity patterns.
# Modest evidence that the shape of bat activity through the night differs among Dark, HPS and LED sites.
# 
# --> Going forward with model mod3

# ***************************************************************
# Temperature
AIC(mod3, mod6)  
#            df      AIC
# mod3 291.6251 8586.042
# mod6 320.4007 9311.157 # ΔAIC ≈ 725 is huge!
# > mod6
# Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, by = fLamp, 
#           k = 10) + s(fNight, bs = "re") + s(fDevice, bs = "re")

#' Temperature and time since sunset were moderately correlated (overall r ≈ -0.31; June r ≈ -0.57), 
#' but temperature was kept in the model because adding it improved model fit substantially.
#' 
#' Checking also if correlation between temperature and TimeSunset is too much:
#' AIC seems to improve substantially after including temperature in the model.
#' Even at moderate correlation (-0.57) temperature seems to explain variation in the data
#' that timeSunset can't explain alone.
#' Additionally concurvity showed elevated values (not too much), and bam.check was passed.
#' --> We keep temperature in the model.
#' 


# ***************************************************************
summary(mod3)
# 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, by = fLamp, 
#                                            k = 10) + s(Temp_mean_z) + s(fNight, bs = "re") + s(fDevice,                                                                                                bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -5.9277     0.5039 -11.762  < 2e-16 ***
#   fLampHPS             1.2380     0.6313   1.961 0.049860 *  
#   fLampLED             0.7879     0.6256   1.259 0.207867    
# fMonthJun            1.1734     0.1671   7.023 2.17e-12 ***
#   fYear2023            0.2715     0.1160   2.340 0.019275 *  
#   fCityH              -1.1314     0.6161  -1.836 0.066292 .  
# fLampHPS:fMonthJun  -0.3095     0.2298  -1.347 0.178106    
# fLampLED:fMonthJun   0.8954     0.2443   3.665 0.000248 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf   Ref.df  Chi.sq p-value    
# s(TimeSunset):fLampdark   7.017    8.020  162.30  <2e-16 ***
#   s(TimeSunset):fLampHPS    6.918    7.906  154.91  <2e-16 ***
#   s(TimeSunset):fLampLED    7.417    8.293  230.52  <2e-16 ***
#   s(Temp_mean_z)            6.154    7.173   52.26  <2e-16 ***
#   s(fNight)               230.360 1004.000  921.02  <2e-16 ***
#   s(fDevice)               21.334   23.000 1221.43  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.249   Deviance explained = 35.2%
# fREML = 4401.6  Scale est. = 1         n = 46272

# ***************************************************************

gam.check(mod3)
# The gam.check:
# k indexes are not small and significant, so all is good.
# The fNight k index is interesting, this indicates that it is important!

# Model rank =  1078 / 1078 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                              k'     edf k-index p-value
# s(TimeSunset):fLampdark    9.00    7.02    0.93    0.26
# s(TimeSunset):fLampHPS     9.00    6.92    0.93    0.22
# s(TimeSunset):fLampLED     9.00    7.42    0.93    0.26
# s(Temp_mean_z)             9.00    6.15    0.95    0.83
# s(fNight)               1008.00  230.36      NA      NA  # seems like fNight is very important. Keep in model!
# s(fDevice)                26.00   21.33      NA      NA

# ***************************************************************

concurvity(mod3, full = TRUE)
# Concurvity: ('multicollinearity')
#' Temperature and time-since-sunset are partially correlated ecological covariates, 
#' as seen also from this. However, their correlation was -0.3, so not major.
#' Also adding it to the model improved the fit (checked separately), so we keep it.

# para s(TimeSunset):fLampdark
# worst       1              0.64845199
# observed    1              0.08636799
# estimate    1              0.53613200
#                s(TimeSunset):fLampHPS s(TimeSunset):fLampLED
# worst                 0.4333296             0.49206875
# observed              0.1217705             0.03579845
# estimate              0.3875898             0.41266118
#             s(Temp_mean_z)  s(fNight) s(fDevice)
# worst         0.9501913 1.00000000          1
# observed      0.6796295 0.01981635          1
# estimate      0.8867479 0.03548521          1












# PREDICTIONS ----
sm <- smooth_estimates(mod3)
draw(sm)


## PREDICTIONS FOR LAMP ----
# Average effects (across the night)
newdat2 <- expand.grid(
  fLamp = levels(Dat_all$fLamp),
  fMonth = levels(Dat_all$fMonth),
  fYear = levels(Dat_all$fYear)[1],   # reference year
  fCity = levels(Dat_all$fCity)[1],   # reference city
  Temp_mean_z = mean(Dat_all$Temp_mean_z, na.rm = TRUE), # average temperature
  TimeSunset = mean(Dat_all$TimeSunset), # average time
  fDevice = levels(Dat_all$fDevice)[1],  # reference device
  fNight = levels(Dat_all$fNight)[1],    # reference night
  AR.start = TRUE                        # required for AR(1)
)

# All factor effects 
newdat2$fLamp  <- factor(newdat2$fLamp, levels = levels(Dat_all$fLamp))
newdat2$fMonth <- factor(newdat2$fMonth, levels = levels(Dat_all$fMonth))
newdat2$fYear  <- factor(newdat2$fYear, levels = levels(Dat_all$fYear))
newdat2$fCity  <- factor(newdat2$fCity, levels = levels(Dat_all$fCity))

# predict accross new data
pred2 <- predict(
  bat_mod_AR1,
  newdata = newdat2,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)") # population-level predictions
)

# Back-transform to probabilities
newdat2$fit <- plogis(pred2$fit)
newdat2$lwr <- plogis(pred2$fit - 1.96 * pred2$se.fit)
newdat2$upr <- plogis(pred2$fit + 1.96 * pred2$se.fit)
newdat2$fMonth <- factor(newdat2$fMonth, levels = c("Jun", "Aug"))


# plot with horizontal 

plotdat <- newdat2

plotdat$ypos <- c(
  ifelse(plotdat$fMonth == "Jun",
         as.numeric(plotdat$fLamp) - 0.10,
         as.numeric(plotdat$fLamp) + 0.10)
)

ggplot(
  plotdat,
  aes(
    x = fit,
    y = ypos,
    colour = fMonth
  )
) +
  geom_segment(
    aes(
      x = lwr,
      xend = upr,
      yend = ypos
    ),
    linewidth = 5,
    alpha = 0.35
  ) +
  geom_point(
    shape = 18,
    size = 4
  ) +
  scale_y_continuous(
    breaks = 1:3,
    labels = c("Unlit", "HPS", "LED")
  ) +
  scale_colour_manual(
    values = c(
      "Jun" = "#A879FF",
      "Aug" = "#013E89"
    )
  ) +
  labs(
    x = "Predicted probability of bat presence",
    y = "Treatment",
    colour = "Month"
  ) +
  theme_classic(base_size = 18)


# #plot
# # without 100%
ggplot(newdat2, aes(x = fLamp, y = fit, colour = fMonth, group = fMonth)) +
  geom_point(size = 5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0., position = position_dodge(width = 0.3)) +
  #geom_line(size = 1, position = position_dodge(width = 0.3)) +
  ylab("Predicted probability of bat presence") +
  xlab("Treatment") +
  scale_x_discrete(
    labels = c(
      "dark" = "Unlit",
      "HPS" = "HPS",
      "LED" = "LED"
    )
  )+
  scale_colour_manual(values = c("Jun" = "#B892FF", "Aug" = "#355070")) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top") +
  labs(colour = "Month")



# with 100%
# ggplot(newdat2, aes(x = fLamp, y = fit, colour = fMonth, group = fMonth)) +
#   geom_point(size = 3, position = position_dodge(width = 0.3)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr),
#                 width = 0., position = position_dodge(width = 0.3), alpha = 0.8) +
#   ylab("Predicted probability of bat presence") +
#   xlab("Lamp type") +
#   scale_colour_manual(values = c("Jun" = "#355070", "Aug" = "#503294")) +  # blue/red
#   scale_linetype_manual(values = c("Jun" = "solid", "Aug" = "dashed")) +
#   theme_classic(base_size = 18) +
#   theme(legend.position = "top",
#         legend.title = element_blank(),
#         axis.text.x = element_text(angle = 0, hjust = 0.5)) +
#   ylim(0, 1)





## PREDICTIONS FOR TIME SINCE SUNSET ----
newdat <- expand.grid(
  fLamp = levels(Dat_all$fLamp),
  fMonth = levels(Dat_all$fMonth),
  fYear = levels(Dat_all$fYear)[1],
  fCity = levels(Dat_all$fCity)[1],
  TimeSunset = seq(min(Dat_all$TimeSunset), max(Dat_all$TimeSunset), length.out = 100),
  fDevice = levels(Dat_all$fDevice)[1],  # <- pick a valid device
  Temp_mean_z = mean(Dat_all$Temp_mean_z, na.rm = TRUE),
  fNight  = levels(Dat_all$fNight)[1]
  
)
newdat$AR.start <- TRUE

# ensure factors match
newdat$fLamp  <- factor(newdat$fLamp, levels = levels(Dat_all$fLamp))
newdat$fMonth <- factor(newdat$fMonth, levels = levels(Dat_all$fMonth))
newdat$fYear  <- factor(newdat$fYear, levels = levels(Dat_all$fYear))
newdat$fCity  <- factor(newdat$fCity, levels = levels(Dat_all$fCity))

pred <- predict(
  bat_mod_AR1,
  newdata = newdat,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)")
)

newdat$fit <- plogis(pred$fit)
newdat$lwr <- plogis(pred$fit - 1.96*pred$se.fit)
newdat$upr <- plogis(pred$fit + 1.96*pred$se.fit)
newdat$fMonth <- factor(newdat$fMonth, levels = c("Jun", "Aug"))


### WITH HOURS

ggplot(newdat,
       aes(TimeSunset, fit,
           colour = fLamp,
           fill = fLamp)) +
  geom_line(linewidth = 1.4) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.2,
              colour = NA) +
  facet_wrap(~ fMonth) +
  scale_colour_manual(values = c(
    "dark" = "#00144B",
    "HPS"  = "#ff7f36",
    "LED"  = "#ffe736"
  ),   labels = c(
    "dark" = "Unlit",
    "HPS" = "HPS",
    "LED" = "LED"
  )) +
  scale_fill_manual(values = c(
    "dark" = "#00144B",
    "HPS"  = "#ff7f36",
    "LED"  = "#ffe736"
  ),   labels = c(
    "dark" = "Unlit",
    "HPS" = "HPS",
    "LED" = "LED"
  )
  ) +
  
  scale_x_continuous(
    breaks = c(1, 13, 25, 37, 49),
    labels = c("0", "2", "4", "6", "8")
  ) +
  labs(
    x = "Hours since sunset",
    y = "Predicted probability of bat presence",
    colour = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 17)


# Colours used in supplementary sensitivity analyses:
# "dark" = "#00144B",
# "HPS"  = "#E1783C",
# "LED"  = "#F9BC4C"

## PREDICTIONS FOR PEAK ACTIVITY ----

newdat %>%
  group_by(fMonth) %>%
  summarise(
    PeakTime = TimeSunset[which.max(fit)]
  )

june_peak_temp <- mean(
  Dat_all$Temp_mean_z[
    Dat_all$fMonth == "Jun" &
      Dat_all$TimeSunset == 16
  ],
  na.rm = TRUE
)



peakdat <- expand.grid(
  fLamp = levels(Dat_all$fLamp),
  fMonth = c("Jun", "Aug"),
  TimeSunset = 16,
  Temp_mean_z = june_peak_temp,
  fYear = levels(Dat_all$fYear)[1],
  fCity = levels(Dat_all$fCity)[1],
  fDevice = levels(Dat_all$fDevice)[1],
  fNight = levels(Dat_all$fNight)[1]
)

peakpred <- predict(
  mod3,
  newdata = peakdat,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)")
)

peakdat$fit <- plogis(peakpred$fit)

peakdat$lwr <- plogis(
  peakpred$fit - 1.96 * peakpred$se.fit
)

peakdat$upr <- plogis(
  peakpred$fit + 1.96 * peakpred$se.fit
)



## POST HOC ----
# Odds-ratio contrasts (this table) → “Which treatments differ on the odds scale?”
# compare to Model summary → “Which parameters shape the logit surface?” AND 
# Response-scale emmeans → “How big and certain are absolute ecological differences?”

emm <- emmeans(
  bat_mod_AR1,
  ~ fLamp | fMonth,
  type = "response"
)
emm
# fMonth = Aug:
#   fLamp    prob      SE    df lower.CL upper.CL
# dark  0.00596 0.00274 45985  0.00241   0.0146
# HPS   0.01668 0.00770 45985  0.00671   0.0408
# LED   0.01383 0.00626 45985  0.00567   0.0333
# 
# fMonth = Jun:
#   fLamp    prob      SE    df lower.CL upper.CL
# dark  0.01900 0.00845 45985  0.00790   0.0450
# HPS   0.03868 0.01730 45985  0.01592   0.0910
# LED   0.09988 0.03990 45985  0.04443   0.2094
# 
#' COMMENT: The emmans is averaged over night. 
#' Looking at the figure we can see tht largest differences are around main activity peak.
#' And that the treatment effect is not constant through the night.


contr <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = "dark",
  by = "fMonth",
  adjust = "none"  
)

summary(contr, infer = TRUE)
# fMonth = Aug:
#   contrast   odds.ratio   SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       2.83 1.82 45985    0.803     9.98    1   1.618  0.1056
# LED / dark       2.34 1.49 45985    0.670     8.17    1   1.332  0.1829
#' COMMENT: In August, bat presence tended to be higher at illuminated sites than
#'  at dark sites, but there was no statistically significant evidence that either 
#'  HPS or LED sites differed from dark sites.

# fMonth = Jun:
#   contrast   odds.ratio   SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       2.08 1.32 45985    0.601     7.19    1   1.155  0.2482
# LED / dark       5.73 3.56 45985    1.697    19.33    1   2.812  0.0049
#' COMMENT: In June, the probability of bat presence was substantially higher 
#' at LED sites than at dark sites. The odds of detecting bats at LED sites 
#' were estimated to be about 5.7 times higher than at dark sites (95% CI: 1.7-19.3).
#' COMMENT: HPS sites also showed higher estimated occurrence than dark sites (OR ≈ 2.1),
#'  but the uncertainty was large and there was no significant evidence for a difference.



# Thinking of the result:
# The model should be fine since the result didn't change after:
# - added temperature,
# - added a night random effect,
# - added a device random effect,
# - accounted for temporal autocorrelation,
# - allowed different activity curves through the night.
# 
# One concern could be that LED sites might systematically differ from dark sites 
# in habitat characteristics you did not measure.
# For example:
# water nearby
# forest edge proximity
# insect abundance
# commuting routes
# distance to roosts

#' However, the lamp effect doesn't seem to be purely a statistical artefact because:
#'   - Pure city effect seems unlikely since Espoo has all three lamptypes and August
#' shows no difference between lamp types
#' - Sites with led and hps could be different (see above) 
#' BUT this should still be reflected in August


# See the differences between how much bats were observed in June versus August:

aggregate(Ept ~ fMonth, Dat_all, mean)
# fMonth        Ept
# 1    Aug 0.01553732
# 2    Jun 0.06342138 # 3.9xmore
aggregate(Ept ~ fMonth + fLamp, Dat_all, mean)
# fMonth fLamp        Ept
# 1    Aug  dark 0.01143378
# 2    Jun  dark 0.03920966
# 3    Aug   HPS 0.03137860
# 4    Jun   HPS 0.08458458
# 5    Aug   LED 0.01151472
# 6    Jun   LED 0.08974359

#' More likely explanation for the lamp effect:
#' There is much more activity in general in June than in August (almost 4x more) (seen also from model output).
#' therefore, it is possible that the effect is visibly only in June since there 
#' is much more bat activity. The same might happen in August if we had enough bats
#' to see the pattern.Another explanation is that for some reason there is more bat 
#' activity at (led) lights in June compared to non-lit areas, and compared to August.
#' Why? Not sure. Ambient lighting? Insect community/composition? Something else?









# SENSITIVITY ANALYSES -----

### 0. PLOT OF THE DATA ----



# PLOT1 
#"Bat occurrence by device and year",

pres_prop <- Dat_all %>%
  group_by(fDevice, Year) %>%
  summarise(
    n = n(),
    presences = sum(Ept, na.rm = TRUE),
    prop = mean(Ept, na.rm = TRUE),
    .groups = "drop"
  )

# plot
ggplot(
  pres_prop,
  aes(x = fDevice,
      y = prop,
      fill = factor(Year))
) +
  geom_col(
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#BF4342",
      "2023" = "#978FA6" # #9389A6 # ADA8B6
    )
  ) +
  scale_y_continuous(
    limits = c(0, 0.155),
    expand = c(0, 0)
  )+
  labs(
    title = "Bat occurrence by site and year",
    x = "Site",
    y = "Proportion of intervals with bat presence",
    fill = "Year"
  )


# PLOT 2.1
# "Bat occurrence per device, lamp type (but not month)
# Is a small nro of devices driving the result?
# pres_prop <- Dat_all %>%
#   group_by(fDevice, fLamp) %>%
#   summarise(
#     prop = mean(Ept),
#     .groups = "drop"
#   )
# 
# ggplot(
#   pres_prop,
#   aes(reorder(fDevice, prop),
#       prop,
#       fill = fLamp)
# ) +
#   geom_col() +
#   coord_flip() +
#   theme_classic() +
#   xlab("Device") +
#   ylab("Bat occurrence proportion")
# 



# Plot 2.2
# "Bat occurrence per device, lamp type and month",
# same as above but per month

pres_prop_month <- Dat_all %>%
  group_by(fDevice, fMonth, fLamp) %>%
  summarise(
    prop = mean(Ept),
    .groups = "drop"
  )

pres_prop_month$fMonth <- factor(
  pres_prop_month$fMonth,
  levels = c("Jun", "Aug")
)

ggplot(
  pres_prop_month,
  aes(
    x = prop,
    y = reorder(fDevice, prop),
    fill = fLamp
  )
) +
  scale_fill_manual(
    values = c(
      "dark" = "#00144B",
      "HPS"  = "#E1783C",
      "LED"  = "#F9BC4C"
    )
  )+
  geom_col() +
  facet_wrap(~ fMonth, scales = "free_y") +
  xlab("Bat occurrence proportion") +
  ylab("Site") +
  theme_classic()+
  theme(
    strip.text = element_text(size = 16, face = "bold"), # facet titles
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    plot.title = element_text(size = 19,face = "bold",hjust = 0.5
    )
  )+
  labs(
    title = "Bat occurrence per site and month",
    x = "Bat occurrence proportion",
    y = "Site",
    fill = "Treatment"
  )



# plot 4

pres_prop <- Dat_all %>%
  group_by(fDevice, Year, fLamp) %>%
  summarise(
    prop = mean(Ept),
    .groups = "drop"
  )



ggplot(
  pres_prop,
  aes(x = factor(Year), y = prop, group = fDevice, colour = fLamp)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ fDevice) +
  scale_colour_manual(
    values = c(
      "dark" = "#00144B",
      "HPS"  = "#E1783C",
      "LED"  = "#F9BC4C"
    )
  ) +
  scale_x_discrete(
    limits = c("2022", "2023")
  ) +
  labs(
    title = "Bat occurrence proportion per site and year",
    x = "Year",
    y = "Proportion of intervals with bat presence",
    colour = "Treatment"
  ) +
  theme_gray() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 19,face = "bold",hjust = 0.5)
  )


### 1. TWO CHANGED DEVICES----

#' Devices/sites "E1L"& "E2L" changed from HPS to LED between 2022 and 2023
#' Here we fit a model with only these two devices/sites to see if they support the full model patterns
#' Note that the sample size is very small
table(Dat_all$fDevice, Dat_all$fLamp)
changed_devices <- c("E1L", "E2L")  

#' NOTE that we cannot separate year effect from lamp effect in this analyses,
#' Although in the full model fYear was not significant, the year effect does exist,
#'  but it was not huge. fYear2023 = 0.27, p = 0.019
xtabs(~ fLamp + fYear, Dat_change) # year and lamp confounding fully


# New dat_all with only changed devices
Dat_change <- subset(
  Dat_all,
  fDevice %in% changed_devices
)
Dat_change <- Dat_change[order(Dat_change$fNight,
                               Dat_change$TimeSunset), ]
Dat_change <- droplevels(Dat_change)
with(Dat_change, table(fDevice, Year, fLamp))

# New AR.start for the subset model
Dat_change$AR.start <- c(
  TRUE,
  Dat_change$fNight[-1] != Dat_change$fNight[-nrow(Dat_change)]
)
AR.start = Dat_change$AR.start

# Re-estimate rho for the subset model
m0.change <- bam(
  Ept ~
    fLamp * fMonth +
    s(TimeSunset, by = fLamp, k = 10) +
    s(Temp_mean_z) +
    s(fNight, bs = "re"),
  family = binomial,
  data = Dat_change,
  method = "fREML",
  discrete = TRUE
)

rsd <- residuals(m0.change, type = "deviance")
mf <- model.frame(m0.change)
rho_vals <- tapply(
  rsd,
  mf$fNight,
  function(x) {
    if(length(x) > 2)
      cor(x[-1], x[-length(x)])
    else
      NA
  }
)
rho_est_change <- mean(rho_vals, na.rm = TRUE)

#  NEW SUBSET MODEL 
mod3.change <-bam(
  Ept ~
    fLamp +
    fDevice + 
    s(TimeSunset, k = 10) +
    s(Temp_mean_z) +
    s(fNight, bs = "re"),
  family = binomial,
  data = Dat_change,
  method = "fREML",
  rho = rho_est_change,
  AR.start = Dat_change$AR.start,
  discrete = TRUE
)


# RAW DATA
aggregate(Ept ~ fDevice + fLamp, Dat_change, mean)
# fDevice fLamp        Ept
# 1     E1L   HPS 0.06593407
# 2     E2L   HPS 0.03479853
# 3     E1L   LED 0.10073260
# 4     E2L   LED 0.10989011
#' Both devices show higher bat occurrence during the LED year.
#' If they didn't move in the same direction - could be noise. 
#' But here they do move both up (~50% and 200%)


# SUBSET MODEL
summary(mod3.change)
# Ept ~ fLamp + fDevice + s(TimeSunset, k = 10) + s(Temp_mean_z) + 
#   s(fNight, bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.4612     0.2255 -15.348   <2e-16 ***
#   fLampLED      0.5649     0.2333   2.421   0.0155 *  
#   fDeviceE2L   -0.3518     0.2195  -1.603   0.1090    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(TimeSunset)   5.923  7.056 110.943  <2e-16 ***
#   s(Temp_mean_z)  3.314  4.072   8.945  0.0656 .  
# s(fNight)      31.677 89.000  54.965   5e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.178   Deviance explained = 22.7%
# fREML = 725.17  Scale est. = 1         n = 4220
# 
#' For these two devices, the odds of detecting Eptesicus
#' were estimated to be about 1.76 (exp(0.5649)) times greater during the
#' LED year than during the HPS year.
#' Note that there might be still a Year effect not seen in this subset model.
#' However, the same direction of difference was found in the full model after accounting 
#' for year, and year effect 2022<2023 was small

# CONFIDENCE INTERVAL FOR THE MODEL ABOVE
tab <- summary(mod3.change)$p.table
out <- data.frame(
  Estimate = tab[, "Estimate"],
  SE       = tab[, "Std. Error"],
  OR       = exp(tab[, "Estimate"]),
  Lower95  = exp(tab[, "Estimate"] - 1.96*tab[, "Std. Error"]),
  Upper95  = exp(tab[, "Estimate"] + 1.96*tab[, "Std. Error"]),
  p        = tab[, "Pr(>|z|)"]
)
round(out, 4)
# > round(out, 4)
#             Estimate     SE     OR Lower95 Upper95      p
# (Intercept)  -3.4612 0.2255 0.0314  0.0202  0.0488 0.0000
# fLampLED      0.5649 0.2333 1.7593  1.1136  2.7794 0.0155
# fDeviceE2L   -0.3518 0.2195 0.7034  0.4575  1.0816 0.1090
#
# LED EFFECT: Odds ratio = 1.76 (95% CI 1.11-2.78, p = 0.016). This is for both lamps together.

###  2 YEARLY MODELS - years separately ----
# +   bat_mod_AR1.22
# > summary(contr, infer = TRUE)
# fMonth = Aug:
#   contrast   odds.ratio    SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       4.32 4.250 21693    0.630    29.66    1   1.490  0.1363
# LED / dark       1.09 1.140 21693    0.138     8.56    1   0.079  0.9373
# 
# fMonth = Jun:
#   contrast   odds.ratio    SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       1.28 1.240 21693    0.191     8.53    1   0.252  0.8010
# LED / dark       0.94 0.957 21693    0.128     6.91    1  -0.061  0.9514
# 
# +  bat_mod_AR1.23
# 
# > summary(contr, infer = TRUE)
# fMonth = Aug:
#   contrast   odds.ratio    SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       2.74  2.42 24278    0.488     15.4    1   1.145  0.2522
# LED / dark       3.37  2.46 24278    0.806     14.1    1   1.665  0.0959
# fMonth = Jun:
#   contrast   odds.ratio    SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       6.58  5.69 24278    1.207     35.9    1   2.177  0.0295
# LED / dark      17.76 12.60 24278    4.408     71.6    1   4.046  0.0001
#
#* COMMENT: What we found is that the lighting effect is much stronger in 2023 than in 2022.
#* In 2023 in June there were stronger effects of both lamp types having an increasing effect
#*  on prob. of occurrence for enil, the effect was stronger for LED but also there for HPS. 
#*  In August again no effect. In 2022 there are no significant effects related to lighting types.
#*  It is also notable, that when we subset the data per year we are also halving the sample size,
#*  and confidence intervals become very big. Therefore we can't say that the effect only exists in 2023.
#*  For example, this interval 0.13–6.91 (June LED / dark 2022) is enormous, and contains both negative,
#*  positive an neutral effects.  This means we do not know what is happening in 2022, 
#*  and that is why the full model results are mainly influenced by 2023.
#*  

Dat_all22<-Dat_all[which(Dat_all$Year=="2022"),]
Dat_all22 <- droplevels(Dat_all22)
Dat_all23<-Dat_all[which(Dat_all$Year=="2023"),]
Dat_all23 <- droplevels(Dat_all23)

#### YEAR 2022 ----
bat_mod22 <- bam( Ept ~ fLamp * fMonth  + fCity + s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_all22, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res22 <- residuals(bat_mod22, type = "pearson")
rho_est22 <- acf(res22, plot = FALSE)$acf[2]
rho_est22


bat_mod_AR1.22 <- bam(
  Ept ~ 
    fLamp * fMonth  + fCity +
    s(TimeSunset, by = fLamp, k = 10) + # treatment specific time-relationship
    s(Temp_mean_z, k=15)+ # temperature, collinearity with time -0.31
    s(fNight, bs="re") + # per night variability
    s(fDevice, bs = "re"),
  family = binomial,
  data = Dat_all22,
  method = "fREML",
  rho = rho_est22,
  AR.start = Dat_all22$AR.start,
  discrete = TRUE
)

summary(bat_mod_AR1.22)
gam.check(bat_mod_AR1.22) # this is ok, k indexes are not low or significant

emm <- emmeans(
  bat_mod_AR1.22,
  ~ fLamp | fMonth,
  type = "response"
)
emm
cld(emm, adjust = "tukey", Letters = letters) 


contr <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = "dark",
  by = "fMonth",
  adjust = "none"  
)

summary(contr, infer = TRUE)


#### YEAR 2023 ----
bat_mod23 <- bam( Ept ~ fLamp * fMonth  + fCity + s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_all23, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res23 <- residuals(bat_mod23, type = "pearson")
rho_est23 <- acf(res23, plot = FALSE)$acf[2]
rho_est23


bat_mod_AR1.23 <- bam(
  Ept ~ 
    fLamp * fMonth  + fCity +
    s(TimeSunset, by = fLamp, k = 10) + # treatment specific time-relationship
    s(Temp_mean_z)+ # temperature, collinearity with time -0.31
    s(fNight, bs="re") + # per night variability
    s(fDevice, bs = "re"),
  family = binomial,
  data = Dat_all23,
  method = "fREML",
  rho = rho_est23,
  AR.start = Dat_all23$AR.start,
  discrete = TRUE
)

summary(bat_mod_AR1.23)
gam.check(bat_mod_AR1.23) # this is ok, k indexes are not low or significant

emm <- emmeans(
  bat_mod_AR1.23,
  ~ fLamp | fMonth,
  type = "response"
)
emm
cld(emm, adjust = "tukey", Letters = letters) 


contr <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = "dark",
  by = "fMonth",
  adjust = "none"  
)

summary(contr, infer = TRUE)



### 3.  MODEL WITHOUT E5L ----

#' E5L has a lot of enil activity
#' What we found was that E5L conributes to the effect,
#' but E5L is not solely responsible for the effect.
#' Without E5L the effect is smaller but still there.
#' OR 5.73 → 3.85, but p is still significant, and 95% CI = 1.09-13.67 excludes 1.
#'Ccontrast table is following:
# fMonth = Aug:
#   contrast   odds.ratio   SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       2.62 1.72 43897    0.721     9.51    1   1.463  0.1434
# LED / dark       2.35 1.56 43897    0.639     8.62    1   1.285  0.1990
# 
# fMonth = Jun:
#   contrast   odds.ratio   SE    df lower.CL upper.CL null t.ratio p.value
# HPS / dark       1.93 1.25 43897    0.541     6.88    1   1.013  0.3109
# LED / dark       3.85 2.49 43897    1.085    13.67    1   2.086  0.0370

Dat_noE5L <- subset(
  Dat_all,
  fDevice != "E5L"
)

bat_mod2 <- bam( Ept ~ fLamp * fMonth + fYear + fCity + s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_noE5L, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res2 <- residuals(bat_mod2, type = "pearson")
rho_est2 <- acf(res2, plot = FALSE)$acf[2]
rho_est2


bat_mod_AR1.noE5L <- bam(
  Ept ~ 
    fLamp * fMonth + fYear + fCity +
    s(TimeSunset, by = fLamp, k = 10) + # treatment specific time-relationship
    s(Temp_mean_z, k = 15)+ # temperature, collinearity with time -0.31
    s(fNight, bs="re") + # per night variability
    s(fDevice, bs = "re"),
  family = binomial,
  data = Dat_noE5L,
  method = "fREML",
  rho = rho_est2,
  AR.start = Dat_noE5L$AR.start,
  discrete = TRUE
)

summary(bat_mod_AR1.noE5L)
gam.check(bat_mod_AR1.noE5L) 

emm <- emmeans(
  bat_mod_AR1.noE5L,
  ~ fLamp | fMonth,
  type = "response"
)
emm

contr <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = "dark",
  by = "fMonth",
  adjust = "none"  
)
summary(contr, infer = TRUE)





##### PREDICTIONS ----
newdat <- expand.grid(
  fLamp = levels(Dat_noE5L$fLamp),
  fMonth = levels(Dat_noE5L$fMonth),
  fYear = levels(Dat_noE5L$fYear)[1],
  fCity = levels(Dat_noE5L$fCity)[1],
  TimeSunset = seq(min(Dat_noE5L$TimeSunset), max(Dat_noE5L$TimeSunset), length.out = 100),
  fDevice = levels(Dat_noE5L$fDevice)[1],  # <- pick a valid device
  Temp_mean_z = mean(Dat_noE5L$Temp_mean_z, na.rm = TRUE),
  fNight  = levels(Dat_noE5L$fNight)[1]
  
)
newdat$AR.start <- TRUE

# ensure factors match
newdat$fLamp  <- factor(newdat$fLamp, levels = levels(Dat_noE5L$fLamp))
newdat$fMonth <- factor(newdat$fMonth, levels = levels(Dat_noE5L$fMonth))
newdat$fYear  <- factor(newdat$fYear, levels = levels(Dat_noE5L$fYear))
newdat$fCity  <- factor(newdat$fCity, levels = levels(Dat_noE5L$fCity))

pred <- predict(
  bat_mod_AR1.noE5L,
  newdata = newdat,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)")
)

newdat$fit <- plogis(pred$fit)
newdat$lwr <- plogis(pred$fit - 1.96*pred$se.fit)
newdat$upr <- plogis(pred$fit + 1.96*pred$se.fit)
newdat$fMonth <- factor(newdat$fMonth, levels = c("Jun", "Aug"))


### WITH HOURS

ggplot(newdat,
       aes(TimeSunset, fit,
           colour = fLamp,
           fill = fLamp)) +
  geom_line(linewidth = 1.4) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.2,
              colour = NA) +
  facet_wrap(~ fMonth) +
  scale_colour_manual(values = c(
    "dark" = "#00144B",
    "HPS"  = "#ff7f36",
    "LED"  = "#ffe736"
  ),   labels = c(
    "dark" = "Unlit",
    "HPS" = "HPS",
    "LED" = "LED"
  )) +
  scale_fill_manual(values = c(
    "dark" = "#00144B",
    "HPS"  = "#ff7f36",
    "LED"  = "#ffe736"
  ),   labels = c(
    "dark" = "Unlit",
    "HPS" = "HPS",
    "LED" = "LED"
  )
  ) +
  
  scale_x_continuous(
    breaks = c(1, 13, 25, 37, 49),
    labels = c("0", "2", "4", "6", "8")
  ) +
  labs(
    x = "Hours since sunset",
    y = "Predicted probability of bat presence",
    colour = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 17)


### 4. Check LUX values ----
Dat_all$logLux <- log1p(Dat_all$Lux)

# LUX
summary(Dat_all$Lux)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.00    3.19   25.33   21.30 2915.84    5848 

hist(Dat_all$Lux, breaks = 100)
hist(log1p(Dat_all$Lux), breaks = 100)

# Month and lux
## Lux differences per month in all sites and lamp types
aggregate(Lux ~ fMonth, Dat_all, mean, na.rm = TRUE)
aggregate(Lux ~ fMonth + fLamp + fCity, Dat_all, mean, na.rm = TRUE)

ggplot(Dat_all,aes(fMonth, logLux)) + geom_boxplot() 
ggplot(subset(Dat_all, fLamp == "dark"),aes(fMonth, Lux)) +
  geom_boxplot() + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))

# plot for appendix
Dat_all$fMonth2 <- factor(Dat_all$fMonth, levels = c("Jun", "Aug"))

ggplot(Dat_all,aes(fMonth2, logLux)) + geom_boxplot() +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 12)
  ) +
  labs(
    title = "Lux per month",
    x = "Month",
    y = " logLux"
  )

# Lux and TimeSunset
## Correlation
cor(Dat_all$Lux, #  -0.2475363
    Dat_all$TimeSunset,
    use = "complete.obs")
cor(Dat_all$logLux,#-0.2363147
    Dat_all$TimeSunset,
    use = "complete.obs")

# Lamp
## We can't add lux to the main model because the lux differs between lamp types so much.
aggregate(Lux ~ fLamp,Dat_all,mean,na.rm = TRUE)
aggregate(Lux ~ fMonth + fLamp,Dat_all,mean,na.rm = TRUE)



### 5. Lux analyses with jun only

Dat_all_jun<-Dat_all[which(Dat_all$Month=="Jun"),]
Dat_all_jun <- droplevels(Dat_all_jun)


##### PREDICTIONS FOR LAMP ----
# Average effects (across the night)

newdat2 <- expand.grid(
  fLamp = levels(Dat_all$fLamp),
  #  fMonth = levels(Dat_all$fMonth),
  fYear = levels(Dat_all$fYear)[1],   # reference year
  fCity = levels(Dat_all$fCity)[1],   # reference city
  TimeSunset = mean(Dat_all$TimeSunset), # average time
  fDevice = levels(Dat_all$fDevice)[1],  # reference device
  AR.start = TRUE                        # required for AR(1)
)

# Make sure factors match exactly
newdat2$fLamp  <- factor(newdat2$fLamp, levels = levels(Dat_all$fLamp))
#newdat2$fMonth <- factor(newdat2$fMonth, levels = levels(Dat_all$fMonth))
newdat2$fYear  <- factor(newdat2$fYear, levels = levels(Dat_all$fYear))
newdat2$fCity  <- factor(newdat2$fCity, levels = levels(Dat_all$fCity))

pred2 <- predict(
  bat_mod_AR1,
  newdata = newdat2,
  type = "link",
  se.fit = TRUE,
  exclude = "s(fDevice)"  # population-level predictions
)

# Back-transform to probabilities
newdat2$fit <- plogis(pred2$fit)
newdat2$lwr <- plogis(pred2$fit - 1.96 * pred2$se.fit)
newdat2$upr <- plogis(pred2$fit + 1.96 * pred2$se.fit)
newdat2$fMonth <- factor(newdat2$fMonth, levels = c("Jun", "Aug"))

# ggplot(newdat2, aes(x = fLamp, y = fit, colour = fLamp, group = fLamp)) +
#   geom_point(size = 3, position = position_dodge(width = 0.2)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, position = position_dodge(width = 0.2)) +
#   geom_line(aes(group = fLamp), position = position_dodge(width = 0.2)) +
#   facet_wrap(~ fMonth) +
#   ylab("Predicted probability of bat presence") +
#   xlab("Lamp type") +
#   theme_classic(base_size = 14) +
#   theme(legend.position = "none")

# without 100%
ggplot(newdat2, aes(x = fLamp, y = fit)) + # colour = fMonth, group = fMonth
  geom_point(size = 5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0., position = position_dodge(width = 0.3)) +
  #geom_line(size = 1, position = position_dodge(width = 0.3)) +
  ylab("Predicted probability of bat presence") +
  xlab("Lamp type") +
  #scale_colour_manual(values = c("Jun" = "#bf8c0a", "Aug" = "#503294")) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top") +
  labs(colour = "Month")




##### SENSITIVITY ANALYSES END -------------




# B) LUX MODEL EPT ----

# Exploring effect of lux to bat presence
ggplot(
  subset(Dat_all, fLamp == "dark" & !is.na(Lux)),
  aes(log1p(Lux), Ept)
) + geom_smooth(method = "gam")

# Subset
Dat_dark <- droplevels(
  subset(Dat_all, fLamp == "dark")
)


# collinearity
Dat_dark$logLux<-log1p(Dat_dark$Lux)
cor(Dat_dark$Temp_mean_z, Dat_dark$logLux, use = "complete.obs") # -0.3085

bat_mod <- bam( Ept ~ fMonth + fYear + fCity +  s(log1p(Lux)) +s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_dark, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res <- residuals(bat_mod, type = "pearson")
rho_est <- acf(res, plot = FALSE)$acf[2]
rho_est


# Model 
mod_dark_lux <- bam(
  Ept ~
    fMonth +
    fYear +
    fCity +
    s(TimeSunset, k = 10) +
    s(log1p(Lux)) +
    s(Temp_mean_z) +
    s(fNight, bs = "re") +
    s(fDevice, bs = "re"),
  family = binomial,
  data = Dat_dark,
  method = "fREML",
  rho = rho_est,
  AR.start = Dat_dark$AR.start,
  discrete = TRUE
)

# Check model
summary(mod_dark_lux)
gam.check(mod_dark_lux) # this is quite ok, although k in lux almost significant

# Autocorrelation
acf(residuals(mod_dark_lux, type = "pearson")) ; acf(residuals(mod_dark_lux, type = "pearson"), plot = FALSE)$acf[2] # 0.06 - 0.11 is ok for the second bar!
res2 <- residuals(mod_dark_lux, type = "pearson")

# DHARMa
sim_res <- simulateResiduals(mod_dark_lux, n = 1000);plot(sim_res)
testZeroInflation(sim_res) #ok


# 'Multicollinearity' 
concurvity(mod_dark_lux)
#' It seems to be impossible to separate the effect of TimeSunset, Lux (and Temperature, but very little),
#' since they change together. Therefore we can't conclude wether the bats are more or
#' less likely to occur under specific lux conditions.
#' However, here we can conclude that lux levels in June are very different from august,
#' and that also at dark sites ept are more present in June than August.

# Partial effect
draw(mod_dark_lux)


## COMPARISON between smothers ----

mod_dark_noTemp <- bam(
  Ept ~
    fMonth +
    fYear +
    fCity +
    s(TimeSunset, k = 10) +
    s(log1p(Lux)) +
    s(fNight, bs = "re") +
    s(fDevice, bs = "re"),
  family = binomial,
  data = Dat_dark,
  method = "fREML",
  rho = rho_est,
  AR.start = Dat_dark$AR.start,
  discrete = TRUE
)


AIC(mod_dark_lux, mod_dark_noTemp)
concurvity(mod_dark_noTemp)
summary(mod_dark_noTemp)
draw(mod_dark_noTemp)
#' Removing temperature made almost no difference in the model
#' AIC difference very small. Most problems seem to arise from
#' Lux - TimeSunset, not temperature.


## PREDICTIONS FOR LUX ----

newdat_lux <- expand.grid(
  fMonth = levels(Dat_dark$fMonth),
  fYear = levels(Dat_dark$fYear)[1],
  fCity = levels(Dat_dark$fCity)[1],
  TimeSunset = mean(Dat_dark$TimeSunset),
  Temp_mean_z = mean(Dat_dark$Temp_mean_z, na.rm = TRUE),
  
  Lux = seq(
    quantile(Dat_dark$Lux, 0.01, na.rm = TRUE),
    quantile(Dat_dark$Lux, 0.99, na.rm = TRUE),
    length.out = 100
  )  ,
  fDevice = levels(Dat_dark$fDevice)[1],
  fNight  = levels(Dat_dark$fNight)[1]
)

pred <- predict(
  mod_dark_lux,
  newdata = newdat_lux,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)")
)


newdat_lux$fit <- plogis(pred$fit)
newdat_lux$lwr <- plogis(pred$fit - 1.96 * pred$se.fit)
newdat_lux$upr <- plogis(pred$fit + 1.96 * pred$se.fit)
newdat_lux$fMonth <- factor(newdat_lux$fMonth, levels = c("Jun", "Aug"))
newdat_lux$logLux <- log1p(newdat_lux$Lux)


# plot
ggplot(newdat_lux, aes(x = logLux, y = fit, colour = fMonth, fill = fMonth)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2, colour = NA
  ) +
  theme_classic(base_size = 17) +
  labs(
    x = "log(1 + Lux)",
    y = "Predicted probability of bat presence at unlit sites",
    colour = "Month",
    fill = "Month"
  )+  scale_colour_manual(values = c("Jun" = "#bf8c0a", "Aug" = "#503294")) 


## COMMENTS on LUX results ----
#' Results/discu paragarph based on lux model:
#' Supplementary analyses only on unlit sites revealed that light levels (lux) 
#' were associated negatively with ept bat occurrence. However, it was impossible to 
#' distinguish the effect of light levels and time of the night, which were highly correlated.
#' However, month remained an important predictor to bat occurrence, meaning that ambient light level
#' was not the only reason for the effect of month seen in the results. Other potential reasons for June 
#' (June odds ratio = exp(2.731) = 15.3)
#' to predict more bats could be seasonal movements, breeding phenology and localities or prey availability among others.

