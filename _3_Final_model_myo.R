
# Adapted from 10.7717/peerj.6876 
# Hierarchical generalized additive models in ecology: an introduction with mgcv, Pedersen et al 2019
# Tried first to run a GS model from the paper but this was not possible with my data (too long run), therefore i have this more simple gam


library(dplyr)
library(ggplot2)
library(DHARMa)
library(gratia)
library(emmeans)
library(mgcv)

# Note: dark = unlit

Dat_all$fNight<-factor(interaction(Dat_all$Device, Dat_all$Date_12, drop = TRUE))
#table(Dat_all$Night)[1:5]
#str(Dat_all)
#bat_test <- Dat_all[Dat_all$Date_12 %in% sample(unique(Dat_all$Date_12), 20), ]

Dat_all <- Dat_all[order(Dat_all$fNight, Dat_all$TimeSunset), ]
Dat_all$AR.start <- c(TRUE, Dat_all$fNight[-1] != Dat_all$fNight[-nrow(Dat_all)])

#run model for rho

bat_modM <- bam( Myo ~ fTreatment * fMonth + fYear + fCity + s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_all, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res <- residuals(bat_modM, type = "pearson")
rho_est <- acf(res, plot = FALSE)$acf[2]
rho_est

# GAM ----
set.seed(123)
bat_modM_AR1 <- bam(
  Myo ~ 
    fTreatment* fMonth + fYear + fCity + # no temperature because collinear with timeSunset
    s(TimeSunset, k = 10) + # in models 1 & 4
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


summary(bat_modM_AR1)
# AIC(bat_modM_AR11,bat_modM_AR1)
gam.check(bat_modM_AR1) # this is all ok
acf(residuals(bat_modM_AR1, type = "pearson")) ; acf(residuals(bat_modM_AR1, type = "pearson"), plot = FALSE)$acf[2] # 0.11 is ok for the second bar!
res2 <- residuals(bat_modM_AR1, type = "pearson")

# plot(Dat_all$TimeSunset, res2, pch = 16, cex = 0.3) #
# abline(h = 0, col = "red")

# residuals
sim_res <- simulateResiduals(bat_modM_AR1, n = 1000)
plot(sim_res)
#testDispersion(sim_res)
testZeroInflation(sim_res) #ok

# plot(
#   fitted(bat_modM_AR1),
#   Dat_all$Myo,
#   pch = 16,
#   cex = 0.3

# 'multicollinearity' 
concurvity(bat_modM_AR1)

# effect
draw(bat_modM_AR11)
smooth_estimates(bat_modM_AR1)

# PREDICTIONS FOR LAMP/Treatment ----
# Average effects (across the night)

newdat2 <- expand.grid(
  fTreatment = levels(Dat_all$fTreatment),
  fMonth = levels(Dat_all$fMonth),
  fYear = levels(Dat_all$fYear)[1],   # reference year
  fCity = levels(Dat_all$fCity)[1],   # reference city
  TimeSunset = mean(Dat_all$TimeSunset), # average time
  fDevice = levels(Dat_all$fDevice)[1],  # reference device
  fNight = levels(Dat_all$fNight)[1],    # reference night
  Temp_mean_z = mean(Dat_all$Temp_mean_z, na.rm = TRUE), # average temperature
  AR.start = TRUE                        # required for AR(1)
)

# Make sure factors match exactly
newdat2$fTreatment  <- factor(newdat2$fTreatment, levels = levels(Dat_all$fTreatment))
newdat2$fMonth <- factor(newdat2$fMonth, levels = levels(Dat_all$fMonth))
newdat2$fYear  <- factor(newdat2$fYear, levels = levels(Dat_all$fYear))
newdat2$fCity  <- factor(newdat2$fCity, levels = levels(Dat_all$fCity))

pred2 <- predict(
  bat_modM_AR1,
  newdata = newdat2,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)")  # population-level predictions
)

# Back-transform to probabilities
newdat2$fit <- plogis(pred2$fit)
newdat2$lwr <- plogis(pred2$fit - 1.96 * pred2$se.fit)
newdat2$upr <- plogis(pred2$fit + 1.96 * pred2$se.fit)
newdat2$fMonth <- factor(newdat2$fMonth, levels = c("Jun", "Aug"))



newdat_treat <- newdat2 %>%
  group_by(fTreatment) %>%
  summarise(
    fit = mean(fit),
    lwr = mean(lwr),
    upr = mean(upr)
  )


# Black plot
ggplot(newdat_treat, aes(x = fTreatment, y = fit)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.) +
  ylab("Predicted probability of bat presence") +
  xlab("Treatment") +
  #ylim(0, 1) +
  theme_classic(base_size = 18)


ggplot(
  newdat_treat,
  aes(x = fTreatment, y = fit,colour = fTreatment)
) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),width = 0) +
  scale_x_discrete(
    labels = c(
      "Dark" = "Unlit",
      "Light"  = "Lit"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Dark" = "#00144B",  # dark blue
      "Light"  = "#F7C42D"   
    ),
    labels = c(
      "Dark" = "Unlit",
      "Light"  = "Lit"
    )
  ) +
  ylab("Predicted probability of bat presence") +
  xlab("Treatment") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none")


# HORIZONTAL style
ggplot(
  newdat_treat,
  aes(
    y = fTreatment,
    x = fit
  )
) +
  geom_segment(
    aes(
      x = lwr,
      xend = upr,
      yend = fTreatment
    ),
    linewidth = 6,
    alpha = 0.35,
    colour = c("black", "yellow")
  ) +
  geom_point(
    aes(colour = fTreatment),
    shape = 18,
    size = 4
  ) +
  scale_colour_manual(
    values = c(
      "Dark"  = "black",
      "Light" = "orange"
    ),
    labels = c(
      "Dark"  = "Unlit",
      "Light" = "Lit"
    )
  ) +
  scale_y_discrete(
    labels = c(
      "Dark"  = "Unlit",
      "Light" = "Lit"
    )
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(
    x = "Predicted probability of bat presence",
    y = "Treatment"
  )




# PREDICTIONS FOR TIME SINCE SUNSET ----

newdat_time <- expand.grid(
  fTreatment = levels(Dat_all$fTreatment),
  fMonth     = levels(Dat_all$fMonth),   # add [1] to the end if you want just one graph for both months
  fYear      = levels(Dat_all$fYear)[1],
  fCity      = levels(Dat_all$fCity)[1],
  TimeSunset = seq(min(Dat_all$TimeSunset),
                   max(Dat_all$TimeSunset),
                   length.out = 100),
  Temp_mean_z = mean(Dat_all$Temp_mean_z, na.rm = TRUE),
  fNight  = levels(Dat_all$fNight)[1],
  fDevice    = levels(Dat_all$fDevice)[1],
  AR.start   = TRUE
)

# Ensure factor levels match model
newdat_time$fTreatment <- factor(newdat_time$fTreatment, levels = levels(Dat_all$fTreatment))
newdat_time$fMonth     <- factor(newdat_time$fMonth, levels = levels(Dat_all$fMonth))
newdat_time$fYear  <- factor(newdat_time$fYear, levels = levels(Dat_all$fYear))
newdat_time$fCity  <- factor(newdat_time$fCity, levels = levels(Dat_all$fCity))


pred_time <- predict(
  bat_modM_AR1,
  newdata = newdat_time,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(fDevice)", "s(fNight)")
)

newdat_time$fit <- plogis(pred_time$fit)
newdat_time$lwr <- plogis(pred_time$fit - 1.96 * pred_time$se.fit)
newdat_time$upr <- plogis(pred_time$fit + 1.96 * pred_time$se.fit)


# ggplot(newdat_time,
#        aes(x = TimeSunset, y = fit, colour = fTreatment)) +
#   geom_line(size = 1.2) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr, fill = fTreatment),
#               alpha = 0.15, colour = NA) +
#   ylab("Predicted probability of bat presence") +
#   xlab("Time since sunset (minutes)") +
# #  ylim(0, 1) +
#   theme_classic(base_size = 18) +
#   theme(legend.position = "top")
# 

ggplot(newdat_time,
       aes(x = TimeSunset, y = fit, colour = fTreatment, fill = fTreatment)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.15, colour = NA) +
  scale_colour_manual(values = c(
    "Dark"  = "black",
    "Light" = "#F7C42D"   # orange
  ),   labels = c(
    "Dark" = "Unlit",
    "Light" = "Lit"
  )) +
  scale_fill_manual(values = c(
    "Dark"  = "black",
    "Light" = "#F7C42D"
  ),   labels = c(
    "Dark" = "Unlit",
    "Light" = "Lit"
  ))+
  
  scale_x_continuous(
    breaks = c(1, 13, 25, 37, 49),
    labels = c("0", "2", "4", "6", "8")
  )  +
  labs(
    x = "Hours since sunset",
    y = "Predicted probability of bat presence",
    colour = "Treatment",
    fill = "Treatment"
  ) +
 # ylim(0, 1) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top")



# separate months





newdat_time$fMonth <- factor(
  newdat_time$fMonth,
  levels = c("Jun", "Aug")
)

ggplot(
  newdat_time,
  aes(
    x = TimeSunset,
    y = fit,
    colour = fTreatment,
    fill = fTreatment
  )
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.15,
    colour = NA
  ) +
  facet_wrap(~ fMonth, nrow = 1) +
  scale_colour_manual(
    values = c(
      "Dark"  = "#00144B",
      "Light" = "#F7C42D"
    ),
    labels = c(
      "Dark"  = "Unlit",
      "Light" = "Lit"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Dark"  = "#00144B", #     "dark" = "#00144B",
      "Light" = "#F7C42D"
    ),
    labels = c(
      "Dark"  = "Unlit",
      "Light" = "Lit"
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
  theme_classic(base_size = 18) +
  theme(
   # legend.position = "top",
    strip.text = element_text(size = 16, face = "bold")
  )






# COMPARISON BETWEEN SMOOTHERS----
#' Here we check four combination of different smoothers / variables to add to the model
#' Including night-to-night variation improved model by a lot,
#' Including separate TimeSunset smooth per treatment improved model a bit (but enough)
#' Therefore, both of these additions are included in the final model (mod3 chosen as the best model based on AIC)

# ***************************************************************
# MODELS:

mod4 <- bam(
  Myo ~ 
    fTreatment* fMonth + fYear + fCity + # no temperature because collinear with timeSunset
    #  s(TimeSunset, by = fTreatment, k = 10) + # treatment specific time-relationship
    s(TimeSunset, k = 10) + # in models 1 & 4
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



# # > mod1
# Myo ~ fTreatment * fMonth + fYear + fCity + s(TimeSunset, k = 10) + 
#   s(Temp_mean_z) + s(fDevice, bs = "re")
# # > mod2
# Myo ~ fTreatment * fMonth + fYear + fCity + s(TimeSunset, by = fTreatment, 
#   k = 10) + s(Temp_mean_z) + s(fDevice, bs = "re")
# # > mod3
# Myo ~ fTreatment * fMonth + fYear + fCity + s(TimeSunset, by = fTreatment, 
#   k = 10) + s(Temp_mean_z) + s(fNight, bs = "re") + s(fDevice, 
# # > mod4                                                                                                  bs = "re")# > mod4
# Myo ~ fTreatment * fMonth + fYear + fCity + s(TimeSunset, k = 10) + 
#   s(Temp_mean_z) + s(fNight, bs = "re") + s(fDevice, bs = "re")
# ***************************************************************
# Different smoothers
AIC(mod1, mod2,mod3, mod4)
#             df      AIC
# mod1  39.30693 2969.433
# mod2  41.91986 2971.788
# mod3 147.61858 2877.449
# mod4 144.90866 2875.096
# 
# Difference between mod3/4 and mod1/2: 
#   mod3 and mod4 compared to m1/m2: adding night-to-night variation improved model a lot
# Difference between mod 3 and mod 4: 
#   mod3: separate TimeSunset smooths for each lamp, mod4: one TimeSunset smooth. AIC difference 2.353
# Not much difference, but mod4 has lower AIC and is more simple.
# --> Going forward with model mod4

# ***************************************************************
# Temperature
AIC(mod4, mod6)  
#            df      AIC
# mod4 144.90866 2875.096
# mod6  31.53019 3159.940 # ΔAIC is big!
# > mod6
# Myo ~ fTreatment * fMonth + fYear + fCity + s(TimeSunset, k = 10) + 
#    s(fDevice, bs = "re")
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




















# # LAMP? ----
# # without 100%  TRY TO DO WITHOUT MONTHS
# ggplot(newdat2, aes(x = fTreatment, y = fit))+#, colour = fMonth, group = fMonth)) +
#   geom_point(size = 5, position = position_dodge(width = 0.3)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr),
#                 width = 0., position = position_dodge(width = 0.3)) +
#   #geom_line(size = 1, position = position_dodge(width = 0.3)) +
#   ylab("Predicted probability of bat presence") +
#   #xlab("Lamp type") +
#   #  scale_colour_manual(values = c("Jun" = "#bf8c0a", "Aug" = "#503294")) +
#   theme_classic(base_size = 18) +
#   theme(legend.position = "top")# +
# # labs(colour = "Month")
# 
# 
# POST HOC ----

summary(bat_modM_AR1)
# Parametric coefficients:
#                             Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                -5.5670     0.4074 -13.666  < 2e-16 ***
#   fTreatmentLight            -1.4109     0.5353  -2.636  0.00839 ** 
#   fMonthJun                   0.3840     0.1950   1.969  0.04892 *  
#   fYear2023                  -0.2812     0.1716  -1.639  0.10131    
#   fCityH                      0.1357     0.5042   0.269  0.78775    
#   fTreatmentLight:fMonthJun   0.5590     0.3803   1.470  0.14160    
# ---
# Approximate significance of smooth terms:
#                       edf   Ref.df Chi.sq  p-value    
#    s(TimeSunset)    6.500    7.623  35.01 1.58e-05 ***
#    s(Temp_mean_z)   6.512    7.555  45.58  < 2e-16 ***
#    s(fNight)      105.535 1006.000 237.12  < 2e-16 ***
#    s(fDevice)      17.286   23.000 312.83  < 2e-16 ***

# Treatment
emm_treat<- emmeans(bat_modM_AR1, ~ fTreatment , type = "response")
plot(emm_treat)
cld(emm_treat, adjust = "tukey", Letters = letters) 
# fTreatment    prob       SE    df lower.CL upper.CL .group
# fMonth = Aug:
#   fTreatment    prob       SE    df lower.CL upper.CL .group
# Light      0.00105 0.000473 46130  0.00038  0.00287  a    
# Dark       0.00427 0.001670 46130  0.00179  0.01020   b   
# 
# fMonth = Jun:
#   fTreatment    prob       SE    df lower.CL upper.CL .group
# Light      0.00268 0.001140 46130  0.00103  0.00693  a    
# Dark       0.00626 0.002420 46130  0.00263  0.01484  a   
#' COMMENT: 
#' Predicted prob of myotis occurrence was 0.17% at illuminated sites
#' And 0.52% at unlit sites   WHEN LOOKING WITHOUT MONTH
#' COMMENT:
#' When looking at specifically both months, it can be seen that occurrence prob
#' is lower always for unlit sites, but this difference was not significant in June.
#' However, interaction of montxtreatment was not significant.
#' We retained the interaction in the model however, because it was central to our 
#' research questions and hypotheses

contr <- contrast(
  emm_treat,
  method = "trt.vs.ctrl",
  ref = "Dark",
  #by = "fMonth",
  adjust = "none"  
)
summary(contr, infer = TRUE)
# fMonth = Aug:
#  contrast     odds.ratio    SE    df lower.CL upper.CL null t.ratio p.value
#  Light / Dark      0.244 0.131 46130   0.0854    0.696    1  -2.636  0.0084
# 
# fMonth = Jun:
#  contrast     odds.ratio    SE    df lower.CL upper.CL null t.ratio p.value
#  Light / Dark      0.427 0.223 46130   0.1535    1.186    1  -1.633  0.1024
#
#' The odds of myotis occurrence at illumnated sites OVERALL
#' were approx 68% lower than at unlit sites (OR = 0.32, 95% CI = 0.12-0.85, p = 0.022).
#' 32.3% of the odds observed at dark sites.
#' 
#' Not significant interaction: limited evidence that the effect of illumination 
#' differed between June and August.
#' 
#' Light / Dark = 0.323, which means:
#' 
#' Odds at Light sites
#' ------------------- = 0.323
#' Odds at Dark sites          #  Denominator is here dark (reference)
#' 
#' Since 0.323 is less than 1:Light sites have lower odds than dark sites.
#' If you want to describe the advantage of Dark sites, you flip the ratio:
#' 1 / 0.323 = 3.10, meaning that The odds of Myotis occurrence were 
#' approximately 3.1 times higher at unlit sites than at illuminated sites.
#' 
#' For eptesicus it was LED / dark = 5.73
#'
#' Odds at LED sites
#' ----------------- = 5.73
#' Odds at dark sites        # #  Denominator was also dark (reference)
#'
#' Since 5.73 is greater than 1: LED sites have higher odds than dark sites.
#' No need to flip it. Meaning that The odds of Eptesicus occurrence were 
#' approximately 5.7 times higher at LED-lit sites than at unlit sites.






# B) LUX MODEL MYO ----

# Exploring effect of lux to bat presence
ggplot(
  subset(Dat_all, fLamp == "dark" & !is.na(Lux)),
  aes(log1p(Lux), Myo)
) + geom_smooth(method = "gam")

# Subset
Dat_dark <- droplevels(
  subset(Dat_all, fLamp == "dark")
)



# collinearity
Dat_dark$logLux<-log1p(Dat_dark$Lux)
cor(Dat_dark$TimeSunset, Dat_dark$logLux, use = "complete.obs") # -0.3085

bat_mod <- bam( Myo ~ fMonth + fYear + fCity +  s(log1p(Lux)) +s(TimeSunset, k = 10) + s(fDevice, bs = "re"), family = binomial, data = Dat_dark, method = "fREML", correlation = corAR1(form = ~ TimeSunset | fNight), discrete = TRUE )
res <- residuals(bat_mod, type = "pearson")
rho_est <- acf(res, plot = FALSE)$acf[2]
rho_est

# Model 
modM_dark_lux <- bam(
  Myo ~
    fMonth +
    fYear +
    fCity +
    s(TimeSunset, k = 15) +
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
summary(modM_dark_lux)
gam.check(modM_dark_lux) # 

# Autocorrelation
acf(residuals(modM_dark_lux, type = "pearson")) ; acf(residuals(modM_dark_lux, type = "pearson"), plot = FALSE)$acf[2] # 0.06 - 0.11 is ok for the second bar!
res2 <- residuals(modM_dark_lux, type = "pearson")

# DHARMa
sim_res <- simulateResiduals(modM_dark_lux, n = 1000);plot(sim_res)
testZeroInflation(sim_res) #ok


# 'Multicollinearity' 
concurvity(modM_dark_lux)
#' It seems to be impossible to separate the effect of TimeSunset, Lux (and Temperature, but very little),
#' since they change together. Therefore we can't conclude wether the bats are more or
#' less likely to occur under specific lux conditions.
#' However, here we can conclude that lux levels in June are very different from august,
#' and that also at dark sites ept are more present in June than August.

# Partial effect
draw(modM_dark_lux)

## COMPARISON between smothers ----

modM_dark_noLux <- bam(
  Myo ~
    fMonth +
    fYear +
    fCity +
    s(TimeSunset, k = 10) +
  #  s(log1p(Lux)) +
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


AIC(modM_dark_lux, modM_dark_noTemp, modM_dark_noLux)
concurvity(modM_dark_noTemp)
summary(modM_dark_noTemp)
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
  modM_dark_lux,
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




