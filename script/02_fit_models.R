# Fit models
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(lubridate)
library(gammit)
library(mgcv)

load("data/paris_decay100221.rda")
load("data/paris_decaycomp100221.rda")
load("data/sero_rev.rda")


# Antibody decay ----
## Marginal -- absolute titers ----
mod1 <- gam(spk_end2 ~ 1 + s(days_post_onset, k=7)+ s(id, bs="re") + 
              sex + age_d , 
            data=paris_sub5, family = gaussian("identity"), method = "REML")

summary(mod1)
gtsummary::tbl_regression(mod1)
# Exponentiated (log2) coefs and 95% CI
# Age 40+: 1.69 (1.27, 2.30)
# Sex F:  1.44 (1.06, 2.0)

chk <- distinct(paris_sub5)

# gam.check(mod1)
# diag(vcov.gam(mod1))
min(paris_sub5$days_post_onset)#0
max(paris_sub5$days_post_onset)#406

newdat_m1 <- data.frame(days_post_onset = rep(0:406,10),
                        sex = factor(c(rep("M", 2035),
                                       rep("F", 2035)), ordered = F),
                        age_d = factor(c(rep("1",2035),
                                         rep("2",2035)), ordered = F),
                        id = factor(rep("60",4070),ordered = F))

p_m1 <- predict_gamm(mod1, newdata = newdat_m1, se = T,re_form = NA)
newdat2_m1 <- data.frame(newdat_m1, p_m1)%>%
  mutate(lower = prediction - 1.96*se,
         upper = prediction + 1.96*se)

newdat2_m1 <- data.frame(newdat_m1, p_m1)%>%
  mutate(lower = prediction - 1.96*se,
         upper = prediction + 1.96*se)

newdat3_m1 <- newdat2_m1%>%
  ungroup()%>%
  group_by(days_post_onset)%>%
  summarise(mn_pred = mean(prediction),
            mn_low = mean(lower),
            mn_up = mean(upper))%>%
  mutate(dif = 		
           11.26869 - mn_pred)

plot(newdat3_m1$days_post_onset, newdat3_m1$dif) # get way shorter halflife...

mn_start <- newdat3_m1%>%
  filter(days_post_onset==0)%>%
  summarise(mn_start= mn_pred)%>%
  select(mn_start)

save(newdat3_m1, file = "data/margdecay_pred.rda")


# Separate splines depending on baseline titer -- absolute titers ----

mod2 <- gam(spk_end2 ~ 1 + s(days_post_onset, k=7, by=spk_binf)+ s(id, bs="re") + 
               sex + age_d + spk_binf, 
             data=paris_sub6, family = gaussian("identity"), method = "REML")

summary(mod2)
# gam.check(mod2)

gtsummary::tbl_regression(mod2)
# Exponentiated (log2) coefs and 95% CI
# Age 40+: 1.69 (1.27, 2.30)
# Sex F:  1.44 (1.06, 2.0)

newdat_m2 <- data.frame(days_post_onset = rep(0:406,10),
                        spk_binf = factor(c(rep("1",2035),
                                                rep("2",2035)), ordered = F),
                        sex = factor(c(rep("M", 2035),
                                       rep("F", 2035)), ordered = F),
                        age_d = factor(c(rep("1",2035),
                                         rep("2",2035)), ordered = F),
                        id = factor(rep("60",4070),ordered = F))

p_m2 <- predict_gamm(mod2, newdata = newdat_m2, se = T,re_form = NA)
newdat2_m2 <- data.frame(newdat_m2, p_m2)%>%
  mutate(lower = prediction - 1.96*se,
         upper = prediction + 1.96*se)


newdat3_m2<- newdat2_m2%>%
  ungroup()%>%
  group_by(days_post_onset, spk_binf)%>%
  summarise(mn_pred = mean(prediction),
            mn_low = mean(lower),
            mn_up = mean(upper))

save(newdat3_m2, file = "data/compdecay_pred.rda")


# Seroreversion----
library(survival)
library(survminer)

## updated definition
km_fit <- survfit(Surv(days_post_onset, as.numeric(sero_rev2)) ~ as.factor(spk_binf), data=sero_revtot)
km <- survdiff(Surv(days_post_onset, as.numeric(sero_rev2)) ~ as.factor(spk_binf), data=sero_revtot)
km

km_base <- ggsurvplot(
  fit = km_fit, 
  pval = T,
  fun = "pct",
  xlab = "Days",
  legend.title = "",
  legend.labs = c("Baseline titer <800", "Baseline titer 800+"),
  conf.int = T,
  censor = T,
  palette = c("#440154FF","#5DC863FF"),
  risk.table = F,
  ggtheme = theme_classic2(base_size=16)
)

km_base  
# 
# ## initial definition
# km_fit1 <- survfit(Surv(days_post_onset, as.numeric(sero_rev1)) ~ as.factor(spk_binf), data=sero_revtot)
# km1 <- survdiff(Surv(days_post_onset, as.numeric(sero_rev1)) ~ as.factor(spk_binf), data=sero_revtot)
# km1
# 
# km_base1 <- ggsurvplot(
#   fit = km_fit1, 
#   pval = T,
#   fun = "pct",
#   xlab = "Days",
#   legend.title = "",
#   legend.labs = c("Baseline titer <800", "Baseline titer 800+"),
#   # ylab = "Hazard of seroreversion",
#   # ylab = "Survival probability",
#   conf.int = T,
#   # surv.median.line = "hv",
#   censor = T,
#   palette = c("#440154FF","#5DC863FF"),
#   risk.table = F,
#   ggtheme = theme_classic2(base_size=16)
# )
# 
# km_base1



