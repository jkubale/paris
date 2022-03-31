# Prep data for analysis

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(lubridate)
library(gammit)

# Read data ----
paris_long <- read_excel("data/paris_long5.xlsx")
paris_part <- read_excel("data/paris_part4.xlsx")
cs_abinfo <- read_excel("data/cs_abinfo.xlsx")


## recode any titer < 80 as 50 (negative)
paris_long$spk_end1 <- as.numeric(ifelse(paris_long$spk_end=="NEG"|
                                           paris_long$spk_end=="Neg"|
                                           paris_long$spk_end=="neg"|
                                           paris_long$spk_end=="Negative"|
                                           paris_long$spk_end=="25"|
                                           paris_long$spk_end=="40"|
                                           paris_long$spk_end=="50",
                                         50,paris_long$spk_end))

paris_long <- paris_long[,-8] # drop original spk_end variable 

## drop unneeded variables
paris_sub <- paris_long%>%
  filter(is.na(spk_end1)==F)%>%
  group_by(id)%>%
  mutate(id = as.factor(id))%>%
  rename(spk_end = spk_end1)%>%
  select(id, day_postbase, sampid, spk_end, auc)

## Identify participants who ever had positive antibodies
have_titer <- paris_sub%>%
  group_by(id)%>%
  mutate(detect = case_when(
    spk_end >50 ~1,
    T~0
  ))%>%
  filter(detect==1)%>%
  distinct(id)

## Subset to only those who had antibodies at some point in the study
paris_sub2 <- left_join(have_titer, paris_sub, by="id")



# Add those that had earlier qualitative antibody 
cs_abinfo <- cs_abinfo%>%
  mutate(onset = case_when(
    is.na(pcr_date)==F~pcr_date,
    is.na(pcr_date) & is.na(symp_date)==F ~ symp_date,
    is.na(pcr_date) & is.na(symp_date) & is.na(ab_date)==F ~ ab_date
  ))%>%
  filter(is.na(onset)==F)%>%
  select(id, onset)

# Any missing from pariticipant data?
chk <- anti_join(cs_abinfo, paris_part, by = "id") # 4 participants missing info id:419, 420, 421, 429 -- will drop for now

paris_part2 <- left_join(paris_part, cs_abinfo, by="id")%>%
  mutate(onset2 = case_when(
    is.na(onset)==F ~ onset,
    is.na(pcr_date_pre)==F~pcr_date_pre,
    is.na(pcr_date_pre) & is.na(pcr_date_post)==F ~ pcr_date_post,
    is.na(pcr_date_pre) & is.na(pcr_date_post) & is.na(symp_onset)==F ~ symp_onset
  ),
  id = as.factor(id))%>%
  filter(is.na(onset)==F, age!="Missing")%>%
  mutate(onset3 = case_when(
    onset2 > entry ~ entry,
    T~ onset2
  ))%>%
  select(-c(onset, onset2))%>%
  rename(onset = onset3)


chk2 <- filter(paris_part2, onset > entry) # should be 0


## add covariates to long form dataset
paris_sub3 <- inner_join(paris_sub2, paris_part2, by = "id") ##142 participants, 812 obs

## Use onset to adjust days
## also creates variable for detect/not detect to use in identifying seroreversion
paris_sub3 <- paris_sub3%>%
  group_by(id)%>%
  mutate(day_ad = as.numeric(entry - onset),
  days_post_onset = day_postbase+ day_ad,
  detect = case_when(
    spk_end >50 ~1,
    T~0
  ))

## Identify first detectable titer for each person and remove any obs before that
first_titer <- paris_sub3%>%
  filter(detect==1)%>%
  group_by(id)%>%
  summarise(mn_tdate = min(day_postbase))

paris_sub4 <- left_join(paris_sub3, first_titer, by="id")

# remove obs before participants became positive and remove those with only 1 obs -- gives 130 participants, 763 obs
paris_sub4 <- paris_sub4%>%
  mutate(del = case_when(
    day_postbase < mn_tdate ~1,
    T~0
  ))%>%
  filter(del==0)%>%
  select(-del)%>%
  filter(age != "Missing")%>%
  group_by(id)%>%
  mutate(obs = row_number(),
         age_d = case_when(
           age == "18-29"|age=="30-39"~1,
           T~2
         ))%>%
  filter(max(obs)>1)%>%
  select(-obs)

base <- paris_sub4%>%
  group_by(id)%>%
  arrange(id, day_postbase)%>%
  mutate(obs = row_number())%>%
  filter(obs==1, is.na(spk_end)==F)%>%
  mutate(base_spk = spk_end,
  spk_bin = case_when(
    base_spk <800 ~1, 
    base_spk >=800 ~2
  ))%>%
  select(id, base_spk, spk_bin)

paris_sub5 <- left_join(paris_sub4, base, by="id")

# finalize variables for regression models and fix any onset dates
paris_sub5 <- paris_sub5%>%
  mutate(sex = as.factor(sex),
         age_d = as.factor(age_d),
         spk_bin = as.factor(spk_bin), 
         spk_end2 = log2(spk_end),
         spk_binf = factor(spk_bin, ordered = F, levels = c("1","2")))

# save(paris_sub5, file = "data/paris_decay100221.rda")
# for model comparing decay by baseline titer
paris_sub6 <- filter(paris_sub5, day_postbase>0)
# save(paris_sub6, file = "data/paris_decaycomp100221.rda")

# Seroreversion ----
## need to fix merging using file from Charles

sero_rev <- read_excel("data/sero_rev.xlsx")%>%
  filter(sero_rev2==2 | sero_rev1==2)%>%
  mutate(id = as.factor(id), 
         spk_binf = as.factor(base_spk),
         sero_rev1 = as.factor(sero_rev1),
         sero_rev2 = as.factor(sero_rev2))%>%
  select(id, spk_binf, days_post_onset, sero_rev1, sero_rev2)
  

no_serorev <- anti_join(paris_sub5, sero_rev, by="id")%>%
  group_by(id)%>%
  summarise(days_post_onset = max(days_post_onset))

covars <- paris_sub5%>%
  distinct(id, .keep_all = T)%>%
  select(id, spk_binf)

no_serorev2 <- left_join(no_serorev, covars, by = "id")%>%
  mutate(sero_rev1 = as.factor(1),
         sero_rev2 = as.factor(1))%>%
  select(id, spk_binf, days_post_onset, sero_rev1, sero_rev2)

sero_revtot <- bind_rows(sero_rev, no_serorev2)

# save(sero_revtot, file = "data/sero_rev.rda")

## Baseline auc----
base_auc <- paris_sub5%>%
  filter(day_postbase==0)

# save(base_auc, file = "data/base_auc.rda")
