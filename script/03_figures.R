# Figures 

# AUC by age violin plots----
base_auc%>%
  mutate(Age = case_when(
    age_d==1 ~ "<40 years",
    age_d==2 ~ "40+ years"
  ))%>%
  ggplot(aes(x = Age, y = auc)) + 
  geom_violin()+
  geom_boxplot(width = .15, outlier.shape = NA) +
  geom_jitter(aes(colour = Age),shape=16, position=position_jitter(width=0.2, height = 0), alpha = 0.3)+
  # scale_color_viridis_d(end = 0.75, labels = c("<40 years", "40+ years"), name = "Age")+
  scale_color_manual(values = c("#440154FF", "#7AD151FF"))+
  # coord_cartesian(xlim = c(1.2, NA), clip = "off")+
  scale_y_continuous(trans = "log", breaks = c(7, 50, 400, 3000))+
  theme_bw()+
  labs(y ="Area under the curve")+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 16),
        # axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        # axis.title.x = element_text(size = 20)
        axis.title = element_text(size = 20))


labs_y <- c("50", "80", "100", "200", "400", "800", "1600", "3200", "6400")
l2brks <- log2(as.numeric(labs_y))

# Antibody decay -- absolute titers -- marginalized over baseline titer ----
load("data/margdecay_pred.rda")

comp_marg_ptlg <-  ggplot(aes(x=days_post_onset,y=mn_pred), data=newdat3_m1)+
  geom_ribbon(aes(ymin = mn_low, ymax=mn_up), alpha=0.5, fill="#2C728EFF")+
  # scale_fill_viridis_d(end=0.8 , name="Baseline titer", labels = c("<800", "800+"))+
  geom_line(color="dark gray")+
  geom_point(inherit.aes = F, data = paris_sub5, aes(x=days_post_onset, y=spk_end2), alpha=0.5, size=1.3)+
  # scale_colour_viridis_d(end=0.8,  name="Baseline titer", labels = c("<800", "800+"))+
  labs(x = "Days since illness", y="Spike Endpoint")+
  # theme_classic()+
  theme_bw()+
  geom_hline(yintercept = log2(80), linetype=2)+
  scale_y_continuous(trans = 'log2',limits =c(5.5,13), breaks = l2brks, labels = labs_y)+
  theme(axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20))
comp_marg_ptlg

# Antibody decay -- absolute titers -- comparing baseline titer ----
load("data/compdecay_pred.rda")

comp_base_ptlg <- ggplot(aes(x=days_post_onset,y=mn_pred, group=spk_binf), data=newdat3_m2) +
  geom_ribbon(aes(ymin = mn_low, ymax=mn_up, fill=spk_binf), alpha=0.5)+
  scale_fill_viridis_d(end=0.8 , name="Baseline titer", breaks = c(1,2), labels = c("<800", "800+"))+
  geom_line(color="dark gray")+
  geom_point(inherit.aes = F, data = paris_sub5, aes(x=days_post_onset, y=spk_end2, colour=spk_binf), alpha=0.5, size=1.3)+
  scale_colour_viridis_d(end=0.8,  name="Baseline titer", breaks = c(1,2), labels = c("<800", "800+"))+
  labs(x = "Days since illness", y="")+
  # theme_classic()+
  theme_bw()+
  geom_hline(yintercept = log2(80), linetype=2)+
  scale_y_continuous(trans = 'log2',limits =c(5.5,13), breaks = l2brks, labels = labs_y)+
  theme(axis.text.y = element_text(size = 18),
        # axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = "top")
comp_base_ptlg





