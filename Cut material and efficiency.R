setwd("...to experiment_data.csv...")
library(tidyverse)
library(lme4)
library(MuMIn)
library(chngpt)

## Load data
dat <- read_csv("experiment_data.csv")


## Summary stats for different materials
time_summary <- dat %>% dplyr::group_by(material) %>% summarize(Mean_time_s = mean(Time_s), SD = sd(Time_s), Time_min_s = min(Time_s), Time_max_s = max(Time_s)) %>% arrange(desc(Mean_time_s))
write.csv(time_summary, "time_by_material.csv")

## Main plots

dat <- dat %>% separate(`plan/profile`, c("Plan_cat", "Profile_cat"), remove = F)
dat <- dat %>% add_column("Plan Curvature" = ifelse(dat$Plan_cat == "high", "High", "Low"))
ggplot(dat, aes(x=sEdge_length_mm, y=Time_s, color=`Plan Curvature`, shape = `Plan Curvature`)) + facet_wrap(~material, ncol = 2, scales = 'free') +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "red"))+
  geom_point(size = 2.2)+ geom_smooth(method = "lm", size = 0.5) + theme_bw() +
  theme(plot.title = element_text(size = 18, face="bold", hjust = 0.5), 
        axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.position = c(0.82, 0.93),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Edge length (Z-score)", y = "Time (s)")


dat <- dat %>% add_column("Profile Curvature" = ifelse(dat$Profile_cat == "high", "High", "Low"))
ggplot(dat, aes(x=sEdge_length_mm, y=Time_s, color=`Profile Curvature`, shape = `Profile Curvature`)) + facet_wrap(~material, ncol = 2, scales = 'free') +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "red"))+
  geom_point(size = 2.2)+ geom_smooth(method = "lm", size = 0.5) + theme_bw() +
  theme(plot.title = element_text(size = 18, face="bold", hjust = 0.5), 
        axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.position = c(0.82, 0.93),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Edge length (Z-score)", y = "Time (s)")

ggplot(dat, aes(x=`artifact size group (cm)`, y=Time_s, color=`Profile Curvature`, shape = `Profile Curvature`)) + facet_wrap(~material, ncol = 2, scales = 'free') +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "red"))+
  geom_boxplot()+ #geom_smooth(method = "lm", size = 0.5) + 
  geom_jitter() +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face="bold", hjust = 0.5), 
        axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.position = c(0.87, 0.93),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Specimen size (cm)", y = "Time (s)")

ggplot(dat, aes(x=`artifact size group (cm)`, y=Time_s, color=`Plan Curvature`, shape = `Plan Curvature`)) + facet_wrap(~material, ncol = 2, scales = 'free') +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "red"))+
  geom_boxplot()+ #geom_smooth(method = "lm", size = 0.5) + 
  geom_jitter() +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face="bold", hjust = 0.5), 
        axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.position = c(0.87, 0.93),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Specimen size (cm)", y = "Time (s)")


### ANALYSIS


## Interaction models
dat$session <- as.factor(dat$session)
dat$sPlan <- scale(dat$`Plan index`)[,1]
dat$sProfile <- scale(dat$`profile index`)[,1]

glm_fit_size_plan <- glmer(Time_s ~ `artifact size group (cm)`*material*sPlan + (1|participant) + (1|session), data = dat, family = "poisson")
a <- summary(glm_fit_size_plan)
write.csv(round(a$coefficients,2), "glm_size_plan.csv")

glm_fit_size_profile <- glmer(Time_s ~ `artifact size group (cm)`*material*sProfile + (1|participant) + (1|session), data = dat, family = "poisson")
a <- summary(glm_fit_size_profile)
write.csv(round(a$coefficients, 2), "glm_size_profile.csv")


glm_fit_edgelength_plan <- glmer(Time_s ~ sEdge_length_mm*material*sPlan + (1|participant) +(1|session), data = dat, family = "poisson")
a <- summary(glm_fit_edgelength_plan)
write.csv(round(a$coefficients, 2), "glm_edgelength_plan.csv")

glm_fit_edgelength_profile <- glmer(Time_s ~ sEdge_length_mm*material*sProfile + (1|participant) + (1|session), data = dat, family = "poisson")
a <- summary(glm_fit_edgelength_profile)
write.csv(round(a$coefficients, 2), "glm_edgelength_profile.csv")


# Thresholds

fit_size=chngptm(formula.1=Time_s~1, formula.2=~artifact_size_num, family="poisson", data=dat, 
                 type="segmented", var.type="bootstrap")

summary(fit_size)
plot(fit_size)

test_size=chngpt.test(formula.null=log(Time_s)~1, formula.chngpt=~artifact_size_num, dat, family="gaussian", type="segmented")
test_size


fit_edge_length=chngptm(formula.1=Time_s~1, formula.2=~edge_length, family="poisson", data=dat, 
                        type="segmented", var.type="bootstrap")

summary(fit_edge_length)
plot(fit_edge_length)

test_length=chngpt.test(formula.null=log(Time_s)~1, formula.chngpt=~edge_length, dat, family="gaussian", type="segmented")
test_length
