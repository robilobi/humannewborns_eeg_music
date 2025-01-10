library(lme4)
library(emmeans)
library(effects) #for lmer
library(ggplot2)
library(car)
library(sjPlot)
library(webshot)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(lubridate)
library(plyr)
library(lubridate)
library(dplyr)
library(ez)
library(schoRsch)  #for formatted output of anova


ls()
rm(list=ls())
dir = "C:/Users/robianco/OneDrive - Fondazione Istituto Italiano Tecnologia/BACHNB/RESULTS_NB/TRF_GT_both+_rnd_sosp/"

setwd(dir)
data = read.csv(paste0(dir, "PredAccTRFforstat.csv"), header = F, sep = ',') 
colnames(data) <- c("sbj", "PredAcc","condition", "model", "melID", "PredAccFull", 'roi') 

data$condition = as.factor(data$condition)
data$model = as.factor(data$model)
data$sbj= as.factor(data$sbj)
data$roi= as.factor(data$roi)
data$melID= as.factor(data$melID)
data$condition = fct_recode(data$condition, "Real" = "1", "Shuffled" = "2")
data$model = fct_recode(data$model, "-so" = '1', "-sp" = '2', "-ioi" = "3", "-iti" = "4", "-mus" = "5")
data$melID = fct_recode(data$melID, "s_01" = "11","s_05" = "12","s_08" = "13","s_10" = "14",
                        "o_05"= "5", "o_08"= "8", "o_03"= "3", "o_06"= "6", "o_04"= "4",
                        "o_09"= "9", "o_02"= "2", "o_01"= "1", "o_07"= "7", "o_10"= "10")
summary(data)
length(unique(data$sbj))
orginaldata = data

ag=with(data, aggregate(cbind(PredAcc ) ~sbj+model+condition, FUN="mean"))
ggplot(ag, aes(x=model, y=PredAcc, color = model, fill = model)) +theme_classic()+
  facet_grid(~condition)+ 
  stat_summary(fun=mean, geom="point", size=4, position=position_dodge(width=0.9))+#coord_cartesian(ylim=c(-0.001,0.001))+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", width=0, alpha=1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at zero
  xlab("") + ylab("Delta r (Full-Reduced)")+theme(legend.position = 'bottom')+
  theme(legend.position = 'bottom',
    text = element_text(size = 16),  # Increase font size
    axis.text = element_text(size = 14),  # Adjust axis text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    strip.text = element_text(size = 14))
ggsave(paste0("PredAccFULLmodel_TRF_ALLCond.pdf"), width = 5, height = 4)

####FOCUS ON MUS RANDOM 
data <- orginaldata %>%
  filter(model %in% c("-mus"))
summary(data)
####### plot full model prediction accuracy (no diff between condition)
m.1= lmer(PredAcc~condition  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum'),
          data)
Anova(m.1, type=3)  #
emt = emmeans(m.1, pairwise~condition)
summary(emt, infer = c(TRUE, TRUE)) #test if different from zero

ag=with(data, aggregate(cbind(PredAcc ) ~sbj+condition, FUN="mean"))
ggplot(ag, aes(x=condition, y=PredAcc, color = condition, fill = condition)) +theme_classic()+
  scale_color_manual( values=c("red", "grey"))+
  scale_fill_manual( values=c("red", "grey"))+
  stat_summary(fun=mean, geom="point", size=4, position=position_dodge(width=0.9))+#coord_cartesian(ylim=c(-0.001,0.001))+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", width=0, alpha=1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at zero
  xlab("") + ylab("Delta r (Full-Reduced)")+theme(legend.position = 'bottom')+
  coord_cartesian(ylim = c(-0.002, 0.0045)) + # Zoom in on y-axis from 0 to 10
  theme(legend.position = 'bottom',
        text = element_text(size = 16),  # Increase font size
        axis.text = element_text(size = 14),  # Adjust axis text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        strip.text = element_text(size = 14))
ggsave(paste0("PredAccFULLmodel_TRF_Full_RvS.pdf"), width = 5, height = 4)
ggsave(paste0("PredAccFULLmodel_TRF_Full_RvS.png"), width = 5, height = 4)


ag=with(data, aggregate(cbind(PredAcc) ~sbj+condition+melID, FUN="mean"))
level_mel =c("o_05", "o_08", "o_03", "o_04", "o_09", "o_06", "o_02", "o_07", "o_01", "o_10", "s_05", "s_01", "s_08", "s_10")
ag <- ag %>% 
  arrange(factor(melID, levels = level_mel))
ag$condition=fct_relevel(ag$condition, "Shuffled", "Real")
ggplot(ag, aes(x =factor(melID, levels = level_mel), y=PredAcc, color = condition, fill = condition)) +theme_classic()+
  scale_color_manual( values=c( "grey", "red"))+
  scale_fill_manual( values=c("grey", "red"))+
  stat_summary(fun=mean, geom="point", size = 4) +
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", width = 0, alpha = 0.6) +
  xlab("") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ylab("Delta r") +
  theme(legend.position = 'bottom', text = element_text(size = 16),  # Increase font size
        axis.text = element_text(size = 14),  # Adjust axis text size
        axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x-axis labels
        axis.title = element_text(size = 16),  # Adjust axis title size
        strip.text = element_text(size = 14)  # Adjust facet label size
  )
ggsave(paste0("PredAccFULLmodel_TRF_byMelCond_AM.pdf"), width = 6, height = 4)
ggsave(paste0("PredAccFULLmodel_TRF_byMelCond_AM.png"), width = 6, height = 4)


ag=with(data, aggregate(cbind(PredAccFull) ~sbj+condition+melID, FUN="mean"))
level_mel =c("o_05", "o_08", "o_03", "o_04", "o_09", "o_06", "o_02", "o_07", "o_01", "o_10", "s_05", "s_01", "s_08", "s_10")
ag <- ag %>% 
  arrange(factor(melID, levels = level_mel))
ag$condition=fct_relevel(ag$condition, "Shuffled", "Real")
ggplot(ag, aes(x =factor(melID, levels = level_mel), y=PredAccFull, color = condition, fill = condition)) +theme_classic()+
  scale_color_manual( values=c( "grey", "red"))+
  scale_fill_manual( values=c("grey", "red"))+
  stat_summary(fun=mean, geom="point", size = 4) +
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", width = 0, alpha = 0.6) +
  xlab("") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ylab("r FULL") +
  theme(legend.position = 'bottom', text = element_text(size = 16),  # Increase font size
        axis.text = element_text(size = 14),  # Adjust axis text size
        axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x-axis labels
        axis.title = element_text(size = 16),  # Adjust axis title size
        strip.text = element_text(size = 14)  # Adjust facet label size
  )
ggsave(paste0("PredAccFULLmodel_TRF_byMel_FULL.pdf"), width = 6, height = 4)
ggsave(paste0("PredAccFULLmodel_TRF_byMel_FULL.png"), width = 6, height = 4)


####RANK FULL and GAIN
am = data
#am=am %>% filter(condition %in% c("Real"))

am$gain = am$PredAcc
am$full = am$PredAcc+am$PredAccFull
am$mean.am <- ave(am$full, am$sbj, FUN = function(x) mean(x, na.rm = TRUE))
am$diff.am <- ave(am$gain, am$sbj, FUN = function(x) mean(x, na.rm = TRUE))
am.maps = unique(am[,c("sbj", "mean.am")])
am.maps$rank = rank(-am.maps$mean.am)
am <- am %>% select(-roi)
m = am %>%
  pivot_longer(! c(sbj, condition, model, melID, mean.am,diff.am))
m=droplevels(m[m$name!='PredAccFull',])
m=droplevels(m[m$name!='PredAcc',])
m = merge(m, am.maps[,c("sbj", "rank")], all.x=T)

# ggplot(m, aes(rank, value, group = name, col = name))+
# #  facet_wrap(~condition)+
#   stat_summary(fun=mean, geom="point")+theme_classic()+
#   scale_color_manual( values=c( "black", "red"))+
#   scale_fill_manual( values=c( "red", "red"))+
#  # stat_summary(aes(x=rank, y=diff.am), fun=mean, geom="bar", fill = "red", color='white')+
#   coord_cartesian(ylim =c(-0.007,0.076))+theme(legend.position="bottom")+
#   stat_summary( fun.data = 'mean_se', geom = "errorbar", width=0, alpha=0.6, color = 'black')
# ggsave("PredAccFULLmodel_rankbysubj_AM.bmp", width = 6, height = 4)
# ggsave("PredAccFULLmodel_rankbysubj_AM.pdf", width = 6, height = 4)

ggplot(m[m$name =='full',], aes(rank, value, group = name, col = name))+
  #  facet_wrap(~condition)+
  stat_summary(fun=mean, geom="point", color ='black')+theme_classic()+
  coord_cartesian(ylim =c(-0.007,0.076))+theme(legend.position="bottom")+
  stat_summary( fun.data = 'mean_cl_boot', geom = "errorbar", width=0, alpha=0.6, color = 'black')
ggsave("PredAccFULLmodel_rankbysubj_AM.bmp", width = 6, height = 4)
ggsave("PredAccFULLmodel_rankbysubj_AM.pdf", width = 6, height = 4)




#####rank by rdiff
am = data
am$gain = am$PredAcc
am$full = am$PredAcc+am$PredAccFull
am$mean.am <- ave(am$full, am$sbj, FUN = function(x) mean(x, na.rm = TRUE))
am$diff.am <- ave(am$gain, am$sbj, FUN = function(x) mean(x, na.rm = TRUE))
am.maps = unique(am[,c("sbj", "diff.am")])
am.maps$rank = rank(-am.maps$diff.am)
am <- am %>% select(-roi)
m = am %>%
  pivot_longer(! c(sbj, condition, model, melID, mean.am,diff.am))
m=droplevels(m[m$name!='PredAccFull',])
m=droplevels(m[m$name!='PredAcc',])
m = merge(m, am.maps[,c("sbj", "rank")], all.x=T)

ggplot(m[m$name =='gain',], aes(rank, value, group = condition, col = condition))+
  facet_grid(condition ~ .)+
  stat_summary(fun=mean, geom="point")+theme_classic()+
  stat_summary( fun.data = 'mean_cl_boot', geom = "errorbar", width=0, alpha=0.6)+
  scale_color_manual( values=c( "red", "grey"))+
  scale_fill_manual( values=c( "red", "grey"))+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  
  # stat_summary(aes(x=rank, y=diff.am), fun=mean, geom="bar", fill = "red", color='white')+
  coord_cartesian(ylim =c(-0.02,0.027))+theme(legend.position="bottom")
ggsave("PredAccFULLmodel_rankbysubj_AMRealvsSHU.bmp", width = 6, height = 4)
ggsave("PredAccFULLmodel_rankbysubj_AMRealvsSHU.pdf", width = 6, height = 4)



##singles subject across real and shuffled
ag=with(data, aggregate(cbind(PredAcc) ~sbj+condition, FUN="mean"))
dw <- ag %>%
  pivot_wider(
    names_from = condition,   # Make columns based on `model`
    values_from = PredAcc # Use `PredAcc` values for the new columns
  )

corr_value <- cor.test(dw$`Real`, dw$`Shuffled`, method = "spearman")
ggplot(dw, aes(x = `Real`  , y = `Shuffled` )) +
  geom_point() +  theme_classic()+
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.01, 0.02),xlim = c(-0.01, 0.03)) + 
  labs(title = "AM gain", x = "REAL", y = "SHUFFLED")
ggsave(paste0("GAIN_RvsS.pdf"), width = 5, height = 4)
ggsave(paste0("GAIN_RvsS.png"), width = 5, height = 4)




####FOCUS ON BOTH + AND MUSICAL PREDICTORS
data <- orginaldata %>%
  filter(model %in% c("-so", "-sp"))
summary(data)

desired_order <- c( "-so", "-sp")
data <- data %>%
  mutate(model = factor(model, levels = desired_order))

####### plot full model prediction accuracy (no diff between condition)
m.1= lmer(PredAcc~condition*model  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum', model='contr.sum'),
          data)
Anova(m.1, type=3)  #
emt = emmeans(m.1, pairwise~model*condition)
summary(emt, infer = c(TRUE, TRUE)) #test if different from zero
tab_model(m.1, file = "stats_TRF.doc")

plot_model(m.1, type = "pred", terms = c( "condition", "model"))+theme_classic()
ag=with(data, aggregate(cbind(PredAcc ) ~sbj+model+condition, FUN="mean"))
ggplot(ag, aes(x=model, y=PredAcc, color = model, fill = model)) +theme_classic()+
  facet_grid(~condition)+  
  scale_color_manual( values=c("forestgreen", "goldenrod1"))+
  scale_fill_manual( values=c("forestgreen", "goldenrod1"))+
  stat_summary(fun=mean, geom="point", size=4, position=position_dodge(width=0.9))+#coord_cartesian(ylim=c(-0.001,0.001))+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", width=0, alpha=1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at zero
  xlab("") + ylab("Delta r (Full-Reduced)")+theme(legend.position = 'bottom')+
  coord_cartesian(ylim = c(-0.002, 0.0045)) + # Zoom in on y-axis from 0 to 10
theme(legend.position = 'bottom',
    text = element_text(size = 16),  # Increase font size
    axis.text = element_text(size = 14),  # Adjust axis text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    strip.text = element_text(size = 14))
 ggsave(paste0("PredAccFULLmodel_TRF_byCond.pdf"), width = 5, height = 4)
ggsave(paste0("PredAccFULLmodel_TRF_byCond.png"), width = 5, height = 4)


####SPLIT BY MODEL to TEST CONDITION 
d <- filter(data, model %in% c("-so"))
m.1= lmer(PredAcc~condition  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum'),
          d)
Anova(m.1, type=3)  #
emt = emmeans(m.1, pairwise~condition)
summary(emt, infer = c(TRUE, TRUE)) #test if different from zero


d <- filter(data, model %in% c("-sp"))
m.1= lmer(PredAcc~condition  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum'),
          d)
Anova(m.1, type=3)  #
emt = emmeans(m.1, pairwise~condition)
summary(emt, infer = c(TRUE, TRUE)) #test if different from zero


#####CHECK DISTRIBUTION
ag=with(data, aggregate(cbind(PredAcc) ~sbj+model+condition, FUN="mean"))

# Plot histograms by condition
ggplot(ag, aes(x = PredAcc, fill = model)) +
  geom_histogram(binwidth = 0.001,alpha = 0.7, position = "identity") +
  facet_wrap(~ condition+model) + # Create separate panels for each condition
  scale_color_manual( values=c("forestgreen", "goldenrod1"))+
  scale_fill_manual( values=c("forestgreen", "goldenrod1"))+
  coord_cartesian(xlim = c(-0.015, 0.015)) + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  labs(title = "",x = "Delta r", y = "Frequency",fill = "Condition") +
  theme_minimal()
ggsave(paste0("Distribution_SpSo_REALSHUFFLED.png"), width = 5, height = 4)
ggsave(paste0("Distribution_SpSo_REALSHUFFLED.pdf"), width = 5, height = 4)


###Correlations real
d = ag[ag$condition=='Real',]
dw <- d %>%
  pivot_wider(
    names_from = model,   # Make columns based on `model`
    values_from = PredAcc # Use `PredAcc` values for the new columns
  )

corr_value <- cor.test(dw$`-so`, dw$`-sp`, method = "spearman")
ggplot(dw, aes(x = `-so`  , y = `-sp` )) +
  geom_point() +  theme_classic()+
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  
  annotate("text", x = 0, y = 0.01, label = paste("rho = ", round(corr_value$estimate, 2), "p =", 
            round(corr_value$p.value, 2)), size = 5, color = "red") +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.01, 0.01),xlim = c(-0.01, 0.012)) + 
  labs(title = "Correlation so  and sp Delta r (REAL)", x = "Delta so", y = "Delta sp")
ggsave(paste0("Correlation_SpSo_REAL.pdf"), width = 5, height = 4)
ggsave(paste0("Correlation_SpSo_REAL.png"), width = 5, height = 4)

##shuff
d = ag[ag$condition=='Shuffled',]
dw <- d %>%
  pivot_wider(
    names_from = model,   # Make columns based on `model`
    values_from = PredAcc # Use `PredAcc` values for the new columns
  )

corr_value <- cor.test(dw$`-so`, dw$`-sp`, method = "spearman")
ggplot(dw, aes(x = `-so`  , y = `-sp` )) +
  geom_point() +  theme_classic()+
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  
  annotate("text", x = 0, y = 0.01, label = paste("rho = ", round(corr_value$estimate, 2), "p =", 
                                                  round(corr_value$p.value, 2)), size = 5, color = "red") +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.01, 0.01),xlim = c(-0.01, 0.012)) + 
  labs(title = "Correlation so  and sp Delta r (SHUFFLED)", x = "Delta so", y = "Delta sp")
ggsave(paste0("Correlation_SpSo_SHUFF.pdf"), width = 5, height = 4)
ggsave(paste0("Correlation_SpSo_SHUFF.png"), width = 5, height = 4)



###########FOCUS ON LOW LEVEL
data <- orginaldata %>%
  filter(model %in% c("-ioi", "-iti"))
summary(data)

desired_order <- c( "-ioi", "-iti")
data <- data %>%
  mutate(model = factor(model, levels = desired_order))

####### plot full model prediction accuracy (no diff between condition)
m.1= lmer(PredAcc~condition*model  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum', model='contr.sum'),
          data)
Anova(m.1, type=3)  #
emt = emmeans(m.1, pairwise~model*condition)
summary(emt, infer = c(TRUE, TRUE)) #test if different from zero
plot_model(m.1, type = "pred", terms = c( "condition", "model"))+theme_classic()
ag=with(data, aggregate(cbind(PredAcc ) ~sbj+model+condition, FUN="mean"))
ggplot(ag, aes(x=model, y=PredAcc, color = model, fill = model)) +theme_classic()+
  facet_grid(~condition)+  
  scale_color_manual( values=c("lightblue", "pink"))+
  scale_fill_manual( values=c("lightblue", "pink"))+
  stat_summary(fun=mean, geom="point", size=4, position=position_dodge(width=0.9))+#coord_cartesian(ylim=c(-0.001,0.001))+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", width=0, alpha=1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at zero
  xlab("") + ylab("Delta r (Full-Reduced)")+theme(legend.position = 'bottom')+
  coord_cartesian(ylim = c(-0.002, 0.0037)) + 
  theme(legend.position = 'bottom',
    text = element_text(size = 16),  # Increase font size
    axis.text = element_text(size = 14),  # Adjust axis text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    strip.text = element_text(size = 14)  # Adjust facet label size
  )
ggsave(paste0("PredAccFULLmodel_TRF_byCond_LOW.pdf"), width = 5, height = 4)
ggsave(paste0("PredAccFULLmodel_TRF_byCond_LOW.png"), width = 5, height = 4)




d <- filter(data, model %in% c("-ioi"))
m.1= lmer(PredAcc~condition  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum'),
          d)
Anova(m.1, type=3)  #

d <- filter(data, model %in% c("-iti"))
m.1= lmer(PredAcc~condition  +(1|sbj) + (1|melID),
          control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 10000)), 
          contrasts = list(condition='contr.sum'),
          d)
Anova(m.1, type=3)  #


#### DISTRUBUTION

ag=with(data, aggregate(cbind(PredAcc) ~sbj+model+condition, FUN="mean"))

ggplot(ag, aes(x = PredAcc, fill = model)) +
  geom_histogram(binwidth = 0.001,alpha = 0.7, position = "identity") +
  facet_wrap(~ condition+model) + # Create separate panels for each condition
  scale_color_manual( values=c("lightblue", "pink"))+
  scale_fill_manual( values=c("lightblue", "pink"))+
  coord_cartesian(xlim = c(-0.02, 0.02)) + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  labs(title = "",x = "Delta r", y = "Frequency",fill = "Condition") +
  theme_minimal()
ggsave(paste0("Distribution_IOIITI_REALSHUFFLED.png"), width = 5, height = 4)
ggsave(paste0("Distribution_IOIITI_REALSHUFFLED.pdf"), width = 5, height = 4)

##real

d = ag[ag$condition=='Real',]
dw <- d %>%
  pivot_wider(
    names_from = model,   # Make columns based on `model`
    values_from = PredAcc # Use `PredAcc` values for the new columns
  )

corr_value <- cor.test(dw$`-ioi`, dw$`-iti`, method = "spearman")
ggplot(dw, aes(x = `-ioi`  , y = `-iti` )) +
  geom_point() +  theme_classic()+
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  
  annotate("text", x = 0, y = 0.01, label = paste("rho = ", round(corr_value$estimate, 2), "p =", 
                                                  round(corr_value$p.value, 2)), size = 5, color = "red") +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.01, 0.01),xlim = c(-0.01, 0.012)) + 
  labs(title = "Correlation IOI  and ITI Delta r (REAL)", x = "Delta IOI", y = "Delta ITI")
ggsave(paste0("Correlation_ioiiti_REAL.pdf"), width = 5, height = 4)
ggsave(paste0("Correlation_ioiiti_REAL.png"), width = 5, height = 4)

##shuff
d = ag[ag$condition=='Shuffled',]
dw <- d %>%
  pivot_wider(
    names_from = model,   # Make columns based on `model`
    values_from = PredAcc # Use `PredAcc` values for the new columns
  )

corr_value <- cor.test(dw$`-ioi`, dw$`-iti`, method = "spearman")
ggplot(dw, aes(x = `-ioi`  , y = `-iti` )) +
  geom_point() +  theme_classic()+
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  
  annotate("text", x = 0, y = 0.01, label = paste("rho = ", round(corr_value$estimate, 2), "p =", 
                                                  round(corr_value$p.value, 2)), size = 5, color = "red") +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.01, 0.01),xlim = c(-0.01, 0.012)) + 
  labs(title = "Correlation ioi  and iti Delta r (SHUFFLED)", x = "Delta ioi", y = "Delta iti")
ggsave(paste0("Correlation_ioiiti_SHUFF.pdf"), width = 5, height = 4)
ggsave(paste0("Correlation_ioiiti_SHUFF.png"), width = 5, height = 4)



####### PLOT MODEL MAIN EFFECT
ag=with(data, aggregate(cbind(PredAcc ) ~sbj+model, FUN="mean"))
ggplot(ag, aes(x=model, y=PredAcc)) +theme_classic()+
  stat_summary(fun=mean, geom="point", size=4)+#coord_cartesian(ylim=c(0,0.002))+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width=0, alpha=0.6)+
  xlab("") + ylab("Delta r (Full-Reduced)")
ggsave(paste0("PredAccFULLmodel_TRF.pdf"), width = 5, height = 4)
ggsave(paste0("PredAccFULLmodel_TRF.png"), width = 5, height = 4)

###### one sample ttest (not corrected)
test_results <- ag %>%
  group_by(model) %>%
  do({
    test <- wilcox.test(.$PredAcc, mu = 0, alternative = "two.side")
    data.frame(
      wilcox_p_value = test$p.value,
      wilcox_statistic = test$statistic
    )
  })


