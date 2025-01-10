
####-------------------------------------
# Project: Bach melodies
##---------------------------------------

library (tidyr)
library(ggplot2)
library(reshape2)
library(pastecs) #for descriptives
library(gridExtra)
library(dplyr)
library(tidyverse)
library(psych)   
library(ggdist) # half violinn plot 

# clear environment
ls()
rm(list=ls())
source("C:/Users/robianco/OneDrive - Fondazione Istituto Italiano Tecnologia/BACHNB/SCRIPTS/functions_RB/fun_plot_geom_flat_violin.R")

####-------------------------------------
setwd("C:/Users/robianco/OneDrive - Fondazione Istituto Italiano Tecnologia/BACHNB/STIMULI/")
####------------------\

## STM IDYOM MODEL 

data=read.table(file="allinfo_stimuli_FINAL_StimMatrix_both+.csv",header=F, sep=',')
colnames(data) <- c("melid", "pitch", "onset", "iti","ioi", "so", "eo", "sp", 'ep') 
data$cond = ifelse(data$melid<11, 1, 2)
data$melid =as.factor(data$melid)
data$melname = fct_recode(data$melid, "s_01" = "11","s_05" = "12","s_08" = "13","s_10" = "14",
                        "o_05"= "5", "o_08"= "8", "o_03"= "3", "o_06"= "6", "o_04"= "4",
                        "o_09"= "9", "o_02"= "2", "o_01"= "1", "o_07"= "7", "o_10"= "10")
data$cond = as.factor(data$cond)
data$cond = fct_recode(data$cond, "Real" = "1","Shuffled" = "2")
data <- data[data$iti != 99, ] #remove first notes of songs



mean_value <- mean(data$iti)
range_value <- range(data$iti)
sd_value <- sd(data$iti)

mean_value <- mean(data$ioi)
range_value <- range(data$ioi)
sd_value <- sd(data$ioi)

########TEST DIFFERENCE BETWEEN CONDITIONS
df_long <- data %>%
  pivot_longer(cols = c(so, sp, eo, ep),
               names_to = "predictor",
               values_to = "value")

ggplot(df_long, aes(predictor,value, fill=cond, color =cond)) + theme_minimal()+
  geom_flat_violin(position = position_nudge(x = 0, y = 0)) +
 # coord_flip()+  
  stat_summary(fun = mean, geom = "point", position = position_nudge(x = -0.1), size = 3) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.05, position = position_nudge(x = -0.1)) +
  scale_fill_manual(values = scales::alpha(c("red", "grey"), alpha = 0.7)) + # Add transparency here
  scale_color_manual(values = c("red", "grey")) +
  xlab("") + ylab(" ")+theme(legend.position="bottom")

#######AGGREGATE BY MELODY
ag=with(data, aggregate(cbind(so, sp, eo, ep) ~melid+cond, FUN="mean")) #aggregate by melody
df_long <- ag %>%
  pivot_longer(cols = c(so, sp, eo, ep), 
               names_to = "predictor", 
               values_to = "value")
desired_order <- c( "so", "eo", "sp", "ep")
df_long <- df_long %>%
  mutate(predictor = factor(predictor, levels = desired_order))

ggplot(df_long, aes(predictor,value, fill=cond, color =cond)) + theme_minimal()+
  geom_flat_violin(position = position_nudge(x = 0, y = 0)) +
  # coord_flip()+  
  stat_summary(fun = mean, geom = "point", position = position_nudge(x = -0.1), size = 3) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.05, position = position_nudge(x = -0.1)) +
  scale_fill_manual(values = scales::alpha(c("red", "grey"), alpha = 0.7)) + # Add transparency here
  scale_color_manual(values = c("red", "grey")) +
  xlab("") + ylab(" ")+theme(legend.position="bottom")+  ggtitle("Grouped by melody, 1 SD")
ggsave(paste0("SurpEntr_bycondition.pdf"), width = 9, height = 4)
ggsave(paste0("SurpEntr_bycondition.png"), width = 9, height = 4)


###plot stim
level_mel =c("s_05", "s_08","s_01",  "s_10", "o_05", "o_08", "o_03","o_04", "o_09", "o_06",  "o_02",  "o_07", "o_01","o_10")
ag=with(data, aggregate(cbind(so, sp,  eo, ep) ~melid+melname+cond, FUN="mean")) #aggregate by melody
ag <- ag %>% 
  arrange(factor(melname, levels = level_mel))
ag$cond=fct_relevel(ag$cond, "Shuffled", "Real")
df_long <- ag %>%
  pivot_longer(cols = c(so, sp,  eo, ep), 
               names_to = "predictor", 
               values_to = "value")
df_long$info <- ifelse(df_long$predictor %in% c("so", "sp"), 'Surprise', 'Entropy')
df_long$dimension <- ifelse(df_long$predictor %in% c("so", "eo"), 'Timing', 'Pitch')
df_long$info=fct_relevel(df_long$info, "Surprise", "Entropy")
df_long$dimension=fct_relevel(df_long$dimension, "Timing", "Pitch")

ggplot(df_long, aes(x =factor(melname, levels = level_mel), y=value, 
                    group = interaction(info, dimension),
                    color = dimension, 
                    fill =  dimension, 
                    lty= info)) + theme_classic()+ facet_wrap(~cond)+
  scale_color_manual( values=c("forestgreen", "goldenrod1"))+
  scale_fill_manual( values=c("forestgreen", "goldenrod1"))+
  stat_summary(fun=mean, geom="point", size = 4) +
  stat_summary(fun=mean, geom="line") +
  xlab("") + ylab("Value") +
  theme( text = element_text(size = 16),  # Increase font size
        axis.text = element_text(size = 14),  # Adjust axis text size
        axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x-axis labels
        axis.title = element_text(size = 16),  # Adjust axis title size
        strip.text = element_text(size = 14)  # Adjust facet label size
  )
ggsave(paste0("S&E_perMelody.pdf"), width = 8, height = 4)
ggsave(paste0("S&E_perMelody.png"), width = 6, height = 6)


ggplot(df_long, aes(x =interaction(info, dimension), y=value, 
                    #group = interaction(info, dimension),
                    color = dimension, 
                    fill =  dimension, 
                    lty= info)) + theme_classic()+facet_wrap(~cond)+
  scale_color_manual( values=c("forestgreen", "goldenrod1"))+
  scale_fill_manual( values=c("forestgreen", "goldenrod1"))+
  stat_summary(fun=mean, geom="point", size = 4) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0) +
  xlab("") + ylab("Value") +
  theme( text = element_text(size = 16),  # Increase font size
         axis.text = element_text(size = 14),  # Adjust axis text size
         axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x-axis labels
         axis.title = element_text(size = 16),  # Adjust axis title size
         strip.text = element_text(size = 14)  # Adjust facet label size
  )
ggsave(paste0("S&E_mean.pdf"), width = 4, height = 4)
ggsave(paste0("S&E_mean.png"), width = 6, height = 6)


# Perform normality test (Shapiro-Wilk test) for each predictor
df_long = as.data.frame(df_long)
g = df_long[df_long$cond == 'Real' & df_long$predictor == 'sp',]
shapiro.test(g$value)
g = df_long[df_long$cond == 'Real' & df_long$predictor == 'so',]
shapiro.test(g$value)
g = df_long[df_long$cond == 'Shuffled' & df_long$predictor == 'ep',]
shapiro.test(g$value)
g = df_long[df_long$cond == 'Shuffled' & df_long$predictor == 'eo',]
shapiro.test(g$value)


########PLOT CORRELATIONS
df <- data[, c(4:9)]

library(corrplot)
cor_results <- corr.test(df)
cor_matrix <- cor_results$r
p_mat <- cor_results$p
pdf("FIG_corr_both+_final.pdf", width = 10, height = 10)  # Open PNG device
# Create a basic plot to initialize the device (in some environments)
# corrplot(cor_matrix, 
#          method = "circle", 
#          type = "upper")
corrplot(cor_matrix, 
         method = "circle", 
         type = "upper", 
         # add = TRUE, 
         # tl.col = "black", 
         tl.srt = 45,  # Text label color and angle
         addCoef.col = "black",  # Color of the correlation coefficients
         number.cex = 0.7,  # Size of the coefficient text
         diag = T)  # Don't show diagonal
dev.off()

