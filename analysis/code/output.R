# Coptotermes homosexual tandem analysis
# output.R
# N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for producing all results
#------------------------------------------------------------------------------#

rm(list = ls())


#------------------------------------------------------------------------------#
{
  # packages
  library(ggplot2)
  library(survminer)
  library(viridis)
  
  require(coxme)
  library(lme4)  
  library(car)
  library(multcomp)
}
#------------------------------------------------------------------------------#

output_tandem_duration_durv()
tandemDuration()

#------------------------------------------------------------------------------#
# Survival analysis of tandem duration
#------------------------------------------------------------------------------#
output_tandem_duration_durv <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  
  # plots
  df_tandem$treatment = factor(df_tandem$treatment, levels=c("FM","FF","MM"))
  ggsurv = ggsurvplot(
    fit      = survfit(Surv(tan_duration, tan_cens) ~ species + treatment, 
                       type = "kaplan-meier", data = df_tandem),
    conf.int = TRUE,
    xlab     = "Duration (sec)", 
    ylab     = "Tandem Prob",
    xlim     = c(0,600),
    legend   = c(0.8,0.2),
    palette  = rep(viridis(2, option = "E"), each=3),
    censor   = FALSE
  )
  ggsurv$plot + scale_x_continuous(breaks = seq(0,600,100)) +
    facet_grid(treatment ~ .) +
    theme(legend.position = "none")
  ggsave("output/plot_tandem_surv.pdf", height = 5, width = 3)

  # statistics
  {
    for(i in c("FM", "FF", "MM")){
      m <- coxme(Surv(tan_duration, tan_cens) ~  species + (1|name), 
                 data = subset(df_tandem, treatment==i))
      print(summary(m))
      print(Anova(m))
    }
    
    for(i in c("CF", "CG")){
      m <- coxme(Surv(tan_duration, tan_cens) ~  treatment + (1|name), 
                 data = subset(df_tandem, species==i))
      print(summary(m))
      print(Anova(m))
      multicomparison<-glht(m, linfct = mcp(treatment = "Tukey"))
      print(summary(multicomparison))
    }
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Comparison of overall tandem period
#------------------------------------------------------------------------------#
output_tandem_duration <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_sum$treatment = factor(df_sum$treatment, levels=c("FM","FF","MM"))
  ggplot(df_sum, aes(x=treatment, y=tandem_total_duration, 
                     fill=species, col=species))+
    geom_boxplot(aes(x=treatment), outlier.shape= NA, 
                 alpha = .75, width = .2, colour = "black") + 
    geom_point(aes(x = as.numeric(treatment)-0.2, fill=species), 
               position = position_jitter(width = 0.05),
               alpha = 0.75, shape = 19, size=0.5)+
    scale_fill_viridis(discrete=T, option = "E") +
    scale_color_viridis(discrete=T, option = "E") +
    scale_y_continuous(breaks = seq(0,1800,600)) +
    ylab("Tandem duration (sec)") +
    xlab("") +
    theme(legend.position  = c(0.8 , 0.9))+
    theme_classic()
  ggsave("output/plot_tandem_duration.pdf", height = 4, width = 4)
  
  # statistics
  y = df_sum$tandem_total_duration / 1800
  df_sum$logit_tandem_prop = log((y+0.01)/(1-y+0.01))

  for(i in c("FM", "FF", "MM")){
    res = t.test(logit_tandem_prop ~ species, 
                 data=subset(df_sum, treatment == i))
    print(res)
  }
  
  cat("----------------------------------------------\n")
  cat("t.test between combinations, for each species\n")
  cat("proportion of tandem period was logit transformed\n")
  cat("by adding 0.01 to avoid infinitive\n")
  cat("t.test(logitTandem ~ Species)\n")
  cat("----------------------------------------------\n")
  
  for(i in c("CF", "CG")){
    res = aov(logit_tandem_prop ~ treatment, 
              data=df_sum[df_sum$species==i,])
    print(summary(res))
    print(TukeyHSD(res))
  }
}
#------------------------------------------------------------------------------#


load("data_fmt/df_tandem_fmt.rda")
ggplot(df_sum, aes(x=species, y=speed_tandem, fill=treatment))+
  geom_boxplot(alpha=.3)

ggplot(df_sep)+
  stat_smooth(aes(x=time, y=speed0), col=viridis(2)[1], alpha=.5)+
  stat_smooth(aes(x=time, y=speed1), col=viridis(2)[2], alpha=.5)+
  facet_grid(species~treatment)
