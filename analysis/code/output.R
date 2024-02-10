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
  library(tidyr)
  
  library(ggplot2)
  library(survminer)
  library(viridis)
  
  require(coxme)
  library(lme4)  
  library(car)
  library(multcomp)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Survival analysis of tandem duration
#------------------------------------------------------------------------------#
plot_tandem_duration_surv <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_tandem$treatment = factor(df_tandem$treatment, levels=c("FM","FF","MM"))
  ggsurv = ggsurvplot(
    fit      = survfit(Surv(tan_duration, tan_cens) ~ species, 
                       type = "kaplan-meier", data = df_tandem),
    facet.by = "treatment",
    xlim     = c(0,600),
    conf.int = TRUE,
    xlab     = "Duration (sec)", 
    ylab     = "Probability of tandem run",
    censor   = FALSE
  )
  ggsurv + 
    scale_x_continuous(breaks = seq(0,600,200)) +
    scale_y_continuous(breaks = seq(0,1,1)) +
    scale_color_viridis(discrete = T, option = "E", labels=c("C. formosanus", "C. gestroi")) +
    scale_fill_viridis(discrete = T, option = "E") +
    facet_wrap(~treatment, 
               nrow = 3,
               ncol = 1,
               strip.position = "left") +
    theme_classic()+
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.position = c(0.8, 0.9),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face = "italic"),
          text = element_text(size = 12)) +
    guides(fill = "none")
  ggsave("output/plot_tandem_surv.pdf", height = 5, width = 3)
}
stat_tandem_duration_surv <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_tandem$treatment = factor(df_tandem$treatment, levels=c("FM","FF","MM"))
  
  cat("---------- Comp tandem duration between species for each treat ----------\n")
  for(i in c("FM", "FF", "MM")){
    cat(paste0("\n----- treat:", i, ", summary() -----\n\n"))
    m <- coxme(Surv(tan_duration, tan_cens) ~  species + (1|name), 
               data = subset(df_tandem, treatment==i))
    print(summary(m))
    cat("\n----- Anova() -----\n\n")
    print(Anova(m))
  }
  cat("-----------------------------------------------------\n\n")
  
  cat("---------- Comp tandem duration among treat for each species ----------\n")
  for(i in c("CF", "CG")){
    cat(paste0("\n----- species:", i, ", summary() -----\n\n"))
    m <- coxme(Surv(tan_duration, tan_cens) ~  treatment + (1|name), 
               data = subset(df_tandem, species==i))
    print(summary(m))
    cat("\n----- Anova() -----\n\n")
    print(Anova(m))
    cat("\n----- glht mcp Tukey --\n\n")
    multicomparison<-glht(m, linfct = mcp(treatment = "Tukey"))
    print(summary(multicomparison))
  }
  cat("-----------------------------------------------------\n\n")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Comparison of overall tandem period
#------------------------------------------------------------------------------#
plot_tandem_duration <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_sum$treatment = factor(df_sum$treatment, levels=c("FM","FF","MM"))
  ggplot(df_sum, aes(x=treatment, y=tandem_total_duration, 
                     fill=species, col=species))+
    geom_boxplot(aes(x=treatment), outlier.shape= NA, 
                 alpha = .75, width = .2, colour = "black") + 
    geom_point(aes(x = as.numeric(treatment)-0.2, fill=species), 
               position = position_jitter(width = 0.05),
               alpha = 0.75, shape = 19, size=0.5)+
    scale_fill_viridis(discrete=T, option = "E", labels=c("C. formosanus", "C. gestroi")) +
    scale_color_viridis(discrete=T, option = "E") +
    scale_y_continuous(breaks = seq(0,1800,600)) +
    scale_x_discrete(labels = c("Female-Male", "Female-Female", "Male-Male")) +
    ylab("Tandem duration (sec)") +
    xlab("") +
    theme_classic()+
    theme(legend.position  = c(0.8 , 0.9),
          legend.title = element_blank(),
          legend.text = element_text(face = "italic"),
          text = element_text(size = 12))+
    guides(col = "none")
  ggsave("output/plot_tandem_duration.pdf", height = 4, width = 4)
}
stat_tandem_duration <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_sum$treatment = factor(df_sum$treatment, levels=c("FM","FF","MM"))
  y = df_sum$tandem_total_duration / 1800
  df_sum$logit_tandem_prop = log((y+0.01)/(1-y+0.01))
  
  cat("---------- Comp tandem prop between species for each treat ----------\n")
  for(i in c("FM", "FF", "MM")){
    cat(paste0("\n----- treat:", i, ", t.test() -----\n\n"))
    res = t.test(logit_tandem_prop ~ species, 
                 data=subset(df_sum, treatment == i))
    print(res)
  }
  cat("-----------------------------------------------------\n\n")

  cat("---------- Comp tandem prop across treats for each species ----------\n")
  for(i in c("CF", "CG")){
    cat(paste0("\n----- treat:", i, ", aov() -----\n\n"))
    res = aov(logit_tandem_prop ~ treatment, 
              data=df_sum[df_sum$species==i,])
    print(res)
    cat(paste0("\n----- treat:", i, ", TukeyHSD() -----\n\n"))
    print(TukeyHSD(res))
  }
  cat("-----------------------------------------------------\n\n")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Comparison of moving speed
#------------------------------------------------------------------------------#
plot_speed_comparison <- function (){
  load("data_fmt/df_tandem_fmt.rda")
  df_sep$treatment = factor(df_sep$treatment, levels=c("FM","FF","MM"))
  # before/after separation
  ggplot(df_sep)+
    stat_smooth(aes(x=time, y=speed0), col=viridis(2)[1], se = FALSE)+
    stat_smooth(aes(x=time, y=speed1), col=viridis(2)[2], se = FALSE)+
    geom_vline(xintercept = 0, linetype = 2, colour="#222222") +
    facet_grid(rows = vars(species),
               cols = vars(treatment),
               switch = "y",
               labeller = labeller(
                 species = c("CF" = "C. formosanus",
                             "CG" = "C. gestroi"),
                 treatment = c("FM" = "Female-Male", 
                               "FF" = "Female-Female", 
                               "MM" = "Male-Male"))
               )+
    coord_cartesian(ylim = c(0, 25)) +
    scale_x_continuous(breaks = seq(-5, 5, 5))+
    scale_y_continuous(breaks = seq(0, 25, 5))+
    xlab("Time (sec)")+
    ylab("Speed (mm/sec)") +
    theme_bw()+
    theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = c(0.8, 0.9),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        strip.text.y = element_text(face="italic"),
        text = element_text(size = 12))
  ggsave("output/plot_move_speed_sep.pdf", height = 4, width = 6)
}
stat_speed_comparison <- function() {
  load("data_fmt/df_tandem_fmt.rda")
  df_sep$treatment = factor(df_sep$treatment, levels=c("FM","FF","MM"))
  
  # statistical analysis
  df_sep_ana <- subset(df_sep, time == 5)
  df_sep_ana$speed_diff <- df_sep_ana$speed0 - df_sep_ana$speed1
  df_sep_ana$sep_event <- paste(df_sep_ana$name, df_sep_ana$sep_event, sep="_")
  
  cat("---------- Speed diff comparison at 5 sec in Cop ges ----------\n")
  r <- lmer(speed_diff ~ treatment + (1|name), 
            data=subset(df_sep_ana, species == "CG"))
  cat("\n----- summary() -----\n\n")
  print(summary(r))
  cat("\n----- Anova() -----\n\n")
  print(Anova(r))
  multicomparison<-glht(r, linfct = mcp(treatment = "Tukey"))
  cat("\n----- glht mcp Tukey -----\n\n")
  print(summary(multicomparison))
  cat("-----------------------------------------------------\n\n")
  
  cat("---------- Speed diff comparison at 5 sec in Cop for ----------\n")
  r <- lmer(speed_diff ~ treatment + (1|name), 
            data=subset(df_sep_ana, species == "CF"))
  cat("\n----- summary() -----\n\n")
  print(summary(r))
  cat("\n----- Anova() -----\n\n")
  print(Anova(r))
  multicomparison<-glht(r, linfct = mcp(treatment = "Tukey"))
  cat("\n----- glht mcp Tukey -----\n\n")
  print(summary(multicomparison))
  cat("-----------------------------------------------------\n\n")

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Survival analysis of tandem duration across studies
#------------------------------------------------------------------------------#
plot_tandem_duration_surv_comp <- function(){
  
  load("data_fmt/df_tandem_fmt.rda")
  df_tandem <- subset(df_tandem, treatment == "FM")
  df_tandem_comp <- data.frame(df_tandem, experiment = "Florda_2021")

  load("data_fmt/df_tandem_fmt_Mizumoto-etal-2020.rda")
  df_tandem_comp = rbind(df_tandem_comp,
                         data.frame(df_tandem, experiment = "Florda_2020"))
  
  ggsurv = ggsurvplot(
    fit      = survfit(Surv(tan_duration, tan_cens) ~ experiment, 
                       type = "kaplan-meier", data = df_tandem_comp),
    facet.by = "species",
    xlim     = c(0,600),
    conf.int = TRUE,
    xlab     = "Duration (sec)", 
    ylab     = "Probability of tandem run",
    censor   = FALSE
  )
  ggsurv + 
    scale_x_continuous(breaks = seq(0,600,200)) +
    scale_y_continuous(breaks = seq(0,1,1)) +
    scale_color_viridis(discrete = T, option = "E") +
    scale_fill_viridis(discrete = T, option = "E") +
    facet_wrap(~species, 
               nrow = 3,
               ncol = 1,
               strip.position = "left") +
    theme_classic()+
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.position = c(0.8, 0.9),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face = "italic"),
          text = element_text(size = 12)) +
    guides(fill = "none")
  ggsave("output/plot_tandem_surv.pdf", height = 5, width = 3)
}
stat_tandem_duration_surv <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_tandem <- subset(df_tandem, treatment == "FM")
  df_tandem_comp <- data.frame(df_tandem, experiment = "Florda_2021")
  
  load("data_fmt/df_tandem_fmt_Mizumoto-etal-2020.rda")
  df_tandem_comp = rbind(df_tandem_comp,
                         data.frame(df_tandem, experiment = "Florda_2020"))
  
  m <- coxme(Surv(tan_duration, tan_cens) ~  experiment + (1|name), 
             data = subset(df_tandem_comp, species=="CF"))
  print(summary(m))
  cat("\n----- Anova() -----\n\n")
  print(Anova(m))
  cat("\n----- glht mcp Tukey --\n\n")
  multicomparison<-glht(m, linfct = mcp(treatment = "Tukey"))
  print(summary(multicomparison))
  
  
  cat("---------- Comp tandem duration between species for each treat ----------\n")
  for(i in c("FM", "FF", "MM")){
    cat(paste0("\n----- treat:", i, ", summary() -----\n\n"))
    m <- coxme(Surv(tan_duration, tan_cens) ~  species + (1|name), 
               data = subset(df_tandem, treatment==i))
    print(summary(m))
    cat("\n----- Anova() -----\n\n")
    print(Anova(m))
  }
  cat("-----------------------------------------------------\n\n")
  
  cat("---------- Comp tandem duration among treat for each species ----------\n")
  for(i in c("CF", "CG")){
    cat(paste0("\n----- species:", i, ", summary() -----\n\n"))
    m <- coxme(Surv(tan_duration, tan_cens) ~  treatment + (1|name), 
               data = subset(df_tandem, species==i))
    print(summary(m))
    cat("\n----- Anova() -----\n\n")
    print(Anova(m))
    cat("\n----- glht mcp Tukey --\n\n")
    multicomparison<-glht(m, linfct = mcp(treatment = "Tukey"))
    print(summary(multicomparison))
  }
  cat("-----------------------------------------------------\n\n")
}
#------------------------------------------------------------------------------#

