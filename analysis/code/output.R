# Coptotermes homosexual tandem analysis
# output.R

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
  library(Hmisc)
  library(patchwork)
  library(magick)
  
  require(coxme)
  library(lme4)  
  library(car)
  library(multcomp)
  library(rstatix)
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Survival analysis of tandem duration
#------------------------------------------------------------------------------#
plot_tandem_duration_surv <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_tandem$treatment = factor(df_tandem$treatment, levels=c("FF","MM","FM"))
  ggsurv = ggsurvplot(
    fit      = survfit(Surv(tan_duration, tan_cens) ~ species, 
                       type = "kaplan-meier", data = subset(df_tandem, leader_clarity >= .6)),
    facet.by = "treatment",
    xlim     = c(0,600),
    conf.int = TRUE,
    xlab     = "Duration (sec)", 
    ylab     = "Probability of tandem run",
    censor   = FALSE,
    linetype = "species"
  )
  ggsurv + 
    scale_x_continuous(breaks = seq(0,600,200)) +
    scale_y_continuous(breaks = seq(0,1,1)) +
    scale_color_viridis(discrete = T, option = "D", end=.5, labels=c("C. formosanus", "C. gestroi")) +
    scale_fill_viridis(discrete = T, option = "D", end=.5) +
    scale_linetype(guide="none") +
    facet_wrap(~treatment, 
               nrow = 2,
               ncol = 2,
               as.table = F,
               strip.position = "left",
               labeller = labeller(
                 treatment = c("FM" = "Female-Male", 
                               "FF" = "Female-Female", 
                               "MM" = "Male-Male"))) +
    theme_bw()+
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.position = c(0.7, 0.8),
          #panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face = "italic", size=12),
          text = element_text(size = 12)) +
    guides(fill = "none")
  ggsave("output/plot_tandem_surv.pdf", height = 3, width = 5)
  ggsave("output/plot_tandem_surv.png", height = 3, width = 5)
}
stat_tandem_duration_surv <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_tandem$treatment = factor(df_tandem$treatment, levels=c("FM","FF","MM"))
  
  cat("---------- Comp tandem duration between species for each treat ----------\n")
  for(i in c("FM", "FF", "MM")){
    cat(paste0("\n----- treat:", i, " -----\n\n"))
    m <- coxme(Surv(tan_duration, tan_cens) ~  species + (1|name), 
               data = subset(df_tandem, treatment==i))
    print(m$call)
    cat("\n----- Anova() -----\n\n")
    print(Anova(m))
  }
  cat("-----------------------------------------------------\n\n")
  
  cat("---------- Comp tandem duration among treat for each species ----------\n")
  for(i in c("CF", "CG")){
    cat(paste0("\n----- species:", i, " -----\n\n"))
    m <- coxme(Surv(tan_duration, tan_cens) ~  treatment + (1|name), 
               data = subset(df_tandem, species==i))
    print(m$call)
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
  {
    {
    p <- list()
    name <- names(df_sum)[4:6]
    title_names <- c("Interaction", "Tandem run", "Non tandem interaction")
    for(i in 1:3){
      if(i == 1){
        p[[i]] <- ggplot(df_sum, aes_string(x = "treatment", y = name[i], 
                           fill = "species", col = "species"))+
          geom_boxplot(aes(x=treatment), outlier.shape= NA, 
                       alpha = .75, width = .2, colour = "black",
                       position = position_dodge(0.2)) + 
          geom_point(aes(x = as.numeric(treatment)-0.25, colour=species), 
                     position = position_jitter(width = 0.05),
                     alpha = 0.5,  size=0.5)
      } else {
        p[[i]] <- ggplot(df_sum, aes_string(x = "treatment", y = name[i], 
                                            fill = "species", col = "species"))+
          geom_boxplot(aes(x=treatment), outlier.shape= NA, 
                       alpha = .75, width = .1, colour = "black", lwd = .5,
                       position = position_dodge(0.2)) + 
          geom_point(aes(x = as.numeric(treatment)-0.25, colour=species), 
                     position = position_jitter(width = 0.05),
                     alpha = 0.5,  size=0.1)
      }
      p[[i]] = p[[i]] +
        scale_fill_viridis(discrete=T, option = "D", end = .5, labels=c("C. formosanus", "C. gestroi")) +
        scale_color_viridis(discrete=T, option = "D", end = .5) +
        scale_y_continuous(breaks = seq(0,1800,600)) +
        scale_x_discrete(labels = c("Female-Male", "Female-Female", "Male-Male")) +
        ggtitle(title_names[i]) +
        theme_classic()+
        theme(axis.title.x=element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(face = "italic"),
              text = element_text(size = 9),
              plot.title = element_text(size = 12))+
        guides(col = "none") 
      if(i == 1){
        p[[i]] <- p[[i]] + theme(legend.position  = c(0.8 , 0.9))+
          ylab("Duration (sec)")
      } else {
        p[[i]] <- p[[i]] + 
          theme(legend.position  = "none",
                axis.title.y=element_blank())
        if(i == 2){
          p[[i]] <- p[[i]] + 
            theme(legend.position  = "none",
                  axis.text.x = element_blank())
        }
      }
    }
  }
  {
    image_tan <- image_read("output/termite.png")
    p_ter <- image_ggplot(image_tan) + theme_classic() +
      scale_y_continuous(breaks = c(32,107), labels = c("C. gestroi", "C. formosanus"))+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(face="italic", 
                                       color = c(viridis(3)[2], viridis(3)[1]), 
                                       size=6),
            axis.text.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "inch"))
  }
  
  (free(p[[1]]) / (p_ter) | p[[2]] / p[[3]]) + 
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face="bold", size = 11))
  ggsave("output/plot_tandem_duration.pdf", height = 4, width = 6.5)
  ggsave("output/plot_tandem_duration.png", height = 4, width = 6.5)
  }
}
stat_tandem_duration <- function(){
  load("data_fmt/df_tandem_fmt.rda")
  df_sum$treatment = factor(df_sum$treatment, levels=c("FM","FF","MM"))
  name <- names(df_sum)[4:6]
  
  for(i in 1:3){
    cat(paste("-##########################################################-\n"))
    cat(paste("----------", name[i], "----------\n"))
    y = df_sum[,name[i]] / 1800
    df_sum$logit_prop = log((y+0.01)/(1-y+0.01))
    
    cat("---------- Comp tandem prop between species for each treat ----------\n")
    for(i in c("FM", "FF", "MM")){
      cat(paste0("\n----- treat:", i, ", t.test() -----\n\n"))
      res = t.test(logit_prop ~ species, 
                   data=subset(df_sum, treatment == i))
      print(res)
      cat(paste0("\n----- cohens_d() -----\n\n"))
      res_effect = cohens_d(logit_prop ~ species, var.equal = FALSE, data=subset(df_sum, treatment == i))
      print(res_effect)
    }
    cat("-----------------------------------------------------\n\n")
  
    cat("---------- Comp tandem prop across treats for each species ----------\n")
    for(i in c("CF", "CG")){
      cat(paste0("\n----- treat:", i, ", aov() -----\n\n"))
      res = aov(logit_prop ~ treatment, 
                data=df_sum[df_sum$species==i,])
      print(res)
      cat(paste0("\n----- treat:", i, ", TukeyHSD() -----\n\n"))
      print(TukeyHSD(res))
      
      cat(paste0("\n----- cohens_d() -----\n\n"))
      res_effect = cohens_d(logit_prop ~ treatment, var.equal = FALSE, 
                            data=df_sum[df_sum$species==i,])
      print(res_effect)
    }
    cat("-----------------------------------------------------\n\n")
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Comparison of moving speed
#------------------------------------------------------------------------------#
plot_speed_comparison <- function (){
  load("data_fmt/df_tandem_fmt.rda")
  df_sep$treatment = factor(df_sep$treatment, levels=c("FM","FF","MM"))
  colors <- c("Follower" = viridis(option = "plasma", 3)[1],
              "Leader"   = viridis(option = "plasma", 3)[2])
  # before/after separation
  ggplot(df_sep)+
    stat_summary(aes(x=time, y=speed_follower, color = "Follower"), 
                 fun = 'mean', 
                 geom = 'line') +
    stat_summary(aes(x=time, y=speed_follower, fill = "Follower"), 
                 fun.data = 'mean_cl_normal', 
                 geom = 'ribbon', alpha = 0.2) +
    stat_summary(aes(x=time, y=speed_leader, color = "Leader"), 
                 fun = 'mean', geom = 'line', linetype = "dashed") +
    stat_summary(aes(x=time, y=speed_leader, fill = "Leader"), 
                 fun.data = 'mean_cl_normal', 
                 geom = 'ribbon', alpha = 0.2) +
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
    coord_cartesian(ylim = c(0, 30)) +
    scale_x_continuous(breaks = seq(-5, 5, 5))+
    scale_y_continuous(breaks = seq(0, 30, 5))+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    xlab("Time (sec)")+
    ylab("Speed (mm/sec)") +
    theme_bw()+
    theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = c(0.075,0.6),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        strip.text.y = element_text(face="italic"),
        text = element_text(size = 12),
        legend.text = element_text(size = 7),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm'))+
    guides(color = "none", linetype = "none")
  ggsave("output/plot_move_speed_sep.pdf", height = 4, width = 6)
  ggsave("output/plot_move_speed_sep.png", height = 4, width = 6)
}
stat_speed_comparison <- function() {
  load("data_fmt/df_tandem_fmt.rda")
  df_sep$treatment = factor(df_sep$treatment, levels=c("FM","FF","MM"))
  
  # statistical analysis
  df_sep_ana <- subset(df_sep, time == 5)
  df_sep_ana$sep_event <- paste(df_sep_ana$name, df_sep_ana$sep_event, sep="_")
  df_sep_ana_long <- pivot_longer(df_sep_ana, cols = starts_with("speed"), 
                                  names_to = "role", 
                                  values_to = "speed")
  df_sep_ana_long$role <- gsub("speed_", "", df_sep_ana_long$role)
  
  df_sep_ana$speed_diff <- df_sep_ana$speed_follower - df_sep_ana$speed_leader
  
  cat("---------- Speed comparison Leader VS Follower at 5 sec ----------\n")
  for(i_species in c("CF", "CG")){
    for(i_pair in c("FM", "FF", "MM")){
      cat(paste("###############################\n
                species:", i_species, "pair:", i_pair,"\n
                ################################"))
      r <- t.test(speed ~ role, pared = TRUE,
             data = subset(df_sep_ana_long, species == i_species & treatment == i_pair))
      print(r)
      r <- cohens_d(speed ~ role, 
               data = subset(df_sep_ana_long, species == i_species & treatment == i_pair))
      print(r)
    }
  }
  cat("-----------------------------------------------------\n\n")
  
  
  cat("---------- Speed diff comparison at 5 sec in Cop ges ----------\n")
  r <- lmer(speed_diff ~ treatment + (1|name), 
            data=subset(df_sep_ana, species == "CG"))
  print(summary(r)$call)
  cat("\n----- Anova() -----\n\n")
  print(Anova(r))
  multicomparison<-glht(r, linfct = mcp(treatment = "Tukey"))
  cat("\n----- glht mcp Tukey -----\n\n")
  print(summary(multicomparison))
  cat("-----------------------------------------------------\n\n")
  
  cat("---------- Speed diff comparison at 5 sec in Cop for ----------\n")
  r <- lmer(speed_diff ~ treatment + (1|name), 
            data=subset(df_sep_ana, species == "CF"))
  print(summary(r)$call)
  cat("\n----- Anova() -----\n\n")
  print(Anova(r))
  multicomparison<-glht(r, linfct = mcp(treatment = "Tukey"))
  cat("\n----- glht mcp Tukey -----\n\n")
  print(summary(multicomparison))
  cat("-----------------------------------------------------\n\n")

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Survival analysis of tandem duration across years
#------------------------------------------------------------------------------#
plot_tandem_duration_surv_comp_years <- function(){
  
  {
    load("data_fmt/df_tandem_fmt.rda")
    df_tandem <- subset(df_tandem, treatment == "FM")
    df_tandem_comp <- data.frame(df_tandem, experiment = "Florda_2021")
    df_sum <- subset(df_sum, treatment == "FM")
    df_sum_comp <- data.frame(df_sum, experiment = "Florda_2021")
  
    load("data_fmt/df_tandem_fmt_Mizumoto-etal-2020_JAE.rda")
    df_tandem_comp = rbind(df_tandem_comp,
                           data.frame(df_tandem, experiment = "Florda_2019"))
    df_sum_comp = rbind(df_sum_comp,
                        data.frame(df_sum, experiment = "Florda_2019"))
    
    load("data_fmt/df_tandem_fmt_Mizumoto-etal-2021_PRSB.rda")
    df_tandem_comp = rbind(df_tandem_comp,
                           data.frame(df_tandem, experiment = "Florda_2020"))
    df_sum_comp = rbind(df_sum_comp,
                        data.frame(df_sum, experiment = "Florda_2020"))
  }
  
  ggsurv = ggsurvplot(
    fit      = survfit(Surv(tan_duration, tan_cens) ~ experiment, 
                       type = "kaplan-meier", data = subset(df_tandem_comp, leader_clarity>.5)),
    facet.by = "species",
    xlim     = c(0,600),
    conf.int = TRUE,
    xlab     = "Duration (sec)", 
    ylab     = "Probability of tandem run",
    censor   = FALSE,
    linetype = "experiment"
  )
  p1 <- ggsurv + 
    scale_x_continuous(breaks = seq(0,600,200)) +
    scale_y_continuous(breaks = seq(0,1,1)) +
    scale_color_viridis(option = "plasma", discrete = T,
                        end = .8,
                        labels=c("Florida 2019", 
                                 "Florida 2020", 
                                 "Florida 2021")) +
    scale_fill_viridis(option = "plasma", discrete = T, 
                       end = .8,
                       labels=c("Florida 2019", 
                                "Florida 2020", 
                                "Florida 2021")) +
    scale_linetype_discrete()+
    facet_wrap(~species, 
               nrow = 3,
               ncol = 1,
               strip.position = "left",
               labeller = labeller(
                 species = c("CF" = "C. formosanus",
                             "CG" = "C. gestroi"))) +
    theme_classic()+
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.position = c(0.8, 0.9),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          strip.text.y = element_text(face="italic"),
          text = element_text(size = 10)) +
    guides(fill = "none", linetype = "none")

  p2 <- ggplot(df_sum_comp, aes(x=experiment, y=tandem_total_duration, 
                          fill=species, col=species))+
    geom_boxplot(aes(x=experiment), outlier.shape= NA, 
                 alpha = .75, width = .2, colour = "black") + 
    geom_point(aes(x = (experiment), fill=species), 
               position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = .2),
               alpha = 1, shape = 19, size=0.5)+
    scale_fill_viridis(discrete=T, option = "D", end=.5, labels=c("C. formosanus", "C. gestroi")) +
    scale_color_viridis(discrete=T, option = "D", end=.5,) +
    scale_y_continuous(breaks = seq(0,1800,600)) +
    scale_x_discrete(labels = c("Florida 2019",#\nMizumoto et al. 2020 JAE",
                                "Florida 2020",#\nMizumoto et al. 2021 PRSB",
                                "Florida 2021"))+#\nThis study")) +
    ylab("Tandem duration (sec)") +
    xlab("") +
    theme_classic()+
    theme(legend.position  = "top",
          legend.title = element_blank(),
          legend.text = element_text(face = "italic"),
          text = element_text(size = 10))+
    guides(col = "none")
  
  p1+p2+plot_annotation(tag_levels = 'A')&
    theme(plot.tag.position = c(0.025, 0.99),
          plot.tag = element_text(face="bold", size = 11))
  ggsave("output/plot_year_comparison.pdf", height = 3.5, width = 7)
  ggsave("output/plot_year_comparison.png", height = 3.5, width = 7)
}
stat_tandem_duration_surv_comp <- function(){
  
  {
    load("data_fmt/df_tandem_fmt.rda")
    df_tandem <- subset(df_tandem, treatment == "FM")
    df_tandem_comp <- data.frame(df_tandem, experiment = "Florda_2021")
    df_sum <- subset(df_sum, treatment == "FM")
    df_sum_comp <- data.frame(df_sum, experiment = "Florda_2021")
    
    load("data_fmt/df_tandem_fmt_Mizumoto-etal-2020_JAE.rda")
    df_tandem_comp = rbind(df_tandem_comp,
                           data.frame(df_tandem, experiment = "Florda_2019"))
    df_sum_comp = rbind(df_sum_comp,
                        data.frame(df_sum, experiment = "Florda_2019"))
    
    load("data_fmt/df_tandem_fmt_Mizumoto-etal-2021_PRSB.rda")
    df_tandem_comp = rbind(df_tandem_comp,
                           data.frame(df_tandem, experiment = "Florda_2020"))
    df_sum_comp = rbind(df_sum_comp,
                        data.frame(df_sum, experiment = "Florda_2020"))
  }
  df_tandem_comp$experiment <- factor(df_tandem_comp$experiment)
  df_tandem_comp$name <- paste(df_tandem_comp$name, df_tandem_comp$experiment, sep="_")
  df_sum_comp$experiment <- factor(df_sum_comp$experiment)
  
  ## C. formosanus
  cat("---------- Comp tandem duration across experiments in CF ----------\n")
  m <- coxme(Surv(tan_duration, tan_cens) ~  experiment + (1|name), 
             data = droplevels(subset(df_tandem_comp, species=="CF")))
  print(m$call)
  cat("\n----- Anova() -----\n\n")
  print(Anova(m))
  
  ## C. gestroi
  cat("---------- Comp tandem duration across experiments in CG ----------\n")
  m <- coxme(Surv(tan_duration, tan_cens) ~  experiment + (1|name), 
             data = subset(df_tandem_comp, species=="CG" & leader_clarity>.5))
  print(m$call)
  cat("\n----- Anova() -----\n\n")
  print(Anova(m))
  cat("\n----- glht mcp Tukey --\n\n")
  multicomparison<-glht(m, linfct = mcp(experiment = "Tukey"))
  print(summary(multicomparison))
  
  ## 
  y = df_sum_comp$tandem_total_duration / 1800
  df_sum_comp$logit_tandem_prop = log((y+0.01)/(1-y+0.01))
  
  cat("---------- Comp tandem prop across experiments in CF ----------\n")
  res = t.test(logit_tandem_prop ~ experiment, 
            data=df_sum_comp[df_sum_comp$species=="CF",])
  print(res)
  cat(paste0("\n----- cohens_d() -----\n\n"))
  res_effect = cohens_d(logit_tandem_prop ~ experiment, var.equal = FALSE, 
                        data=droplevels(df_sum_comp[df_sum_comp$species=="CF",]))
  print(res_effect)
  cat("-----------------------------------------------------\n\n")
  
  cat("---------- Comp tandem prop across experiments in CG ----------\n")
  res = aov(logit_tandem_prop ~ experiment, 
               data=df_sum_comp[df_sum_comp$species=="CG",])
  print(Anova(res))
  
  cat(paste0("\n----- cohens_d() -----\n\n"))
  res_effect = cohens_d(logit_tandem_prop ~ experiment, var.equal = FALSE, 
                        data=droplevels(df_sum_comp[df_sum_comp$species=="CG",]))
  print(res_effect)
  cat("-----------------------------------------------------\n\n")
  
  ## Comparison of species
  
  cat("---------- Comparison of species Cox ----------\n")
  m <- coxme(Surv(tan_duration, tan_cens) ~  species + (1|experiment/name), 
             data = df_tandem_comp)
  print(m$call)
  cat("\n----- Anova() -----\n\n")
  print(Anova(m))
  
  cat("---------- Comparison of species LMM ----------\n")
  res = lmer(logit_tandem_prop ~ species + (1|experiment), 
            data=df_sum_comp)
  print(Anova(res))
  
  
}
#------------------------------------------------------------------------------#

      