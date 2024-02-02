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
  {
    library(ggplot2)
    library(survminer)
    library(viridis)
    
    require(coxme)
    library(lme4)  
    library(car)
    library(multcomp)
  }
}
#------------------------------------------------------------------------------#

tandemDurationSurv()
tandemDuration()

#------------------------------------------------------------------------------#
# Survival analysis of tandem duration
#------------------------------------------------------------------------------#
tandemDurationSurv <- function(){
  load("data/rda/df_tandem.rda")
  
  # plots
  {
    df.tan.time$Treat = factor(df.tan.time$Treat, levels=c("FM","FF","MM"))
    theme_set(theme_minimal())
    ggsurv = ggsurvplot(
      fit      = survfit(Surv(Tan.time, Cens) ~ Species + Treat, 
                         type = "kaplan-meier", 
                         data = df.tan.time),
      conf.int = TRUE,
      xlab     = "Duration (sec)", 
      ylab     = "Tandem Prob",
      xlim     = c(0,600),
      legend   = c(0.8,0.2),
      palette  = rep(viridis(2, option = "E"), each=3)
    )
    ggsurv$plot + scale_x_continuous(breaks = seq(0,600,100)) +
      facet_grid(Treat ~ .) 
    pdfName <- paste0("output/TandemDurationSurv.pdf")  
    ggsave(pdfName, height = 5, width = 3)
    
    ggsurv$plot + scale_x_continuous(breaks = seq(0,600,100)) +
      facet_grid(Species ~ .) 
    
  }
  
  # statistics
  {
    sink("output/TandemDurationSurv.txt")
    cat("----------------------------------------------\n")
    cat("Cox mixed effect model, for each combinations\n")
    cat("----------------------------------------------\n")
    
    for(i in c("FM", "FF", "MM")){
      cat(paste0(i, " ----------------------------------------\n"))
      cat("coxme(Surv(Tan.time, Cens) ~ Species + (1|Video))\n")
      m <- coxme(Surv(Tan.time, Cens) ~  Species + (1|Video), 
                 data = df.tan.time[df.tan.time$Treat==i,])
      print(summary(m))
      print(Anova(m))
      cat("\n-----------------------------------------------\n")
    }
    
    cat("----------------------------------------------\n")
    cat("Cox mixed effect model, for each species\n")
    cat("----------------------------------------------\n")
    
    for(i in c("CF", "CG")){
      cat(paste0(i, " ----------------------------------------\n"))
      cat("coxme(Surv(Tan.time, Cens) ~ Treat + (1|Video))\n")
      m <- coxme(Surv(Tan.time, Cens) ~  Treat + (1|Video), 
                 data = df.tan.time[df.tan.time$Species==i,])
      print(summary(m))
      print(Anova(m))
      multicomparison<-glht(m,linfct=mcp(Treat="Tukey"))
      print(summary(multicomparison))
      cat("\n-----------------------------------------------\n")
    }
    sink()
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Comparison of overall tandem period
#------------------------------------------------------------------------------#
tandemDuration <- function(){
  load("data/rda/df_tandem.rda")
  
  # statistics
  {
    df.pair$Treat = factor(df.pair$Treat, levels=c("FM","FF","MM"))
    y = df.pair$Tandem / 1800
    logit.y = log((y+0.01)/(1-y+0.01))
    df.pair$logitTandem = logit.y
    
    fname <- paste0(odir,"TandemDuration.txt")
    sink(fname)
    cat("----------------------------------------------\n")
    cat("t.test between species, for each combinations\n")
    cat("proportion of tandem period was logit transformed\n")
    cat("by adding 0.01 to avoid infinitive\n")
    cat("t.test(logitTandem ~ Species)\n")
    cat("----------------------------------------------\n")
    
    for(i in c("FM", "FF", "MM")){
      cat(paste0(i, " ----------------------------------------\n"))
      res = t.test(logitTandem ~ Species, data=df.pair[df.pair$Treat==i,])
      print(res)
    }
    
    cat("----------------------------------------------\n")
    cat("t.test between combinations, for each species\n")
    cat("proportion of tandem period was logit transformed\n")
    cat("by adding 0.01 to avoid infinitive\n")
    cat("t.test(logitTandem ~ Species)\n")
    cat("----------------------------------------------\n")
    
    for(i in c("CF", "CG")){
      cat(paste0(i, " ----------------------------------------\n"))
      res = aov(logitTandem ~ Treat, data=df.pair[df.pair$Species==i,])
      print(summary(res))
      print(TukeyHSD(res))
    }
    sink()
  }
  
  # plot
  {
    ggplot(df.pair, aes(x=Treat, y=Tandem, fill=Species, col=Species))+
      geom_boxplot(aes(x=Treat), outlier.shape= NA, 
                   alpha = .75, width = .2, colour = "black") + 
      geom_point(aes(x = as.numeric(Treat)-0.2, fill=Species), 
                 position = position_jitter(width = 0.05),
                 alpha = 0.75, shape = 19, size=0.5)+
      scale_fill_viridis(discrete=T, option = "E") +
      scale_color_viridis(discrete=T, option = "E") +
      scale_y_continuous(breaks = seq(0,1800,600)) +
      ylab("Tandem duration (sec)") +
      xlab("") +
      theme(legend.position  = c(0.8 , 0.9))
    pdfName <- paste0(odir, "TandemDuration.pdf")  
    ggsave(pdfName, height = 4, width = 4)
  }
}
#------------------------------------------------------------------------------#
