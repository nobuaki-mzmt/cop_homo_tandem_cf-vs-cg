## Data analysis for Coptotermes homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file plots trajectories to check tracking errors
#------------------------------------------------------------------------------#

rm(list = ls())
scale_trajectory()

#------------------------------------------------------------------------------#
{
  library(data.table)
  library(stringr)
  
  library(ggplot2)
  library(viridis)
  
  ## parameters
  error_speed_thersh = 20
  fps_original <- 30
  fps_analysis <- 5
  dish_mm <- 150
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# reads all csv files and 1) output raw trajectories in .png file
# 2) convert them in 5FPS and R data.frame saved in .rda files
#------------------------------------------------------------------------------#
scale_trajectory <- function(Plot = T, Dataframe = T){
  options(warn = 0)
  
  # data
  df_meta <- data.frame(fread("data_raw/bodylength/scale_bodylength.csv", header=T))
  
  rawdata <-  list.files("data_raw/trajectories", full.names = TRUE, pattern = ".csv")
  dataname <- list.files("data_raw/trajectories", full.names = FALSE, pattern = ".csv")
  
  df_tra <- data.frame()
  df_ind  <- data.frame()
  for(v in 1:length(rawdata)){
    
    # file info
    d <- data.frame(fread(rawdata[v], header=T))
    d <- d[1:(30*60*30+1),]
    
    species = substr(dataname[v], 1, 2)
    treat = substr(dataname[v], 4, 5)
    id = substr(dataname[v], 6, 7)
    name = paste(species, treat, id, sep="_")
    print(paste(v, "/", length(rawdata), "->", name))

    # scale
    d[,1] <- d[,1]/fps_original
    d <- d[seq(1,54000,fps_original/fps_analysis),]
    
    
    scale_factor = dish_mm/df_meta[df_meta$name == name,"scale"]
    
    body_length <- df_meta[df_meta$name == name, c("bodyLength0", "bodyLength1")] * scale_factor
      
    d[,2:5] <- round(d[,2:5] * scale_factor, 4)
    colnames(d)[1] <- "time"
    
    if(Dataframe){
      df_tra <- rbind(df_tra, data.frame(name, d))
      df_ind <- rbind(df_ind, data.frame(name, body_length))
    }
    
    # plot
    if(Plot){
      ggplot(d) + 
        geom_path(aes(x=x0,y=y0), col = viridis(2)[1])+
        geom_path(aes(x=x1,y=y1), col = viridis(2)[2])+
        theme_bw() +
        theme(aspect.ratio = 1) +
        ggtitle(paste(v, name))
      ggsave(paste0("output/trajectory/", name, ".png"), width = 4, height = 4)
    }
    
    # speed error check
    error_speed_thersh = 20
    error_check <- abs(diff(d$x0)) > 15 | abs(diff(d$x1)) > 15 | abs(diff(d$y0)) > 15 | abs(diff(d$y1)) > 15
    if( sum(error_check)>0) {
      print(d[error_check,])
    }
  }  
  save(df_tra, df_ind, file = "data_fmt/data_trajectory.rda", compress="xz")
}
#------------------------------------------------------------------------------#
