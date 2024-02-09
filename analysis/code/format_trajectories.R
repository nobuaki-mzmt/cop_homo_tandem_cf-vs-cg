## Data analysis for Coptotermes homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file plots trajectories to check tracking errors
#------------------------------------------------------------------------------#

rm(list = ls())
scale_trajectory(Plot=F)
tandem.detect()

#------------------------------------------------------------------------------#
{
  library(data.table)
  library(stringr)
  library(dplyr)
  
  library(ggplot2)
  library(viridis)
  
  ## parameters
  error_speed_thersh = 20
  fps_original <- 30
  fps_analysis <- 5
  dish_mm <- 150
  
  min.sep.sec <- 2
  min.tandem.sec <- 5
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
  
  df_tra <- tibble()
  df_ind  <- tibble()
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
      df_tra <- bind_rows(df_tra, tibble(name, d))
      df_ind <- bind_rows(df_ind, tibble(name, body_length))
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
    error_speed_thersh = 30
    error_check <- abs(diff(d$x0)) > 15 | abs(diff(d$x1)) > 15 | abs(diff(d$y0)) > 15 | abs(diff(d$y1)) > 15
    if( sum(error_check)>0) {
      print(d[error_check,])
    }
  }  
  save(df_tra, df_ind, file = "data_fmt/data_trajectory.rda", compress="xz")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads tandem vector (0/1), and
# smooth it so that tandem shorter than min.sec is treated as non tandem
#------------------------------------------------------------------------------#
tandem.smoothing <- function(vec, min.sec){ 
  if(sum(vec)>0){
    timing <- which(vec)[c(T, diff(which(vec))>1)]
    end    <- which(vec)[c(diff(which(vec))>1,T)]
    for(fi in 1:length(timing)){
      if(length( vec[timing[fi]:end[fi]]) < min.sec ){
        vec[timing[fi]:end[fi]] <- F
      }
    }
  }
  return(vec)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# This function reads allTrajectoryData.rda files, then
# 1) Determine if a pair is doing tandem or not for each frame.
# 2) save data.frame for each tandem event (),
#    and for summary for each pair (df_sum_tandem.rda)
#------------------------------------------------------------------------------#
tandem.detect <- function(){
  
  load("data_fmt/data_trajectory.rda")
  
  name_list <- df_ind$name
  
  df_tandem = df_sep = df_sum <- NULL
  
  for(i in 1:length(name_list)){
    print(paste(i, "/", length(name_list), "->", name_list[i]))
    
    df_temp <- subset(df_tra, name==name_list[i])
    df_temp$speed0 <- c(NA, sqrt(diff(df_temp$x0)^2 + diff(df_temp$y0)^2))*fps_analysis
    df_temp$speed1 <- c(NA, sqrt(diff(df_temp$x1)^2 + diff(df_temp$y1)^2))*fps_analysis
    sumbodylength <- sum(as.numeric(subset(df_ind, name==name_list[i])[,2:3]))
    ind_dis = sqrt( (df_temp$x0 - df_temp$x1)^2 + (df_temp$y0 - df_temp$y1)^2 )
    
    interaction <- ind_dis < sumbodylength
    non.interaction <- !interaction
    
    # remove short separation
    tandem <- !tandem.smoothing(non.interaction, min.sep.sec*5) 
    # remove short tandem
    tandem <- tandem.smoothing(tandem, min.tandem.sec*5) 
    
    if(sum(tandem)>0){
      tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
      tan.sta <- which(tandem)[c(T, diff(which(tandem))>1)]
    }
    
    tan_duration <- tan.end - tan.sta
    sep_duration <- tan.sta[-1] - tan.end[-length(tan.end)]
    
    tan_cens <- rep(1, length(tan_duration))
    tan_cens[tan.sta == 1 | tan.end == 9000] <- 0
    
    speed_tandem = (df_temp$speed0[tandem]+df_temp$speed1[tandem])/2
    df_sum_temp <- data.frame(
      name =  name_list[i],
      species = str_split(name_list[i], "_")[[1]][1],
      treatment = str_split(name_list[i], "_")[[1]][2],
      tandem_total_duration = sum(tandem) / fps_analysis,
      speed_tandem = mean(speed_tandem, na.rm=T)
    )

    df_tandem_temp <- data.frame(
      name =  name_list[i],
      species = str_split(name_list[i], "_")[[1]][1],
      treatment = str_split(name_list[i], "_")[[1]][2],
      tan_duration = tan_duration/fps_analysis,
      tan_cens
    )
    
    df_sep_temp <- NULL
    if(length(sep_duration)>0){
      sep_thresh = 25
      for(i_sep in 1:length(sep_duration)){
        if(tan.end[i_sep]>(min.tandem.sec*fps_analysis)+1){
          if(sep_duration[i_sep] > sep_thresh){
            sep_range <- -(min.tandem.sec*fps_analysis):sep_thresh+tan.end[i_sep]
            speed0 = (df_temp$speed0[sep_range])
            speed1 = (df_temp$speed1[sep_range])
            if(mean(speed0) < mean(speed1)){
              temp = speed0 
              speed0 = speed1
              speed1 = temp
            }
            df_sep_temp_temp <- data.frame(
              name =  name_list[i],
              species = str_split(name_list[i], "_")[[1]][1],
              treatment = str_split(name_list[i], "_")[[1]][2],
              time = (-(min.tandem.sec*fps_analysis):sep_thresh)/fps_analysis,
              sep_event = i_sep,
              speed0, speed1
            )
            df_sep_temp <- rbind(df_sep_temp, df_sep_temp_temp)
          }
        }
      }
    }
    
    df_sum <- rbind(df_sum, df_sum_temp)
    df_tandem <- rbind(df_tandem, df_tandem_temp)
    df_sep <- rbind(df_sep, df_sep_temp)
  }  
  
  save(df_sum, df_tandem, df_sep, file = "data_fmt/df_tandem_fmt.rda")
  
}
#------------------------------------------------------------------------------#
