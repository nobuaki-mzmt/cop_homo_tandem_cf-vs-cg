##
## Data analysis for Coptotermes homosexual tandem analysis
## 

#------------------------------------------------------------------------------#
# This file plots trajectories to check tracking errors
#------------------------------------------------------------------------------#

rm(list = ls())
scale_trajectory(Plot=F)
scale_trajectory(dataset = "Mizumoto-etal-2021_PRSB", dish_mm = 138)
scale_trajectory(dataset = "Mizumoto-etal-2020_JAE", dish_mm = 138)
tandem.detect()
tandem.detect(dataset = "Mizumoto-etal-2021_PRSB")
tandem.detect(dataset = "Mizumoto-etal-2020_JAE")
#------------------------------------------------------------------------------#

{
  library(data.table)
  library(stringr)
  library(dplyr)
  
  library(ggplot2)
  library(viridis)
  
  ## parameters
  
  min.sep.sec <- 2
  min.tandem.sec <- 5
  tandem_leader_thresh <- .5
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# reads all csv files and 1) output raw trajectories in .png file
# 2) convert them in 5FPS and R data.frame saved in .rda files
#------------------------------------------------------------------------------#
scale_trajectory <- function(Plot = T, Dataframe = T, dataset = NULL,
                             error_speed_thersh = 20, fps_original = 30,
                             fps_analysis = 5, dish_mm = 150){
  options(warn = 0)
  
  # data
  data_loc_met = "data_raw/bodylength/scale_bodylength.csv"
  data_loc_tra = "data_raw/trajectories"
  plot_loc_tra = "output/trajectory/"
  data_loc_fmt = paste0("data_fmt/data_trajectory.rda")
  if(length(dataset) > 0){
    data_loc_met = paste0("data_raw/bodylength_", dataset, "/scale_bodylength.csv")
    data_loc_tra = paste0("data_raw/trajectories_", dataset)
    plot_loc_tra = paste0("output/trajectory_", dataset, "/")
    data_loc_fmt = paste0("data_fmt/data_trajectory_", dataset, ".rda")
  }
  df_meta <- data.frame(fread(data_loc_met, header=T))
  
  rawdata <-  list.files(data_loc_tra, full.names = TRUE, pattern = ".csv")
  dataname <- list.files(data_loc_tra, full.names = FALSE, pattern = ".csv")
  
  df_tra <- tibble()
  df_ind  <- tibble()
  for(v in 1:length(rawdata)){
    
    # file info
    d <- data.frame(fread(rawdata[v], header=T))
    
    species = substr(dataname[v], 1, 2)
    treat = substr(dataname[v], 4, 5)
    id = substr(dataname[v], 6, 7)
    name = paste(species, treat, id, sep="_")
    print(paste(v, "/", length(rawdata), "->", name, "row length=", dim(d)[1]))
    
    if(dim(d)[1] < (fps_original*60*30+1)){
      print(paste("row length is smaller than 54001, ", dim(d)[1]))
      cut_length = (dim(d)[1])
    }else {
      cut_length = (fps_original*60*30+1)
    }
    
    if(cut_length - dim(d)[1] > 6){
      for(i in 2:5){
        d[,i] = stats::filter(d[,i], rep(1/3, 3), sides = 2)
      }
      d <- d[seq(3,cut_length+2, fps_original/fps_analysis),]
    } else{
      for(i in 2:5){
        temp <- stats::filter(d[,i], rep(1/3, 3), sides = 2)
        temp[1:2] = d[1:2,i]
        temp[((dim(d)[1]-1):dim(d)[1])] = d[((dim(d)[1]-1):dim(d)[1]),i]
        d[,i] = temp
      }
      d <- d[seq(1,cut_length, fps_original/fps_analysis),]
    }

     # scale
    d[,1] <- d[,1]/fps_original
    
    
    
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
      ggsave(paste0(plot_loc_tra, name, ".png"), width = 4, height = 4)
    }
    
    # speed error check
    error_speed_thersh = 30
    error_check <- abs(diff(d$x0)) > 15 | abs(diff(d$x1)) > 15 | abs(diff(d$y0)) > 15 | abs(diff(d$y1)) > 15
    if( sum(error_check, na.rm=T)>0) {
      print(d[error_check,])
    }
  }  
  save(df_tra, df_ind, file = data_loc_fmt, compress="xz")
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
tandem.detect <- function(dataset = NULL, fps_analysis = 5){
  
  if(length(dataset)>0){
    load(paste0("data_fmt/data_trajectory_", dataset, ".rda"))
  } else{
    load("data_fmt/data_trajectory.rda")
  }
  
  name_list <- df_ind$name
  
  df_tandem = df_sep = df_sum <- NULL
  
  for(i in 1:length(name_list)){
    print(paste(i, "/", length(name_list), "->", name_list[i]))
    
    if(name_list[i] == "CF_FM_42" ||
       name_list[i] == "CF_MM_14" ||
       name_list[i] == "CF_MM_25" ||
       name_list[i] == "CG_FF_25" ||
       name_list[i] == "CG_FF_26" ||
       name_list[i] == "CG_FF_33" ||
       name_list[i] == "CG_FF_34" ||
       name_list[i] == "CG_FM_32" ||
       name_list[i] == "CG_FM_33" ||
       name_list[i] == "CG_FM_38" ||
       name_list[i] == "CG_FM_39" ||
       name_list[i] == "CG_FM_40" ||
       name_list[i] == "CG_MM_03" ||
       name_list[i] == "CG_MM_38" ||
       name_list[i] == "CG_MM_39"
    ){
      next;
    }
    #(CGFF22-24, maybe but not)
    #(CGFM7, maybe but not)
    #(CGFM29-30, maybe but not)
    #(CGFM46, maybe but not)
    #(CGMM5, maybe but not)

    df_temp <- subset(df_tra, name==name_list[i])
    df_temp$speed0 <- c(NA, sqrt(diff(df_temp$x0)^2 + diff(df_temp$y0)^2))*fps_analysis
    df_temp$speed1 <- c(NA, sqrt(diff(df_temp$x1)^2 + diff(df_temp$y1)^2))*fps_analysis
    sumbodylength <- sum(as.numeric(subset(df_ind, name==name_list[i])[,2:3]))
    ind_dis = sqrt( (df_temp$x0 - df_temp$x1)^2 + (df_temp$y0 - df_temp$y1)^2 )
    
    interaction <- ind_dis < sumbodylength
    non.interaction <- !interaction
    
    # remove short separation
    interaction <- !tandem.smoothing(non.interaction, 5*5) 
    # remove short tandem
    interaction <- tandem.smoothing(interaction, min.tandem.sec*5) 
    
    tandem = interaction
    
    # leader role
    # because we did not identify leader and follower during the tracking,
    # we estimated the leader role from the relative position and moving direction.
    # if dot product of these two vectors >0 -> partner is in their heading direction (follower)
    # if dot product of these two vectors <0 -> leader
    dx = diff(df_temp$x0)
    dy = diff(df_temp$y0)
    rx = (df_temp$x1 - df_temp$x0)[1:(dim(df_temp)[1]-1)]
    ry = (df_temp$y1 - df_temp$y0)[1:(dim(df_temp)[1]-1)]
    df_temp$leader0 = c((dx * rx + dy * ry) < 0, NA)
    
    dx = diff(df_temp$x1)
    dy = diff(df_temp$y1)
    rx = (df_temp$x0 - df_temp$x1)[1:(dim(df_temp)[1]-1)]
    ry = (df_temp$y0 - df_temp$y1)[1:(dim(df_temp)[1]-1)]
    df_temp$leader1 = c((dx * rx + dy * ry) < 0, NA)
    
    ## for each tandem event
    if(sum(tandem)>0){
      tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
      tan.sta <- which(tandem)[c(T, diff(which(tandem))>1)]
    }
    
    tan_duration <- tan.end - tan.sta
    sep_duration <- tan.sta[-1] - tan.end[-length(tan.end)]
    
    tan_cens <- rep(1, length(tan_duration))
    tan_cens[tan.sta == 1 | tan.end == 9000] <- 0
    
    speed_tandem = (df_temp$speed0[tandem]+df_temp$speed1[tandem])/2
    
    
    leader_clarity <- rep(0, length(tan_duration))
    if(length(tan_duration)>0){
      for(i_tan in 1:length(tan_duration)){
        tan_range <- tan.sta[i_tan]:tan.end[i_tan]
        leader0 = df_temp$leader0[tan_range]
        leader1 = df_temp$leader1[tan_range]
        leader_clarity[i_tan] = 
          abs(sum(leader0, na.rm = T) - sum(leader1, na.rm = T))/
          tan_duration[i_tan]
        if(leader_clarity[i_tan] < tandem_leader_thresh){
          tandem[tan_range] <- F
        }
      }
    }
    non_tandem_interaction <- interaction & !tandem
    
    ## rotation
    {
      dx0 = diff(df_temp$x0)
      dy0 = diff(df_temp$y0)
      dx1 = diff(df_temp$x1)
      dy1 = diff(df_temp$y1)
      
      leader.vec.length <- sqrt(dx0^2 + dy0^2)
      x0_vec <- dx0/leader.vec.length
      y0_vec <- dy0/leader.vec.length
      
      leader.vec.length <- sqrt(dx1^2 + dy1^2)
      x1_vec <- dx1/leader.vec.length
      y1_vec <- dy1/leader.vec.length
      
      x.center <- (df_temp$x0 + df_temp$x1)/2
      y.center <- (df_temp$y0 + df_temp$y1)/2
      
      lv <- length(df_temp$x0)
      
      x0.rotation.vec <- df_temp$x0[2:lv-1] - x.center[2:lv-1]
      y0.rotation.vec <- df_temp$y0[2:lv-1] - y.center[2:lv-1]
      leader.rotation.vec.length <- sqrt(x0.rotation.vec^2+y0.rotation.vec^2)
      x0.rotation.vec <- x0.rotation.vec/leader.rotation.vec.length
      y0.rotation.vec <- y0.rotation.vec/leader.rotation.vec.length
      
      x1.rotation.vec <- df_temp$x1[2:lv-1] - x.center[2:lv-1]
      y1.rotation.vec <- df_temp$y1[2:lv-1] - y.center[2:lv-1]
      leader.rotation.vec.length <- sqrt(x1.rotation.vec^2+y1.rotation.vec^2)
      x1.rotation.vec <- x1.rotation.vec/leader.rotation.vec.length
      y1.rotation.vec <- y1.rotation.vec/leader.rotation.vec.length
      
      x0.ang.moment <- x0_vec*x0.rotation.vec - y0_vec*y0.rotation.vec
      y0.ang.moment <- x0_vec*y0.rotation.vec + y0_vec*x0.rotation.vec
      x1.ang.moment <- x1_vec*x1.rotation.vec - y1_vec*y1.rotation.vec
      y1.ang.moment <- x1_vec*y1.rotation.vec + y1_vec*x1.rotation.vec
      
      leader.ang.moment.length <- sqrt(x0.ang.moment^2+y0.ang.moment^2)
      follower.ang.moment.length <- sqrt(x1.ang.moment^2+y1.ang.moment^2)
      
      rotation <- sqrt(
        (x0.ang.moment + x1.ang.moment)^2+
          (y0.ang.moment + y1.ang.moment)^2
      )/2
      
      while(sum(is.na(rotation))>0){
        rotation[which(is.na(rotation))] <- rotation[which(is.na(rotation))-1]
      }
    }
    competition <- non_tandem_interaction & ((c(rotation, F) > 0.75) & (ind_dis < sumbodylength*0.2))
    # remove short competition
    competition <- !tandem.smoothing(!competition, 2*5) 
    competition <- tandem.smoothing(competition, 2*5) 
    
    competition.sta <- NULL
    if(sum(competition)>0){
      competition.end <- which(competition)[c(diff(which(competition))>1,T)]
      competition.sta <- which(competition)[c(T, diff(which(competition))>1)]
      print(competition.sta/5)
      print((competition.end-competition.sta)/5)
    }
    
    # data summarize  
    df_sum_temp <- data.frame(
      name =  name_list[i],
      species = str_split(name_list[i], "_")[[1]][1],
      treatment = str_split(name_list[i], "_")[[1]][2],
      interaction_total_duration = sum(interaction) / fps_analysis,
      tandem_total_duration      = sum(tandem) / fps_analysis,
      non_tandem_interaction_total_duration = sum(non_tandem_interaction) / fps_analysis,
      rotation = sum(competition),
      num_comp = length(competition.sta),
      speed_tandem = mean(speed_tandem, na.rm=T)
    )
    
    df_tandem_temp <- data.frame(
      name =  name_list[i],
      species = str_split(name_list[i], "_")[[1]][1],
      treatment = str_split(name_list[i], "_")[[1]][2],
      tan_duration = tan_duration/fps_analysis,
      tan_cens,
      leader_clarity
    )
    
    df_sep_temp <- NULL
    if(length(sep_duration)>0){
      sep_thresh = 25
      for(i_sep in 1:length(sep_duration)){
        if(leader_clarity[i_sep] > tandem_leader_thresh){
          if(tan.end[i_sep]>(min.tandem.sec*fps_analysis)+1){
            if(sep_duration[i_sep] > sep_thresh){
              
              sep_range <- -(min.tandem.sec*fps_analysis):sep_thresh+tan.end[i_sep]
              tan_range <- tan.sta[i_sep]:tan.end[i_sep]
              speed0 = (df_temp$speed0[sep_range])
              speed1 = (df_temp$speed1[sep_range])
              leader0 = df_temp$leader0[tan_range]
              leader1 = df_temp$leader1[tan_range]
              
              if( sum(leader0) > sum(leader1) ){
                speed_leader   = speed0
                speed_follower = speed1
              } else {
                speed_leader   = speed1
                speed_follower = speed0
              }
              
              df_sep_temp_temp <- data.frame(
                name =  name_list[i],
                species = str_split(name_list[i], "_")[[1]][1],
                treatment = str_split(name_list[i], "_")[[1]][2],
                time = (-(min.tandem.sec*fps_analysis):sep_thresh)/fps_analysis,
                sep_event = i_sep,
                speed_leader, speed_follower
              )
              df_sep_temp <- rbind(df_sep_temp, df_sep_temp_temp)
            }
          }
        }
      }
    }
    
    df_sum <- rbind(df_sum, df_sum_temp)
    df_tandem <- rbind(df_tandem, df_tandem_temp)
    df_sep <- rbind(df_sep, df_sep_temp)
  }  
  
  
  if(length(dataset)>0){
    save(df_sum, df_tandem, df_sep, 
         file = paste0("data_fmt/df_tandem_fmt_", dataset, ".rda"))
  } else{
    save(df_sum, df_tandem, df_sep, file = "data_fmt/df_tandem_fmt.rda")
  }
}
#------------------------------------------------------------------------------#
