# Coptotermes homosexual tandem analysis
# preprocess.R
# N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for preprocess all data for movement analysis
#------------------------------------------------------------------------------#

rm(list = ls())

#------------------------------------------------------------------------------#
{
  # packages
  library(data.table)
  library(stringr)
  
  # constants
  {
    min.sep.sec                 <- 2
    min.tandem.sec              <- 5
    threshold.moving.dis.tandem <- 30
    termite.body.length         <- c(9, 9, 8, 8)
    names(termite.body.length)  <- c("CF-F", "CF-M", "CG-F", "CG-M")
  }
  
  process.all()
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
process.all <- function(){
  tandem.detect()
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# This function reads allTrajectoryData.rda files, then
# 1) Determine if a pair is doing tandem or not for each frame.
# 2) save data.frame for each tandem event (),
#    and for summary for each pair (df_sum_tandem.rda)
#------------------------------------------------------------------------------#
tandem.detect <- function(){
  
  df = readRDS("data/rda/AllData-scaled.rda")
  
  # calcurate step length
  df$sl <- c(NA, sqrt( diff(df$x)^2 + diff(df$y)^2))
  df$sl[df$time == 0] <- NA
  
  ind.names <- unique(df$name) 
  df.tan.time = res.sep.time = df.pair <- NULL
  
  for(i in 1:length(ind.names)){
    print(paste(i, "/", length(ind.names), "->", ind.names[i]))
    
    df.temp <- df[df$name == ind.names[i],]

    df.ind1 <- df.temp[df.temp$ind == 1,]
    df.ind2 <- df.temp[df.temp$ind == 2,]
    
    x.ind1  = df.ind1$x
    y.ind1  = df.ind1$y
    x.ind2  = df.ind2$x
    y.ind2  = df.ind2$y
    sl.ind1 = df.ind1$sl
    sl.ind2 = df.ind2$sl
    ind.dis = sqrt( (x.ind1 - x.ind2)^2 + (y.ind1 - y.ind2)^2 )
    
    # definition of tandem running
    {
      # 1. within "interaction.threshold" for more than "min.tandem.sec"
      #    also separation less than "min.tandem.sec" do not count as separation
      #    interaction.threshold = 1.3 * bodylength (averaged between partners)
      {
        treat   = str_sub(ind.names[i], 4,5)
        species = str_sub(ind.names[i], 1,2)
        if( species == "CF"){
          if( treat == "FF"){
            interaction.threshold = termite.body.length[1] * 1.3
          } else if( treat == "MM"){
            interaction.threshold = termite.body.length[2] * 1.3
          } else {
            interaction.threshold = mean(termite.body.length[1:2]) * 1.3
          }
        } else {
          if( treat == "FF"){
            interaction.threshold = termite.body.length[3] * 1.3
          } else if( treat == "MM"){
            interaction.threshold = termite.body.length[4] * 1.3
          } else {
            interaction.threshold = mean(termite.body.length[3:4]) * 1.3
          }
        }
        interaction     <- ind.dis <= interaction.threshold
        non.interaction <- !interaction
        
        # remove short separation
        tandem          <- !tandem.smoothing(non.interaction, min.sep.sec*5) 
        
        # remove short tandem
        tandem          <- tandem.smoothing(tandem, min.tandem.sec*5) 
        
        
      }
      
      # 2. tandem: both need to move more than "threshold.moving.dis.tandem" during interaction
      if(sum(tandem)>0){
        tan.timing <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan.end    <- which(tandem)[c(diff(which(tandem))>1,T)]
        for(j in 1:length(tan.timing)){
          if(sum(sl.ind1[tan.timing[j]:tan.end[j]], na.rm=T)< threshold.moving.dis.tandem){
            tandem[tan.timing[j]:tan.end[j]] <- F
          }
          if(sum(sl.ind2[tan.timing[j]:tan.end[j]], na.rm=T)< threshold.moving.dis.tandem){
            tandem[tan.timing[j]:tan.end[j]] <- F
          }
        }
      }
    }
    
    #---------------------------
    # may be removed in the final version (no plan to use this)
    # get scheme
    {
      scheme <- tandem
      scheme[interaction] <- "i"
      scheme[non.interaction] <- "r"
      scheme[tandem] <- "t"
      
      if(sum(tandem)>0){
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        if(tan.end[1]==length(tandem)){  }else {
          sep.begin <- tan.end[tan.end<length(tandem)]+1
          
          ## sep until next interaction
          tan.timing <- which(interaction)[c(T, diff(which(interaction))>1)]
          if(length(tan.timing)==1){ sep.end <- NULL} else {
            sep.end <- tan.timing[2:(length(tan.timing))]-1}
          if(length(sep.begin) - length(sep.end) == 1){ sep.end <- c(sep.end,dim(tandem)[1])}
          for(j in 1:length(sep.begin)){
            for(k in 1:length(sep.end)){
              if(sep.end[k]>sep.begin[j]){
                scheme[sep.begin[j]:sep.end[k]] <- "s"
                break;
              }
            }
          }
        }
      }
    }
    df.scheme = data.frame(
      frame = df.temp$frame,
      sl1 = sl.ind1,
      sl2 = sl.ind2,
      name = ind.names[i],
      scheme
    )
    ggplot(df.scheme) +
      geom_path(aes(x=frame, y=sl1, col=2))+
      geom_path(aes(x=frame, y=sl2, col=3))+
      ggtitle(ind.names[i])
    #---------------------------
    
    # tandem survival time
    {
      if(sum(tandem)>0){
        tan.timing <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan.end    <- which(tandem)[c(diff(which(tandem))>1,T)]
        df.tandem <-  data.frame(Video    = df.temp$name[1],
                                 Species  = species,
                                 Treat    = treat,
                                 Tan.time = (tan.end - tan.timing + 1)/5,
                                 Cens     = !(tan.timing[j]==1 || tan.end[j]==length(tandem))
                                 )
        df.tan.time <- rbind(df.tan.time, df.tandem)
      }
    }
    
    #---------------------------
    # may be removed in the final version (no plan to use this)
    ## sep search survival time
    {
      if(sum(scheme=="s")>0){
        Video = Species = Treat = Sep.time = Cens <- NULL
        sep.timing <- which(scheme=="s")[c(T, diff(which(scheme=="s"))>1)]
        sep.end <- which(scheme=="s")[c(diff(which(scheme=="s"))>1,T)]
        for(j in 1:length(sep.timing)){
          Sep.time = c(Sep.time, sep.end[j] - sep.timing[j] + 1)
          if(sep.timing[j]==1 || sep.end[j]==length(tandem)){
            Cens = c(Cens, 0)
          } else {
            Cens = c(Cens, 1)
          }
          Video = c(Video, df.temp$name[1])
          Species = c(Species, species)
          Treat = c(Treat, treat)
        }
        res.sep.time <- rbind(res.sep.time, data.frame(Video, Species, Treat, Sep.time=Sep.time/5, Cens))
      }
    }
    #---------------------------
    
    # pair base tandem data
    {
      df.pair <- rbind(df.pair, 
                       data.frame(Video   = df.temp$name[1],
                                  Species = species,
                                  Treat   = treat,
                                  Tandem  = sum(tandem)/5))
    }
    
  }
  
  save(df.tan.time, df.pair, file = "data/rda/df_tandem.rda")
  
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