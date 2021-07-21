## Data analysis for Coptotermes homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file plots trajectories to check tracking errors
#------------------------------------------------------------------------------#

rm(list = ls())
raw_trajectory()
scaled_trajectory()

#------------------------------------------------------------------------------#
{
  today <- Sys.Date()
  
  library(data.table)
  library(stringr)
  
  library(ggplot2)
  require(grid)
  require(gridExtra)

  ## parameters
  fps <- 5
  arena.size <- 150 # mm
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# reads all csv files and 1) output raw trajectories in .png file
# 2) convert them in 5FPS and R data.frame saved in .rda files
#------------------------------------------------------------------------------#
raw_trajectory <- function(Plot = T, Dataframe = T){
  options(warn = 0)

  # data
  rawdata <- list.files("data/raw", full.names = TRUE, pattern = ".csv")
  dataname <- list.files("data/raw", full.names = FALSE, pattern = ".csv")
  
  df <- data.frame()
  
  # plots
  for(v in 1:length(rawdata)){
    
    # file info
    d <- data.frame(fread(rawdata[v], header=T))
    d <- d[1:(30*60*30+1),]
    
    species = substr(dataname[v], 1, 2)
    treat = substr(dataname[v], 4, 5)
    id = substr(dataname[v], 6, 7)
    name <- paste0(species, "-", treat, "-", id)
    print(paste(v, "/", length(rawdata), "->", name))

    # plot
    if(Plot){
      png(paste0("output/trajectory/raw/", str_replace(dataname[v], "csv", "png")))
      par(mfrow=c(1,1), pin=c(3,3))
      plot( d[,2], d[,3], type="l", col=1, xlim=c(0,1200), ylim=c(0,1200),
            main=paste(name),
            xlab="x(mm)", ylab="y(mm)")
      dev.off()
    }
    
    # datafrmae
    d[,1] <- d[,1]/30
    d <- d[seq(1,54000,6),]
    if(Dataframe){
      dftemp1 <- data.frame(name, time = d[,1], x=d[,2], y=d[,3])
      df <- rbind(df,dftemp1)
      dftemp2 <- data.frame(name, time = d[,1], x=d[,4], y=d[,5])
      df <- rbind(df,dftemp2)
    }
  }
  save(df, file = "data/rda/AllData.rda")
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# assuming that termites moved throughout the arena, scale the movement range as arena size
# scaling is performed for each observed arena (pool tandem and separation for each individual)
# also check speed time development to detect tracking error
#------------------------------------------------------------------------------#
scaled_trajectory <- function(Plot = T, Dataframe = T){
  
  load("data/rda/AllData.rda")
  
  ind.names <- unique(df$name) 
  df <- na.omit(df)
  
  for(v in 1:length(ind.names)){
    print(paste(v, "/", length(ind.names), "->", ind.names[v]))
    
    df.temp <- df[df$name == ind.names[v],]
    
    x <- df.temp$x
    y <- df.temp$y
    xL <- max(x) - min(x)
    yL <- max(y) - min(y)
    df.temp$x <- (x-min(x))/xL * arena.size
    df.temp$y <- (y-min(y))/yL * arena.size
    
    # plot
    if(Plot){
      p1 <- ggplot(data=df.temp, aes(x=x, y=y, col=as.factor(ind), group=as.factor(ind))) +
        geom_path() + 
        theme_bw() +
        theme(legend.position = "none", aspect.ratio = 1)+ 
        xlim(c(0,140)) + ylim(c(0,140)) +
        ggtitle(paste(v, ind.names[v]))
      df.temp$speed = c(NA, sqrt( diff(df.temp$x)^2 + diff(df.temp$y)^2 ) )
      p2 <- ggplot(data=df.temp, aes(x=time, y=speed, col=as.factor(ind), group=as.factor(ind))) +
        geom_path() + 
        theme_bw() +
        theme(legend.position = "none", aspect.ratio = 1)+ 
        coord_cartesian(xlim=c(0,1800), ylim=c(0, 30)) +
        ggtitle(paste(v, ind.names[v]))
      ggsave(paste0("output/trajectory/scaled/", df.temp$name[1], ".png"), 
             plot = grid.arrange(p1, p2, nrow = 1), width = 6, height = 3)
    }
    
    # data.frame
    df[df$name == ind.names[v],] <- df.temp[,1:(dim(df.temp)[2]-1)]
    if(length(df.temp[df.temp$speed > 20  & df.temp$time>0,1])>0){
      print(df.temp[df.temp$speed > 20  & df.temp$time>0,])
    }
    
  }  
  save(df, file = "data/rda/AllData-scaled.rda")
}
#------------------------------------------------------------------------------#
