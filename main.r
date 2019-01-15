#############
##Libraries##
#############
library("ggplot2")

#################
##  Load Data  ##
#################

##List of CSV data frames
folder <- readline("Where is the data from the spectrophotometer?")
dat <- list() 

##List of files in current directory (ensure directory is correct)
fils <- list.files(folder)
fils
##Loads CSV files
for(i in 1:length(fils)){
  ##Reads the CSV file
  dat[[i]] <- read.csv(fils[[i]], head=T, sep=",")
  ##Removes whitespace between rows, this decreases object size.
  dat[[i]] <- dat[[i]][complete.cases(dat[[i]]),]
  
}


#############
##FUNCTIONS##
#############

my.model <- function(){
  models <- list()
  for(i in 1:length(dat)){
    models[[i]] <- glm(Intensity~Time, data = dat[[i]])
  }
  return(models)
}

my.LowestAICModel <- function(model){
  AIC_Val <- c()
  index_Val <- c()
  
  AICs <- data.frame(AIC_Val, index_Val) ##df of AIC and index of model
  
  ##Creates a data frame of AIC values with their respective index value
  for(i in 1:length(model)){
    AICs <- rbind(AICs, data.frame(AIC_Val=AIC(model[[i]]), index_Val=i))
  }
  
  ##Finds lowest AIC value
  lowestAIC <- min(AICs$AIC_Val)
  
  ##Defines variable ret
  ret <- NULL
  
  ##Compares the AIC_Val to the lowest AIC, if the same then ret = index value of the model
  for(i in 1:length(model)){
    ret[AICs$AIC_Val==lowestAIC & AICs$index_Val==i] <- i
  }
  
  ##For some reason complete.cases must be used as the second for loop creates a vector contain NA's this is then confuses the computer
  ##if the NA's are not removed.
  ret <- ret[complete.cases(ret)]
  return(paste("Model ", ret, " has the lowest AIC"))
}

my.graph <- function(){
  
  mods <- my.model()
  for(i in 1:length(dat)){
    ggplot(dat[[i]], aes(Time, Intensity, colour=Intensity)) + 
      geom_point(alpha=0.2, show.legend = F) +
      xlab("Time (min)") +
      ylab("Intensity (A.U)") +
      theme_classic() +
      xlim(0, 16) +
      ylim(25, 250) 
  }

}

