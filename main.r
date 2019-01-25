#############
##Libraries##
#############
library("ggplot2")

#################
##  Load Data  ##
#################

##List of CSV data frames
#folder <- readline("Where is the data from the spectrophotometer?  ")
in_dat <- list() 
dat <- list()

##List of files in current directory (ensure directory is correct)
#fils <- list.files(folder, full.names = TRUE)
files <- list.files(pattern = "\\.csv$")

##Loads CSV files
in_dat <- lapply(files, read.csv)
dat <- lapply(in_dat, function(x) x[complete.cases(x),])

##Variables are set to NULL once they are no longer required.
in_dat <- NULL
files <- NULL


#############
##FUNCTIONS##
#############

my.model <- function(){
  models <- list()
  for (i in 1:length(dat)) {
    models[[i]] <- glm(Intensity~Time, data = dat[[i]])
  }
  return(models)
}

my.LowestAICModel <- function(model){
  AIC_Val <- c()
  index_Val <- c()
  
  AICs <- data.frame(AIC_Val, index_Val) ##df of AIC and index of model
  
  ##Creates a data frame of AIC values with their respective index value
  for (i in 1:length(model)) {
    AICs <- rbind(AICs, data.frame(AIC_Val = AIC(model[[i]]), index_Val = i))
  }
  
  ##Finds lowest AIC value
  lowestAIC <- min(AICs$AIC_Val)
  
  ##Defines variable ret
  ret <- NULL
  
  ##Compares the AIC_Val to the lowest AIC,
  ##if the same then ret = index value of the model
  for (i in 1:length(model)) {
    ret[AICs$AIC_Val == lowestAIC & AICs$index_Val == i] <- i
  }
  
  ##For some reason complete.cases must be used as the second for loop creates a vector contain NA's this is then confuses the computer
  ##if the NA's are not removed.
  ret <- ret[complete.cases(ret)]
  return(paste("Model ", ret, " has the lowest AIC"))
}

my.graph <- function(){
  NormDat <- my.NormDat()
  GGP <- list()
    for (i in 1:length(NormDat[[1]])) {
      GGP[[i]] <- ggplot(NormDat[[1]][[i]], aes(Time, Intensity, colour = Intensity)) +
        geom_point(alpha = 0.2, show.legend = F) +
        xlab("Time (min)") +
        ylab("Intensity (A.U)") +
        theme_classic() +
        ylim(0,ceiling(NormDat[[2]] * 1.2))
    }
  do.call(grid.arrange, GGP)
  GGP <- NULL
}

my.NormDat <- function(){

  returnVal <- c()
  listDF <- list()
  normaliseInt <- list()
  ndat <- list()

  ndat <- lapply(dat, '[', c('Intensity'))
  yMins <- lapply(ndat, min)

  for (i in 1:length(yMins)) {
    normaliseInt[[i]] <- lapply(ndat[[i]], function(x) x/yMins[[i]])
    listDF[[i]] <- data.frame("Time" = dat[[i]]['Time'], 
                              "Intensity" = normaliseInt[[i]])
  }
  
  yMax <- max(sapply(unlist(normaliseInt), max))
  returnVal[[1]] <- listDF
  returnVal[[2]] <- yMax
  ndat <- NULL
  return(returnVal)
}