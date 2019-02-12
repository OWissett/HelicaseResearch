#############
##Libraries##
#############
library("ggplot2")
library("gridExtra")
library("itsmr")
library("pracma")
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

my.NormDat <- function(){
  
  ##Defines local vairables
  returnVal <- c()
  listDF <- list()
  normaliseInt <- list()
  ndat <- list()
  
  ##Creates a list of vectors containing the intensities
  ndat <- lapply(dat, '[', c('Intensity'))
  
  ##Creates a list of minimum alues from each vector in ndat
  yMins <- lapply(ndat, min)
  
  ##Loop creates new data frame with normalised intensities against time
  for (i in 1:length(yMins)) {
    normaliseInt[[i]] <- lapply(ndat[[i]], function(x) x/yMins[[i]])
    listDF[[i]] <- data.frame("Time" = dat[[i]]['Time'], 
                              "Intensity" = normaliseInt[[i]])
  }
  
  ##Nulls ndat once no longer needed (this may be redundant...)
  ##May possibly decrease memory used...
  ndat <- NULL
  
  ##Finds max light intenisity out of all values in all vectors, returns an int
  yMax <- max(sapply(unlist(normaliseInt), max))
  
  ##Assigns the values that arer returned by the funtion to the retrun array
  returnVal[[1]] <- listDF
  returnVal[[2]] <- yMax
  
  return(returnVal)
}

my.SMA <- function(x, k) {      # k is the span, x is the list of vectors
  # x = my.NormDat() ##Loads data frames
  # x <- x[[1]] ##Removes unneeded element
  erg = list()
  for (n in 1:length(x)) {
    erg[[n]] <- data.frame("Time" = x[[n]]$Time, 
                           "Intensity" = smooth.ma(x[[n]]$Intensity, k))
  }
  return(erg)
}

#########
# my.model <- function(){
#   models <- list()
#   normDat <- my.NormDat()
#   for (i in 1:length(normDat[[1]])) {
#     models[[i]] <- glm(Intensity~Time, data = normDat[[1]][[i]])
#   }
#   return(models)
# }

# my.LowestAICModel <- function(model){
#   AIC_Val <- c()
#   index_Val <- c()
#   
#   AICs <- data.frame(AIC_Val, index_Val) ##df of AIC and index of model
#   
#   ##Creates a data frame of AIC values with their respective index value
#   for (i in 1:length(model)) {
#     AICs <- rbind(AICs, data.frame(AIC_Val = AIC(model[[i]]), index_Val = i))
#   }
#   
#   ##Finds lowest AIC value
#   lowestAIC <- min(AICs$AIC_Val)
#   
#   ##Defines variable ret
#   ret <- NULL
#   
#   ##Compares the AIC_Val to the lowest AIC,
#   ##if the same then ret = index value of the model
#   for (i in 1:length(model)) {
#     ret[AICs$AIC_Val == lowestAIC & AICs$index_Val == i] <- i
#   }
#   
#   ##For some reason complete.cases must be used as the second for loop creates a vector contain NA's this is then confuses the computer
#   ##if the NA's are not removed.
#   ret <- ret[complete.cases(ret)]
#   return(paste("Model ", ret, " has the lowest AIC"))
# }

#########

##Creates are graph of the data sets 
my.graph <- function(){
  NormDat <- my.NormDat()
  GGP <- list()
    for (i in 1:length(NormDat[[1]])) {
      GGP[[i]] <- ggplot(NormDat[[1]][[i]], 
                         aes(Time, Intensity, colour = Intensity)) +
        geom_point(alpha = 0.2, show.legend = F) +
        xlab("Time (min)") +
        ylab("Intensity (A.U)") +
        theme_classic() +
        ylim(0,ceiling(NormDat[[2]] * 1.2))
    }
  do.call(grid.arrange, GGP)
  GGP <- NULL
}

##Creates a graph with the SMA line drawn too. K = SMA span. 175-200 works best
my.SMAgraph <- function(k){
  NormDat <- my.NormDat()
  SMANormDat <- my.SMA(k)
  ldf <- list()
  GGP <- list()
  for (i in 1:length(SMANormDat)) {
    ##Creates a list of ggplot plots
    GGP[[i]] <- ggplot(NormDat[[1]][[i]], aes(Time, Intensity)) +
      geom_point(alpha = 0.2, show.legend = T, color = "grey") +
      geom_line(data = SMANormDat[[i]], color = "red") +
      xlab("Time (min)") +
      ylab("Intensity (A.U)") +
      theme_classic() +
      ylim(0,ceiling(NormDat[[2]] * 1.2))
  }
  ##Arranges all plots within GGP on a single grob
  do.call(grid.arrange, GGP)
  GGP <- NULL
}

my.gradient <- function(f) {
  grads <- list()
  for (i in 1:length(f)) {
    grads[[i]] <- gradient(f[[i]]$Intensity)
  }
  return(grads)
}

my.Average <- function(){
  
  normDat <- my.NormDat()
  normDat[[2]] <- NULL
  
  return()
}