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
dat <- list()

##List of files in current directory (ensure directory is correct)
#fils <- list.files(folder, full.names = TRUE)
files <- list.files(pattern = "\\.csv$")

##Loads CSV files
dat <- lapply(files, read.csv)

##Fixes Column Names to be standarised
# dat <- lapply(dat, function(x) x[[3]] <- NULL)
dat <- lapply(dat, function(x) setNames(x, c('Time', 'Intensity')))

dat <- lapply(dat, function(x) x[complete.cases(x),])

##Variables are set to NULL once they are no longer required.
in_dat <- NULL
files <- NULL


#################
##  FUNCTIONS  ##
#################

my.NormDat <- function(){

  returnVal <- c()
  listDF <- list()
  normaliseInt <- list()
  ndat <- list()
  
  ##Creates a list of vectors containing the intensities
  ndat <- lapply(dat, '[', c('Intensity'))
  # ndat <- lapply(ndat, function(x){x <- as.numeric(x[[1]])})
  
  ##Creates a list of minimum alues from each vector in ndat
  yMins <- lapply(ndat, min)
  
  ##Loop creates new data frame with normalised intensities against time
  for (i in 1:length(yMins)) {
    normaliseInt[[i]] <- sapply(ndat[[i]], function(x) x/yMins[[i]])
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

my.SMA <- function(k) {      # k is the spand
  x = my.NormDat() ##Loads data frames
  x <- x[[1]] ##Removes unneeded element
  erg = list()
  for (n in 1:length(x)) {
    erg[[n]] <- data.frame("Time" = x[[n]]$Time, 
                           "Intensity" = smooth.ma(x[[n]]$Intensity, k))
  }
  return(erg)
}


##########

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

##Returns nls for each dataset
my.nls <- function(Yo){
  
  list_nls <- list()
  df <- my.NormDat()
  df <- df[[1]]
  
  for (i in 1:length(df)) {
    
    ##Set y
    y <- as.vector(df[[i]]$Intensity)
    x <- as.vector(df[[i]]$Time)
    
    ##a = Plateau-Yo 
    ##K is the rate constant
    ##Yo is the y value at the experimental time t=0
    
    ##Exponential Model for the graph
    m1 <- nls(y~a*(1-exp(K*x)) + Yo)
    
    list_nls[[i]] <- cor(y, predict(m1))
  }
  return(list_nls)
}

##Creates are graph of the data sets 
my.graph <- function(){
  NormDat <- my.NormDat()
  GGP <- list()
    for (i in 1:length(NormDat[[1]])) {
      GGP[[i]] <- ggplot(NormDat[[1]][[i]], 
                         aes(Time, Intensity, colour = Intensity)) +
        geom_point(alpha = 0.2, show.legend = F) +
        xlab("Time (min)") +
        ylab("Normalised Intensity") +
        theme_classic() +
        ylim(0,ceiling(NormDat[[2]] * 1.2))
    }
  do.call(grid.arrange, GGP)
  GGP <- NULL
}

# ##Plot a graph of the gradients of a data set.
# my.difgraph <- function(k){
#   NormDat <- data.frame("Time" = my.NormDat()[[1]]my.gradient(k)
#   GGP <- list()
#   for (i in 1:length(NormDat)) {
#     GGP[[i]] <- ggplot(NormDat[[i]], 
#                        aes(Time, Intensity, colour = Intensity)) +
#       geom_point(alpha = 0.2, show.legend = F) +
#       xlab("Time (min)") +
#       ylab("Normalised Intensity") +
#       theme_classic()
#   }
#   do.call(grid.arrange, GGP)
#   GGP <- NULL
# }

##Creates a graph with the SMA line drawn too. K = SMA span. 175-200 works best
my.SMAgraph <- function(k){
  GGP <- list()
  index <- 1
  if (is.vector(k) && is.numeric(k)) {
    NormDat <- my.NormDat()
    for (a in 1:length(k)) {
      SMANormDat <- my.SMA(k[[a]])
      for (i in 1:length(SMANormDat)) {
        ##Creates a list of ggplot plots
        GGP[[index]] <- ggplot(NormDat[[1]][[i]], 
                                   aes(Time, Intensity)) +
          geom_point(alpha = 0.2, aes(color = "myGrey")) +
          geom_line(data = SMANormDat[[i]], aes(color = "myRed")) +
          xlab("Time (min)") +
          ylab("Intensity") +
          theme_classic() +
          ylim(0,ceiling(NormDat[[2]] * 1.2)) +
          scale_color_manual(name = paste("Repeat: ", i),
                             values = c(myGrey = "grey", myRed = "red"),
                             labels = c("Normalised \nIntensity", 
                                        paste("SMA (Span = ", k[[a]], ")")))
        index <- index + 1
      }
    }
    # do.call(grid.arrange, GGP)
    grid.arrange(grobs = GGP)
  }else if (is.numeric(k)) {
    NormDat <- my.NormDat()
    SMANormDat <- my.SMA(k)
    for (i in 1:length(SMANormDat)) {
      ##Creates a list of ggplot plots
      GGP[[i]] <- ggplot(NormDat[[1]][[i]], aes(Time, Intensity)) +
        geom_point(alpha = 0.2, aes(color = "myGrey")) +
        geom_line(data = SMANormDat[[i]], aes(color = "myRed")) +
        xlab("Time (min)") +
        ylab("Intensity") +
        theme_classic() +
        ylim(0,ceiling(NormDat[[2]] * 1.2)) +
        scale_color_manual(name = "",
                           values = c(myGrey = "grey", myRed = "red"),
                           labels = c("Normalised \nIntensity", 
                                      paste("SMA (Span = ", k, ")")))
    }
    
    ##Arranges all plots within GGP on a single grob
    do.call(grid.arrange, GGP)
    GGP <- NULL
    
  } else {
    print("Please ensure that K is either numeric or a numerical vector")
    break
  }
}

##Returns gradients for given list of DFs, f.
my.gradient <- function(f) {
  grads <- list()
  for (i in 1:length(f)) {
    grads[[i]] <- gradient(f[[i]]$Intensity)
  }
  return(grads)
}

##Returns a dataframe of the Average between all OK dataframes. *Incomplete*
my.Average <- function(f, k = NULL){
  f <- lapply(f, function(x){unlist(x$Intensity)})
  fsum <- c(0)
  for (i in 1:length(f)){
    if (k == NULL){
      fsum = fsum + f[[i]]
    }else if(k[[i]]){
      fsum = fsum + f[[i]]
    }
  }
  return(fsum/length(f))
}