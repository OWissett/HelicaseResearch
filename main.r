########################
##      Libraries     ##
########################

##Nice little function that I found that will install packages if they are not 
##installed.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

libs <- c("ggplot2", 
          "gridExtra",
          "itsmr",
          "pracma",
          "git2r",
          "car",
          "minpack.lm")

ipak(libs)

########################
##      Load Data     ##
########################

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


########################
##  DATA  PROCESSING  ##
########################

##Normalises data by dividing by lowest intensity
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

##Normalises data by dividing by max intensity resulting in data being [0, 1]
my.NormDat1 <- function(){
  
  returnVal <- c()
  listDF <- list()
  normaliseInt <- list()
  ndat <- list()
  
  ##Creates a list of vectors containing the intensities
  ndat <- lapply(dat, '[', c('Intensity'))
  # ndat <- lapply(ndat, function(x){x <- as.numeric(x[[1]])})
  
  ##Creates a list of max values from each vector in ndat
  yMax <- lapply(ndat, max)
  
  ##Loop creates new data frame with normalised intensities against time
  for (i in 1:length(yMax)) {
    normaliseInt[[i]] <- sapply(ndat[[i]], function(x) x/yMax[[i]])
    listDF[[i]] <- data.frame("Time" = dat[[i]]['Time'], 
                              "Intensity" = normaliseInt[[i]])
  }
  
  ##Nulls ndat once no longer needed (this may be redundant...)
  ##May possibly decrease memory used...
  ndat <- NULL
  
  return(listDF)
}


##Creates a smooth moving average for the data set
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

##Returns nls for each dataset
my.nls <- function(){
  
  list_nls <- list()
  df <- my.NormDat()
  df <- df[[1]]
  
  for (i in 1:length(df)) { 
    
    y <- as.vector(df[[i]]$Intensity)
    x <- as.vector(df[[i]]$Time)
    dat2 <- data.frame(x = x, y = y)
    
    eq1 <- function(x, a, K, Y0){Y0 + (a - Y0)*(1 - exp(-K*x))}
    
    eq2 <- function(x, a, PropFast, Kfast, Kslow, Y0){
      SpanFast <- (a - Y0)*PropFast
      SpanSlow <- (a - Y0)*(1 - PropFast)
      return(Y0 + SpanFast*(1 - exp(-Kfast*x)) + SpanSlow*(1 - exp(-Kslow*x)))
    }
    
    eq3 <- function(x, A, B, k.1, k.2, Y0){
      Y0 + A*(1 - exp(-k.1*x)) + B*(1 - exp(-k.2*x))}
     ## Derived by convoluting eq2 then analysing resulting function (assuming that Y0 ~ 0)
    eq4 <- function(x, k.1, k.2){1 + (k.2/(k.2-k.1))*exp(-k.2*x) + (k.1/(k.2-k.1))*exp(-k.1*x)}
    
    ##a = Plateau
    ##K is the rate constant
    ##Y0 is the y value at the experimental time t=0
    
    ###############
    ##Exponential##
    ###############
    
    ##Linearised Fit to find starting values.
    lin_fit1 <- nls(log(y) ~ log(eq1(x,a,K,Y0)), 
                   dat2,
                   start = c(a = max(dat2$y)/2,
                             K = 0.1,
                             Y0 = 1))
    
    ##Exponential Model for the graph
    nls_fit1 <- nls(y ~ eq1(x, a, K, Y0), 
                   dat2, 
                   start = c(a = coef(lin_fit1)[[1]],
                             K = coef(lin_fit1)[[2]],
                             Y0 = coef(lin_fit1)[[3]]))
    
    #################
    ##BiExponential##
    #################
    
    # a0 <- max(dat$y) / 2
    # 
    # lin_fit2 <- lm(log(y - a0) ~ x, data = dat2)
    # start <- list(b = exp(coef(lin_fit2)[1]), 
    #               b = coef(lin_fit2)[2],
    #               )
    
    # lin_fit2 <- nlsLM(log(y) ~ log(eq2(x,a,PropFast,Kfast,Kslow,Y0)),
    #                dat2,
    #                start = c(a = coef(nls_fit1)[[1]],
    #                          PropFast = 0.5,
    #                          Kfast = coef(nls_fit1)[[2]],
    #                          Kslow = coef(nls_fit1)[[2]],
    #                          Y0 = coef(nls_fit1)[[3]]),
    #                lower = c(0,0,0,0,0),
    #                upper = c(max(dat2$y), 1, Inf, Inf, Inf))
    # 
    # 
    # ##Biexponential Model for the graph
    # nls_fit2 <- nlsLM(y ~ eq2(x,a,PropFast,Kfast,Kslow,Y0),
    #                 dat2,
    #                 start = c(a = coef(lin_fit2)[[1]],
    #                           PropFast = coef(lin_fit2)[[2]],
    #                           Kfast = coef(lin_fit2)[[3]],
    #                           Kslow = coef(lin_fit2)[[4]],
    #                           Y0 = coef(lin_fit2)[[5]]))
    
    lin_fit3 <- nlsLM(log(y) ~ log(eq2(x, A, B, k.1, k.2, Y0)),
                      dat2,
                      start = c(A = max(dat2$y)/4,
                                B = max(dat2$y)/4,
                                k.1 = coef(nls_fit1)[[2]],
                                k.2 = coef(nls_fit1)[[2]],
                                Y0 = coef(nls_fit1)[[3]]),
                      lower = c(-Inf,-Inf,0,0,0),
                      upper = c(max(dat2$y), 1, Inf, Inf, Inf))
    
    
    # ##Biexponential Model for the graph
    # nls_fit3 <- nlsLM(y ~ eq3(x, A, B, k.1, k.2, Y0),
    #                   dat2,
    #                   start = c(a = coef(lin_fit3)[[1]],
    #                             PropFast = coef(lin_fit3)[[2]],
    #                             Kfast = coef(lin_fit3)[[3]],
    #                             Kslow = coef(lin_fit3)[[4]],
    #                             Y0 = coef(lin_fit3)[[5]]))

    list_nls[[i]] <- list("fit1" = nls_fit1, "fit3" = lin_fit3)
  }
  return(list_nls)
}

########################
## GRAPHING FUNCTIONS ##
########################

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

##Plot a graph of the gradients of a data set. *INCOMPLETE*
my.difgraph <- function(k){
  NormGrad <- data.frame("Time" = my.NormDat()[[1]]$Time,
                        "Gradient" = my.gradient(k))
  GGP <- list()
  for (i in 1:length(NormDat)) {
    GGP[[i]] <- ggplot(NormDat[[i]],
                       aes(Time, Intensity, colour = Intensity)) +
      geom_point(alpha = 0.2, show.legend = F) +
      xlab("Time (min)") +
      ylab("Normalised Intensity") +
      theme_classic()
  }
  do.call(grid.arrange, GGP)
  GGP <- NULL
}

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

##Creates a graph with the NLS line drawn.
my.NLSgraph <- function(){
  
  df <- my.NormDat()[[1]]
  GGP <- list()
  for (i in 1:length(df)) {
    
    ##Creates local dataframe to be used in calculating the NLS
    y <- as.vector(df[[i]]$Intensity)
    x <- as.vector(df[[i]]$Time)
    dat2 <- data.frame(x = x, y = y)
    eq1 <- function(x, a, K, Y0){Y0 + (a-Y0)*(1-exp(-K*x))}
    
    ##a = Plateau
    ##K is the rate constant
    ##Y0 is the y value at the experimental time t=0
    
    ##Linearised Fit to find starting values.
    lin_fit <- nls(log(y) ~ log(eq1(x,a,K,Y0)), 
                   dat2,
                   start = c(a = max(dat2$y)/2,
                             K = 0.1,
                             Y0 = 1))
    
    ##Creates the graphic plot
    GGP[[i]] <- ggplot(df[[i]],
                       aes(Time, Intensity)) +
      geom_point(alpha = 0.2, show.legend = F) +
      geom_smooth(method="nls",
                  formula = y~Y0 + (a-Y0)*(1-exp(-K*x)),
                  method.args = list(start = c(a = coef(lin_fit)[[1]],
                                               K = coef(lin_fit)[[2]],
                                               Y0 = coef(lin_fit)[[3]])),
                  se=F,
                  colour = "red",
                  aes(linetype = 2)) +
      xlab("Time (min)") +
      ylab("Normalised Intensity") +
      theme_classic()
  }
  do.call(grid.arrange, GGP)
  GGP <- NULL
}

########################
##  STATS  FUNCTIONS  ##
########################

##Returns gradients for given list of DFs, f.
my.gradient <- function(f) {
  grads <- list()
  for (i in 1:length(f)) {
    grads[[i]] <- gradient(f[[i]]$Intensity)
  }
  return(grads)
}

##Returns a dataframe of the Average between all OK dataframes. *Incomplete*
my.Average <- function(f, test_vector = NULL){
  f <- lapply(f, function(x){unlist(x$Intensity)})
  fsum <- c(0)
  for (i in 1:length(test_vector)){
    if (is.null(test_vector)){
      fsum = fsum + f[[i]]
    }else if(k[[i]] == T){
      fsum = fsum + f[[i]]
    }
  }
  return(fsum/length(f))
}

##Returns the model with the lowest AIC value from a list of models
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

##Test data sets rates and data span
my.StatsTest <- function(){
  list_nls <- my.nls()
  
  for(i in 1:length(list_nls)){
    co <- coef(summary(list_nls[[i]]))
    
    
  }
  
  
}




