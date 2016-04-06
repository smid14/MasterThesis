#Third Party Libraries
require(vars)
require(fastSOM)
require(TTR)
require(xts)
library(dygraphs)
library(igraph)


# File that makes the calculation for the static and dynamic connectedness Analysis:

#-----------------------------Part 1: Reads and Transform the data------------------------------------------------------------

#----Reads the Data------
#---crisis: checks if the shorter or longer sample should be calculated. Standard is longer sample
#---RETURN: an xts data structure with all bond yields
readData <- function(crisis = FALSE){
  #---Reads the bond yields
  aus <- read.csv(file="/Input/Australia.csv", stringsAsFactors=FALSE, sep = ",")
  aut <- read.csv(file="/Input/Austria.csv", stringsAsFactors=FALSE, sep = ",")
  bel <- read.csv(file="/Input/Belgium.csv", stringsAsFactors=FALSE, sep = ",")
  can <- read.csv(file="/Input/Canada.csv", stringsAsFactors=FALSE, sep = ",")
  den <- read.csv(file="/Input/Denmark.csv", stringsAsFactors=FALSE, sep = ",")
  fin <- read.csv(file="/Input/Finland.csv", stringsAsFactors=FALSE, sep = ",")
  fra <- read.csv(file="/Input/France.csv", stringsAsFactors=FALSE, sep = ",")
  ger <- read.csv(file="/Input/Germany.csv", stringsAsFactors=FALSE, sep = ",")
  gre <- read.csv(file="/Input/Greece.csv", stringsAsFactors=FALSE, sep = ",")
  ire <- read.csv(file="/Input/Ireland.csv", stringsAsFactors=FALSE, sep = ",")
  ita <- read.csv(file="/Input/Italy.csv", stringsAsFactors=FALSE, sep = ",")
  jpn <- read.csv(file="/Input/Japan.csv", stringsAsFactors=FALSE, sep = ",")
  net <- read.csv(file="/Input/Netherlands.csv", stringsAsFactors=FALSE, sep = ",")
  nor <- read.csv(file="/Input/Norway.csv", stringsAsFactors=FALSE, sep = ",")
  por <- read.csv(file="/Input/Portugal.csv", stringsAsFactors=FALSE, sep = ",")
  spa <- read.csv(file="/Input/Spain.csv", stringsAsFactors=FALSE, sep = ",")
  swe <- read.csv(file="/Input/Sweden.csv", stringsAsFactors=FALSE, sep = ",")
  usa <- read.csv(file="/Input/Usa.csv", stringsAsFactors=FALSE, sep = ",")
  uk <- read.csv(file="/Input/Uk.csv", stringsAsFactors=FALSE, sep = ",")
  
  #---Just a data structure transformation
  xts.aus <- xts(aus$Value, as.Date(aus$Date))
  xts.aut <- xts(aut$Value, as.Date(aut$Date))
  xts.bel <- xts(bel$Value, as.Date(bel$Date))
  xts.can <- xts(can$Value, as.Date(can$Date))
  xts.den <- xts(den$Value, as.Date(den$Date))
  xts.fin <- xts(fin$Value, as.Date(fin$Date))
  xts.fra <- xts(fra$Value, as.Date(fra$Date))
  xts.ger <- xts(ger$Value, as.Date(ger$Date))
  xts.gre <- xts(gre$Value, as.Date(gre$Date))
  xts.ire <- xts(ire$Value, as.Date(ire$Date))
  xts.ita <- xts(ita$Value, as.Date(ita$Date))
  xts.jpn <- xts(jpn$Value, as.Date(jpn$Date))
  xts.net <- xts(net$Value, as.Date(net$Date))
  xts.nor <- xts(nor$Value, as.Date(nor$Date))
  xts.por <- xts(por$Value, as.Date(por$Date))
  xts.spa <- xts(spa$Value, as.Date(spa$Date))
  xts.swe <- xts(swe$Value, as.Date(swe$Date))
  xts.usa <- xts(usa$Value, as.Date(usa$Date))
  xts.uk <- xts(uk$Value, as.Date(uk$Date))
  
  #---Identifies if a data set with the crisis countries e.g. shorter sample should be create or not 
  if (crisis == FALSE){
    data <- merge(xts.aus, xts.aut, xts.bel, xts.can, xts.den, xts.fin,
                  xts.fra, xts.ger, xts.ire, xts.ita, xts.jpn, 
                  xts.net, xts.nor, xts.spa, xts.swe, xts.usa, xts.uk, join  = "inner")
    colnames(data) <- c('AUS', 'AUT', 'BEL', 'CAN', 'DEN', 'FIN', 'FRA', 'GER' , 'IRE', 'ITA', 'JPN', 'NET', 'NOR', 'SPA', 'SWE', 'USA', 'UK')
  } else {
    data <- merge(xts.aus, xts.aut, xts.bel, xts.can, xts.den, xts.fin,
                  xts.fra, xts.ger, xts.gre, xts.ire, xts.ita, xts.jpn, 
                  xts.net, xts.nor, xts.por, xts.spa, xts.swe, xts.usa, xts.uk, join  = "inner")
    colnames(data) <- c('AUS', 'AUT', 'BEL', 'CAN', 'DEN', 'FIN', 'FRA', 'GER', 'GRE', 'IRE', 'ITA', 'JPN', 'NET', 'NOR', 'POR' ,'SPA', 'SWE', 'USA', 'UK')
  }
  return (data)
}

#------Calculates the Bond Price of the Data Yield data Set
#---data: the xts data structure with bond yields from readData
#---RETURN: an xts data structure with all bond prices
calculatePrice <- function(data){
  price <- apply.daily(data, FUN = function(x) {(100/((1+(x/100))^10))})
  return (price)
}

#------Calculates the Return and Volatility of the Country Price Set
#---data: the xts data structure with bond yields or prices
#---RETURN: a list with two xts elements: first is return, second realized volatility
calculateReturnAndVolatility <- function(data){
  # Deletes NA Values from the Price transformation
  data <- na.exclude(data)
  # Calculates the Return
  ret <- ROC(data, n=1)  
  # Calculates the Return shown in Equation
  ret.normalized <- apply.daily(ret,FUN = function(x){x*100})
  # Calculates Realized Volatility
  vola <- apply.daily(ret.normalized, FUN = function(x){x^2})
  # Deletes NA Values
  vola <- na.exclude(vola)
  ret.normalized <- na.exclude(ret.normalized)
  # Returns a List with one Return Data Set and one Volatility Data Set
  return (list(ret.normalized, vola))
}

#------------------------------------------END Part 1------------------------------------------------------------------------- 

#----------------------------Part 2: Estimation of Static and Dynamik Connectedness:------------------------------------------

#-------Calculates the Connectedness Table for a given Forecast Error Variance Matrix-----------------------------------------
#---spilloverTable: matrix with the spillover table from fevd or cholesky Identification
#---RETURN: Connectedness Table with 'To' and 'From' values
CalculateFullTable <- function(spilloverTable){
  #----Read the Names------
  colname <- colnames(spilloverTable)
  rowname <- rownames(spilloverTable)
  #---Define Matrix and set diagonal to 0--------------------------------------------------------------------------------------
  spilloverMatrix <- as.matrix(spilloverTable)
  diag(spilloverMatrix) <- 0
  #---Calculate the Row and Column Sum-----------------------------------------------------------------------------------------
  From <- rowSums(spilloverMatrix)
  To <- colSums(spilloverMatrix)
  #---Total Sum----------------------------------------------------------------------------------------------------------------
  total <- sum(From)
  To <- c(To, total)
  #---Construct the connectedness table (ret)----------------------------------------------------------------------------------
  ret <- as.matrix(spilloverTable)
  colnames(ret) <- colname
  rownames(ret) <- rowname
  ret <- cbind(ret, From)
  ret <- rbind(ret, To)
  return (ret)
}

#------Random Sample of Cholesky Ordering---
#---data: the estimation set for the VAR model
#---n.ahead: forecast horizon for the forecasts, 10 is the default value
#---lag: the lag level for the VAR, 1 is the default value
#---RETURN: the forecast error variance matrix
SotAverageAnalysis <- function(data, n.ahead=10, lag = 1){
  # Estimate the VAR Model
  EstRes.2c <- VAR(data, type="const", p = lag)
  # Get the covariance matrix
  Sigma <- summary(EstRes.2c)$covres
  # Calculates Phi (Impact Multipliers):  
  Phi <- Phi(EstRes.2c, nstep=n.ahead)
  # Return of result of a sample of Cholesky random orderings
  res <- sot_avg_est(Sigma, Phi, ncores=0)
  return (res)
}  

#-----Calculates the GVD Identification----
#---data: the estimation set for the VAR model
#---n.ahead: forecast horizon for the forecasts, 10 is the default value
#---normalize: if the GVD normalized or not, TRUE is default
#---lag: the lag level for the VAR, 1 is the default value
#---RETURN: the forecast error variance matrix
fevd_generalised <- function(data, n.ahead=10, normalize=TRUE, lag = 1) {
  #---------VAR Model estimation-----------
  model <- VAR(data, type="const", p = lag)
  #--------Check if right package is used---
  if (class(model) != "varest") {
    return("The model class is not varest!")
  }
  #print (summary(model))
  #-------Get Variable Names------------------- 
  names <- summary(model)$names
  #-------Estimate the A H-head MA ceofficients
  A <- Phi(model, n.ahead)
  #-------Calculates the Residuals-------------
  epsilon <- residuals(model)
  #-------Calculates the Covariance matrix-----
  #S <- summary(model)$covres
  #cat ('Summary: ', S, '\n')
  Sigma <- t(epsilon)%*%epsilon / (model$obs - model$p*model$K)
  
  
  cat ('Sigma: ', Sigma, '\n')
  #-------Definition of the numerator----------
  gi <- array(0, dim(A))
  #-------Diagonal elements--------------------
  sigmas <- sqrt(diag(Sigma))
  cat ('sigmas: ', sigmas, '\n')
  
  #-------Loop to construct the numerator
  for (j in 1:dim(A)[3]) {
    cat('A - ',j,' : ',A[,,j],'\n')
    gi[,,j] <- t( t( A[,,j]%*%Sigma ) / sqrt(sigmas) )
    cat('Gi - ',j,' : ', gi[,,j], '\n')
  }
  #-------Definition and Loop to construct denumerator
  d <- array(0, dim(A)[c(2,3)])
  for (j in 1:dim(d)[2]) {
    d[,j] <- diag(A[,,j]%*%Sigma%*%t(A[,,j]))
  }
  #------Make the square and sum of numerator-------
  num <- apply(gi^2,1:2,sum)
  cat ('Num: ', num, '\n')
  #------Make the sum of each row thus each variable---
  den <- c(apply(d,1,sum))
  cat ('Den: ', den, '\n')
  #------Calculate GVD --------------------------------
  fevd <- num/den
  #-------Set the names--------------------------------
  #colnames(fevd) <- names
  #rownames(fevd) <- names
  #-------Check if it should be normalized-------------
  if (normalize) {
    cat ('Sum Fevd: ', apply(fevd,1 ,sum))
    return(fevd/apply(fevd, 1, sum)) # hier berechnet man die anteile. 
  } else {
    return(fevd)
  }
}

#-----Method to calculate the Rolling Dynamik estimaton-----------
#---data: the whole return or volatility data set
#---window_size: the window size for the analysis
#---n.ahead: forecast horizon for the analysis
#---method: if Cholesky or GVD shoud be used, GVD is default
#---RETURN: (n,n,T) array, where the row and columns represent the forecast error table for every timepoint and T the whole time dimension in the third dimension
Rolling_Window <- function(data, window_size, n.ahead=10, method = "GFEVD"){
  
  #--------Set the names for the connectedness table--------------
  names <- colnames(data)
  names.length <- length(names)
  names.row <- names
  names.column <- names
  names.row[names.length+1] <- 'To'
  names.column[names.length+1] <- 'From'
  begin <- 1 + window_size
  
  #---------Calculate the number of dynamik estimations-----------
  number.models <- nrow(data)
  n <- ncol(data)
  #----------The Cholesky average analysis gives back min/max and avg values of the connectedness table-------
  avg_array <- array(dim = c(n+1,n+1,number.models))
  rownames(avg_array) <- names.row
  colnames(avg_array) <- names.column
  min_array <- array(dim = c(n+1,n+1,number.models))
  rownames(min_array) <- names.row
  colnames(min_array) <- names.column
  max_array <- array(dim = c(n+1,n+1,number.models))
  rownames(max_array) <- names.row
  colnames(max_array) <- names.column
  
  #---------Read out the dates of the dataset----------------------------------------
  df.date <- data.frame(Start = as.Date(character()), End = as.Date(character())) 
  index.data <- index(data)
  
  #---------Main Loop to calculate the analysis-------------------------------------------
  #---------Iterate from the end to the beginning-----------------------------------------
  for (i in nrow(data):begin){
    #---Calculate both dates of the dynamic analysis--------------------------------------
    j = i-window_size
    #---Get the dates
    start <- index.data[j]
    end <- index.data[i]
    df.date[i,1] <- start
    df.date[i,2] <- end
    #---Check which Identification scheme should be used-----------------------------------
    if (method == "GFEVD"){
      #---Calculate GVD---------------------------------------------------------------------
      result <- fevd_generalised(data[j:i,], n.ahead=n.ahead, normalize=TRUE)
      #---Calculate Full Connectedness Table for each Time Point----------------------------
      avg_array[,,i] <- CalculateFullTable(result)
    } else{
      #---Identification and Connectedness Table calculation for Cholesky Ident.------------
      result <- SotAverageAnalysis(data[j:i,], n.ahead)
      avg_array[,,i] <- CalculateFullTable(result$Average) 
      min_array[,,i] <- CalculateFullTable(result$Minimum)
      max_array[,,i] <- CalculateFullTable(result$Maximum)  
    }
  }
  #---Return the Rolling Dynamic Analysis---------------------------------------------------
  output <- list(df.date, avg_array, min_array, max_array)
  return (output)
}

#-----Static Connectedness Table ---------------------------------#
#---forecast.error.matrix: The static forecast errror variance decomposition matrix
#---RETURN: The Connectedness table
ConnectednessTable <- function(forecast.error.matrix){
  n <- dim(forecast.error.matrix)[1]
  full.matrix <- CalculateFullTable(forecast.error.matrix*100)
  tc <- full.matrix['To','From']/n
  full.matrix['To','From'] <- tc
  return (full.matrix)
}


#-------------------------END Part 2-------------------------------------------------------------------------------------------

#-------------------------Part 3: Print and Slice Method-----------------------------------------------------------------------

#-----Calculates the Connectedness Table for a Result of the Rolling Forecast Error Variance table over the whole T samples-----------
#---spilloverTable: the (n,n,T) array of the Rolling_Window approach
#---RETURN: a (n,n,T) array with the corresponding Connectedness Tables
CalculateRollingConnectednessTable <- function(spilloverTable){
  n <- dim(spilloverTable)[1]
  m <- dim(spilloverTable)[2]
  l <- dim(spilloverTable)[3]
  ret <- array(dim=c(n+1,m+1,l))
  for (i in 1:l){
    spilloverMatrix <- as.matrix(spilloverTable[,,i])
    col.countries <- colnames(spilloverMatrix)
    row.countries <- rownames(spilloverMatrix)
    own.weights <- diag(spilloverMatrix)
    diag(spilloverMatrix) <- 0
    from <- rowSums(spilloverMatrix)
    to <- colSums(spilloverMatrix)
    total <- sum(from)
    to <- c(to, total)
    diag(spilloverMatrix) <- own.weights
    spilloverMatrix <- cbind(spilloverMatrix, from)
    spilloverMatrix <- rbind(spilloverMatrix, to)
    col.names <-c(col.countries, 'From')
    row.names <-c(row.countries, 'To')
    rownames(spilloverMatrix) <- row.names
    colnames(spilloverMatrix) <- col.names
    ret[,,i] <- as.matrix(spilloverMatrix)
  }
  rownames(ret) <- row.names
  colnames(ret) <- col.names
  return (ret)
}

#------Extracts the total Connectedness of a (n,n,T) over T samples-----
#---data: the (n,n,T) data array from CalculateRollingConnectednessTable
#---time: the date object (T) of the Rolling 
#---division: parameter if the totalConnectedness should be divided by other number than the number of countries. Usually not 
totalConnectedness <- function(data, time, division = NULL){
  if (is.null(division) == FALSE){
    countries = division
  } else {
    n = dim(data)[1]
    countries = n-1
  }
  end.date <- time[,2]
  print (head(end.date))
  end.date <- na.exclude(end.date)
  #end.date <- rev(end.date)
  end.date <- as.data.frame(end.date)
  total.connectedness <- data['To', 'From',]
  total.connectedness <- total.connectedness*100
  total.connectedness <- na.exclude(total.connectedness)
  total.connectedness <- as.data.frame(total.connectedness)
  time.series <- xts(total.connectedness, as.Date(end.date$end.date))
  colnames(time.series) <- 'Total Connectedness'
  normalized <- time.series/countries
  colnames(normalized) <- 'Total Normalized Connectedness'
  ret <- merge(time.series, normalized)
  return (ret)
}


#--------Calculates the Quantiles for the Network Plots
quantiles <- function(x){
  return (quantile(x, c(0.25, 0.5, 0.75, 0.9)))
}



#---------Reads out the Network (Connectedness Table) for a specific time point
#---matrix.data: the array (n,n,T) from Rolling Dynamic Analysis
#---matrix.date: the date object (T) of the Rolling Dynamic Analysis
#---time.point: the date of the netwoek in YY-MM-DD format
findNetworkfromTimePoint <- function(matrix.data, matrix.date, time.point){
  sub <- subset(matrix.date, End == time.point)
  index <- rownames(sub)
  index <- as.integer(index)
  connectednessTable <- matrix.data[,,index]
  return (connectednessTable)
}

#---Calculate the internal Connectedness of a group of countries. 
#---data: the full connectedness table of the Rolling Window approach
#---time: the date dataframe of the Rolling Window approach
#---group: is a vector of countries
#---goup.id: is a text for the column header e.g. member of the groups
#---RETURN: the xts dataset with the group connectedness measure
calculate_Group_Connectedness <- function(data, time, group, group.id){
  
  # division=N is the dimension of the full data sample thus for the division through 1/N
  division <- dim(data)[1]
  
  # Selects the sub network and calulates the Connectedness table for it
  euro.countries <- data[group,group,]
  euro.countries <- CalculateRollingConnectednessTable(euro.countries)
  
  # Extracts the total Connectedness over all time periods
  TC.euro.countries <- totalConnectedness(euro.countries, time, division)
  
  # Selects the percentage share of the Time Series
  TC.euro.countries.normalized  <- TC.euro.countries[,'Total.Normalized.Connectedness']
  colnames(TC.euro.countries.normalized) <- 'Normalized Connectedness EURO Countries'
  return (TC.euro.countries.normalized)
}


#----Calculates the Flow between two Groups of Countries-----------------------------------------------------------------------
#---data: the full connectedness table of the Rolling Window approach
#---time: the date dataframe of the Rolling Window approach
#---countries.block.1: a vector of countries
#---countries.block.2: a vector of disjoint countries compared to countries.block.1
#---name.block.i: identifier for the country groups
#---RETURN: an xts dataset that contains both flow connectedness measures
calculate_Flow_Connectedness <- function(data, time, countries.block.1, countries.block.2, name.block.1, name.block.2){
  division <- dim(data)[1]
  
  # Set the transition matrix to with the dimension and data of the full network
  transition <- data[1:19,1:19,]
  
  # Set the structur of the internal countries to 0. 
  transition[countries.block.1, countries.block.1,] <- 0
  transition[countries.block.2, countries.block.2,] <- 0
  
  # First just Copy to new variable, second, the flow connectedness measures. See the example in  the master Thesis
  to.countries.block.2 <- transition
  to.countries.block.2[countries.block.1, countries.block.2,] <- 0
  to.countries.block.1 <- transition
  to.countries.block.1[countries.block.2, countries.block.1,] <- 0
  to.countries.block.1 <- CalculateRollingConnectednessTable(to.countries.block.1)
  to.countries.block.2 <- CalculateRollingConnectednessTable(to.countries.block.2)
  
  # Calculate the total Connectedness of every block
  TC.to.countries.block.1 <- totalConnectedness(to.countries.block.1, time, division)
  TC.to.countries.block.2 <- totalConnectedness(to.countries.block.2, time, division)
  TC.to.countries.block.1 <- TC.to.countries.block.1[,2]
  TC.to.countries.block.2 <- TC.to.countries.block.2[,2]
  
  # Set the labels of the blocks
  label.to.block.1 <- paste('TC from', name.block.2, 'to', name.block.1, sep= " ")
  label.to.block.2 <- paste('TC from', name.block.1, 'to', name.block.2, sep= " ")
  colnames(TC.to.countries.block.1) <- label.to.block.1
  colnames(TC.to.countries.block.2) <- label.to.block.2
  
  # Join the results to one data set and return the results 
  res <- merge(TC.to.countries.block.1, TC.to.countries.block.2)
  return (res)
}


#----Makes a Network Graph with the 75% and 90% 
#---data.ret: The input connectedness table for a specific date
#---name: The title for the network
#---RETURN: NONE, just prints the graph to the Plots area
visualizeGraph <- function(data.ret, name){
  
  #Return Data:
  row <- nrow(data.ret)
  col <- ncol(data.ret)
  matrix.ret <- data.ret[1:(row-1),1:(col-1)]
  
  #Get the right definition like in the master thesis
  matrix.ret <- t(matrix.ret*100)
  
  # Save the diagonal elements to set the size of the nodes
  nodeSize <- diag(matrix.ret)
  
  # Set diag elements to 0 to become a adjacency matrix
  diag(matrix.ret) <- 0
  
  #Create the full return network
  network.ret <- graph.adjacency(matrix.ret, weighted = TRUE, mode="directed")
  E(network.ret)$width <- E(network.ret)$weight/4
  E(network.ret)$arrow.size <- 0.3
  V(network.ret)$size <- 20
  
  #Make the quantiles for the ret
  edges <- E(network.ret)$weight 
  quantiles <- quantiles(edges)
  #cat (quantiles)
  
  network.new <- network.ret
  
  #Set the network edge according to the edge quantile
  E(network.new)$color <- ifelse(E(network.new)$weight > quantiles[[4]], "#252525", ifelse(E(network.new)$weight > quantiles[[3]], "#737373", "#F0F0F0"))
  
  # 75% Quantile Return: 
  cut.off <- quantiles[[3]]
  net.quant.75 <- delete.edges(network.new, E(network.new)[weight < cut.off])
  E(net.quant.75)$arrow.size <- 0.6
  E(net.quant.75)$width <- E(net.quant.75)$weight/20
  plot(net.quant.75, vertex.label.font = 1.3, vertex.label.cex = 1.3, layout=layout.auto, edge.curved = .4, vertex.label.color = 'black')
  tit <- paste('75% Quantiles', name, sep = " | ")
  title(tit)
    
  # 90% Quantile: 
  cut.off <- quantiles[[4]]
  net.quant.90 <- delete.edges(network.new, E(network.new)[weight < cut.off])
  E(net.quant.90)$arrow.size <- 0.6
  E(net.quant.90)$width <- E(net.quant.90)$weight/20
  plot(net.quant.90, vertex.label.font = 1.3, vertex.label.cex = 1.3, layout=layout.auto, edge.curved = .6, vertex.label.color = "black")
  tit <- paste('90% Quantiles', name, sep=" | ")
  title(tit)
}




#-------------------------End Part 3-------------------------------------------------------------------------------------------






#--------------------------Programm Control------------------------------------------------------------------------------------
calculateModel <- function(n.ahead = 10, lag = 1){
  
  #---Get the Data----------------------------------
  data <- readData(crisis = TRUE)
  prices <- calculatePrice(data)
  ret.vola <- calculateReturnAndVolatility(prices)
  ret <- ret.vola[[1]]
  vola <- ret.vola[[2]]
  
  #---Calculate Static Model-------------------------
  ret.table.fevd <- fevd_generalised(ret, n.ahead = n.ahead, lag)
  vola.table.fevd <- fevd_generalised(vola, n.ahead = n.ahead, lag)
  
  #---Calculate Static Spillover Table---------------
  total.connectedness.table.return <- ConnectednessTable(ret.table.fevd) 
  total.connectedness.table.vola <- ConnectednessTable(vola.table.fevd)
  
  #---Calculate Dynamic Model------------------------
  ret.dynamic <- Rolling_Window(ret, 300, 10, 'GFEVD')
  vola.dynamic <- Rolling_Window(vola, 300, 10, 'GFEVD')
  
  #---Calculate Total Connectedness-------------------
  ret.dynamic.total <- totalConnectedness(ret.dynamic[[2]],ret.dynamic[[1]])
  
  vola.dynamic.total <- totalConnectedness(vola.dynamic[[2]], vola.dynamic[[1]])
  
  
  
  
  #---------------------------------------------------EURO - NON EURO --------------------------------------------------------------------------------------------------------------------------
  #---Calculate Return Split Euro - Non Euro------------------
  ret.dynamic.euro <- calculate_Group_Connectedness(ret.dynamic[[2]], ret.dynamic[[1]], c('SPA', 'ITA', 'AUT', 'GER', 'FIN', 'BEL', 'NET', 'FRA', 'IRE', 'GRE', 'POR'), 'Euro')
  ret.dynamic.non.euro <- calculate_Group_Connectedness(ret.dynamic[[2]], ret.dynamic[[1]], c('UK', 'DEN', 'SWE', 'CAN', 'JPN', 'NOR', 'AUS', 'USA'), 'Non Euro')  
  ret.dynamic.euro.split <- merge(ret.dynamic.euro, ret.dynamic.non.euro)
  
  #---Calculate Return Flow Euro - Non Euro--------------------
  ret.euro.flow <- calculate_Flow_Connectedness(ret.dynamic[[2]], ret.dynamic[[1]], c('SPA', 'ITA', 'AUT', 'GER', 'FIN', 'BEL', 'NET', 'FRA', 'IRE', 'GRE', 'POR'), c('UK', 'DEN', 'SWE', 'CAN', 'JPN', 'NOR', 'AUS', 'USA'), 'Euro', 'Non Euro')
  
  
  
  #---Calculate Vola Split Euro - Non Euro------------------
  vola.dynamic.euro <- calculate_Group_Connectedness(vola.dynamic[[2]], vola.dynamic[[1]], c('SPA', 'ITA', 'AUT', 'GER', 'FIN', 'BEL', 'NET', 'FRA', 'IRE', 'GRE', 'POR'), 'Euro')
  vola.dynamic.non.euro <- calculate_Group_Connectedness(vola.dynamic[[2]], vola.dynamic[[1]], c('UK', 'DEN', 'SWE', 'CAN', 'JPN', 'NOR', 'AUS', 'USA'), 'Non Euro')  
  vola.dynamic.euro.split <- merge(vola.dynamic.euro, vola.dynamic.non.euro)
  
  #---Calculate Vola Flow Euro - Non Euro--------------------
  vola.euro.flow <- calculate_Flow_Connectedness(vola.dynamic[[2]], vola.dynamic[[1]], c('SPA', 'ITA', 'AUT', 'GER', 'FIN', 'BEL', 'NET', 'FRA', 'IRE', 'GRE', 'POR'), c('UK', 'DEN', 'SWE', 'CAN', 'JPN', 'NOR', 'AUS', 'USA'), 'Euro', 'Non Euro')
  #----------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  #--------------------------------------------------CORE CRISIS --------------------------------------------------------------------------------------
  
  #---Calculate Return Crisis Split -------------------------
  ret.dynamic.core <- calculate_Group_Connectedness(ret.dynamic[[2]], ret.dynamic[[1]], c('AUT', 'GER', 'FIN', 'BEL', 'NET','FRA'), 'Core')
  ret.dynamic.crisis <- calculate_Group_Connectedness(ret.dynamic[[2]], ret.dynamic[[1]], c('SPA', 'ITA', 'IRE', 'GRE', 'POR'), 'Crisis')
  ret.core.crisis <- merge(ret.dynamic.core, ret.dynamic.crisis)
  
  #---Calculate Return Crisis Flow ---------------------------
  ret.core.crisis.flow <- calculate_Flow_Connectedness(ret.dynamic[[2]], ret.dynamic[[1]], c('AUT', 'GER', 'FIN', 'BEL', 'NET','FRA'), c('SPA', 'ITA', 'IRE', 'GRE', 'POR'), 'Core', 'Crisis')
  
  
  
  #---Calculate Vola Crisis Split -------------------------
  vola.dynamic.core <- calculate_Group_Connectedness(vola.dynamic[[2]], vola.dynamic[[1]], c('AUT', 'GER', 'FIN', 'BEL', 'NET','FRA'), 'Core')
  vola.dynamic.crisis <- calculate_Group_Connectedness(vola.dynamic[[2]], vola.dynamic[[1]], c('SPA', 'ITA', 'IRE', 'GRE', 'POR'), 'Crisis')
  vola.core.crisis <- merge(vola.dynamic.core, vola.dynamic.crisis)
  
  #---Calculate Vola Crisis Flow ---------------------------
  vola.core.crisis.flow <- calculate_Flow_Connectedness(vola.dynamic[[2]], vola.dynamic[[1]], c('AUT', 'GER', 'FIN', 'BEL', 'NET','FRA'), c('SPA', 'ITA', 'IRE', 'GRE', 'POR'), 'Core', 'Crisis')
  #-----------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #--Viszualize Return dynamic Network-------------------------------------------
   ret.snapshot.1 <- findNetworkfromTimePoint(ret.dynamic[[2]], ret.dynamic[[1]], '2009-11-02')
   ret.snapshot.2 <- findNetworkfromTimePoint(ret.dynamic[[2]], ret.dynamic[[1]], '2012-02-03')
   ret.snapshot.3 <- findNetworkfromTimePoint(ret.dynamic[[2]], ret.dynamic[[1]], '2012-10-19')
   ret.snapshot.4 <- findNetworkfromTimePoint(ret.dynamic[[2]], ret.dynamic[[1]], '2015-07-15')
  
   ret.snapshots <- list(ret.snapshot.1, ret.snapshot.2, ret.snapshot.3, ret.snapshot.4)  

  
  
  #--Viszualize Volatility dynamic Network-----------------------------------------
   vola.snapshot.1 <- findNetworkfromTimePoint(vola.dynamic[[2]], vola.dynamic[[1]], '2009-11-02')
   vola.snapshot.2 <- findNetworkfromTimePoint(vola.dynamic[[2]], vola.dynamic[[1]], '2012-02-03')
   vola.snapshot.3 <- findNetworkfromTimePoint(vola.dynamic[[2]], vola.dynamic[[1]], '2012-10-19')
   vola.snapshot.4 <- findNetworkfromTimePoint(vola.dynamic[[2]], vola.dynamic[[1]], '2015-07-15')

   vola.snapshots <- list(vola.snapshot.1, vola.snapshot.2, vola.snapshot.3, vola.snapshot.4)

  
  #Save the results in a list: 
  ret.list <- list(ret.dynamic.total[,2], ret.dynamic.euro.split, ret.euro.flow, ret.core.crisis, ret.core.crisis.flow)
  vola.list <- list(vola.dynamic.total[,2], vola.dynamic.euro.split, vola.euro.flow, vola.core.crisis, vola.core.crisis.flow)
  
  
  return (list(total.connectedness.table.return,total.connectedness.table.vola, ret.list, vola.list, ret.snapshots, vola.snapshots))
}





