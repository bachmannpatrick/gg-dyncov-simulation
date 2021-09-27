####
## Scenario
## XX
####

# set working directory
#setwd("~/pnbd-dyncov-simulation")

# Clean Workspace
#rm(list= ls())

#Preparation ----
#Load libaries
library(CLVTools)
library(data.table)
#library(ggplot2)
library(lubridate)
library(optimx)
library(reshape2)
library(foreach)
library(doParallel)
library(forecast)

# source the functions that do all the magic
  #pnbd-simulation
  source("functions/pnbd_Simulation_Dynamic.R")
  #gg-simulation
  source("functions/sim-dyncov.R")
  #gg-model
  source("functions/ll-shifted.r")
  source('functions/ll-shifted_afixed.r')
  source("functions/ll-standard.r")
  #source("functions/expectation-shifted.r")
  #source("functions/condexpectation-shifted.r")
  #source("functions/egg_plotSpending.r")
  #source("functions/egg_plotExpectation.r")
  #source("functions/egg_getMAETimeSeries.r")
  #source("functions/egg_ExpectationTable.r")
  #source("functions/egg_DECTtable.r")
  #source("functions/pnbd_dyncov_DECT_egg.r")
  #source("functions/egg_getCLV.r")
  #source("functions/egg_pvalue.r")
  #source("functions/egg_getMAECondExpectation.r")
  #source("functions/egg_getMAEExpectation.r")
  #source("functions/egg_getMAERealTransactions.r")
  source("functions/gg_ECME.r")
  #source("functions/gg_estimateModel.r")
  #source("functions/gg_expectation.r")




#Measure computation time
start.time <- Sys.time()

# Options ----
  #Set a seed to ensure consistent results
  set.seed(12345)
  
  #lubridate options
  #set week start to Monday ubridate.week.start = 1
  #set week start to Sunday ubridate.week.start = 0
  options(lubridate.week.start = 0)
  
  # Setting number
  setting <- "001"
  
  # Number of simulations
  sim.no <- 100
  
  # Number of customers to be simulated
  n = 2500 # max: 11000 customers due to exogenous contextual factors
  
  # Duration of estimation period in weeks
  estimation.duration <- 104
  
  # Specify the type of contextual factors (continuous/binary)
  cov.type = "continuous" #"continuous" or "binary"
  
  ## GG:Base parameters ----
  alpha.gg = 2
  r.gg = 4.5
  b.gg = 0
  beta.gg =0.8
  a.gg = 0.0
  ## GG:Covariate parameters ----
  gamma1 = 0
  
  #Gamma-Gamma: remove first transaction
  #Should we consider the first transaction for the spending model?
  #TRUE: No, we do not use it (function default if not defined)
  #FALSE: Yes, we use it
  remove.first.transaction = FALSE
  
  
  
  #available covariates (defined exogenous): "direct.marketing", "high.season", "low.season", "gender"
  covariate.names.gg <- c("direct.marketing") #"high.season", "low.season", "gender"
  
  # Create an empty list to store the results
  results.param <- list()
  
  # Loop through simulation ----
  no.seed <- 1
  no <- 1
  
  while(no <= sim.no){
    ##Load Data ----
    # Load the exogenous transactions and covariates (they are not generated!)
    # Currently there is data available for 2500 customers
    load("data/SimDataCovariates.RData")
    #only use relevant covariates
    covariates.dynamic <- covariates.dynamic[, list(Id, Cov.Date, high.season, low.season, direct.marketing, gender)]
    #fix the Date variable
    mydata[, Date:=ymd(Date)]
    #make sure direct marketing is a dummy variable
    if(cov.type == "binary"){
      covariates.dynamic[direct.marketing>0, direct.marketing:=1]
    }
    
    
    
    print(paste("Starting simulation round", no, "with seed", no.seed, "; Setting:", setting))
    # Change the seed
    set.seed(no.seed)
    
    # Extract all required info form the CLVTools object
    # start-date of the estimation period
    estimation.start <- min(mydata$Date)
    # end-date of the estimation period
    estimation.end <- max(mydata$Date)

    #combine the parameters
    params.true <- c("alpha"= log(alpha.gg),
                     "r"=log(r.gg),
                     "b" = b.gg,
                     "beta" = log(beta.gg),
                     "a"= a.gg,
                     "direct.marketing" = gamma1
                     #high.season" = high.season.gg,
                     #"low.season" = low.season.gg,
                     #"gender"= gender.gg
                     )
    
    # get the weekly date for the transactions in order to allocate the correct weekly covariates
    mydata[, date.weekly:=floor_date(Date, unit = "week")]
    
    # create a table combining the transactions with the weekly covariates
    transactions.cov <- merge(mydata, covariates.dynamic, by.x =c("Id","date.weekly"), by.y = c("Id", "Cov.Date"), all.x = T)
    
    #Add a transaction counter for every customer
    transactions.cov[, trans.no:=(seq_len(.N)), by=list(Id)]
    
    #simulate the prices in cov.table.dyncov
    cov.table.dyncov<-sim.dyncov(params=params.true, cov.table=transactions.cov, names.cov=covariate.names.gg)
    
    
    #We need to add the simulated prices to the transactional data in the CLVTools objects for the two PNBD models
    mydata <- cov.table.dyncov[, list(Id, Date, Price)]
    
    #remove all info we do not need
    transactions.cov <- copy(cov.table.dyncov[,list(Id, date.weekly, Date, Price, high.season, direct.marketing, low.season, gender)])
    
    #Add a transaction counter for every customer
    transactions.cov[, trans.no:=(seq_len(.N)), by=list(Id)]
    
    # eGG: Estimation ----
    #Start parameters
    params1 <- c("alpha"= log(alpha.gg),
                 "r"= log(r.gg),
                 "b" = b.gg,
                 "beta" = log(beta.gg),
                 "a" = a.gg,  #Note: Has to be defined for ECME,
                 "direct.marketing"= gamma1
                 #"high.season"= 0.1,
                 #"gender"=0.1
                 #"low.season"= 0.1,
    )
    
    LL.shift(params=params1, cov.table = transactions.cov, names.cov=covariate.names.gg, remove.first.transaction = remove.first.transaction)
    
    # ECME Algorithm ----
    resultECME <- ggECME(params=params1, cov.table = transactions.cov, names.cov=covariate.names.gg, remove.first.transaction = remove.first.transaction)
    params <- resultECME
    
    # Store parameters
    results.param[[no]] <- data.table(alpha=exp(params[1]), r=exp(params[2]), b=params[3], beta=exp(params[4]), a=params[5], gamma1=params[6])
    
    save.image(file=paste0("simdata/",setting,"-SimData-",cov.type,"-",no,"-", gamma1,".RData"))
    
    #Increase counter
    no.seed <- no.seed+1
    no <- no+1
  }
  
  #combine the model results in data.tables
  results.param.df <- rbindlist(results.param)
  
  save("results.param.df",  file=paste0("results/",setting,"-SimResults-",cov.type, gamma1,".RData"))
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
