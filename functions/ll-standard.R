
LL.standard <- function(params, cov.table, cbs=NULL){
  
  # Gamma parameter:
  # alpha > 0
  # r > 0
  # beta > 0
  # process parameter:
  # 0 <= a < 1
  # b >=0
  
  require(data.table)
  require(lubridate)
  require(doFuture)
  require(foreach)
  
  alpha <- exp(params["alpha"])
  r     <- exp(params["r"])
  b <- 0 #exp(params["b"])
  beta <- exp(params["beta"])
  a <- 0 #1/(1+exp(params["a"]))
  vec.gamma <- c("MAIL_COV_1"=0,"HIGH_SEASON"=0,"CHANNEL"=0) #params[6:length(params)]
  
  
  #fix parameter
  #a <- 0
  
  #names of covariates
  names.cov <- c("MAIL_COV_1", "HIGH_SEASON","CHANNEL")
  
  # bring gammas to same order as columns
  vec.gamma  <- vec.gamma[names.cov]
  
  # extract the names
  cols.gamma  <-  paste0("gamma.",names(vec.gamma))
  
  # Add gammas to tables
  cov.table[,  (cols.gamma)  := as.data.table(t(vec.gamma))]
  
  
  #get the maximum number of transactions per customer
  cov.table[, max.x:=max(trans.no), by=Id]
  
  # Ignore zero repeaters!
  cov.table<- cov.table[!(max.x==0), ]
  
  
  
  #At
  cov.table[,  At := exp(rowSums(cov.table[,  .SD, .SDcols = cols.gamma]  * cov.table[,  .SD, .SDcols = names.cov]))]
  
  #x
  if (is.null(cbs)){
  cbs <- cov.table[, list(x=max(max.x)), by=Id]
  }
  
  #list of all customers
  customers <- unique(cov.table$Id)
  
  #B1
  B1 <- cov.table[trans.no==1, At] * (b + cov.table[trans.no==1, Price])
  
  #loop thought the customers
  #for(i in 1:length(customers)){
  parts <- foreach(i=1:length(customers), .packages = c("data.table"))%dopar%{
    cust.no <- customers[i]
    
    #get max x for the customer
    xi <- unique(cov.table[Id==cust.no, max.x])
    
    #B
    if(xi>1){
      Bi<-sum(sapply(seq(from = 2, to=xi), function(j){
        Aj<- cov.table[Id==cust.no & trans.no==j, At]
        return(Aj* (cov.table[Id==cust.no & trans.no==j, Price] - a*cov.table[Id==cust.no & trans.no==(j-1), Price] + b)) 
      }))
      #B[i] <-B1[i]+ Bi
      B <-B1[i]+ Bi
      
    } else{
      #B[i] <-B1[i]
      B <-B1[i]
    }
    
    #First part of the LL
    F1 <- sum(sapply(seq(from = 1, to=xi), function(j){
      
      Aj<- cov.table[Id==cust.no & trans.no==j, At]
      
      return(beta * log(Aj) + (beta-1) * log(cov.table[Id==cust.no & trans.no==(j), Price] - a * ifelse(j==1, 0, cov.table[Id==cust.no & trans.no==(j-1), Price]) + b))
    }))
    
    return(data.table(F1=F1, B=B))
    #return(list(F1=F1, B=B))
  }
  
  
  
  parts <- rbindlist(parts)
  F1 <- parts$F1
  B <- parts$B
  
  #Individual Likelihood
  LL <- F1 + r*log(alpha) + lgamma(beta*cbs[, x] +r) - (beta*cbs[, x]+r)*log(B + alpha) - lgamma(r) - cbs[, x] * lgamma(beta)
  
  f <- -sum(LL)
  if(is.nan(f)){
    f <- 10^10
  }
  
  return(f)
}

