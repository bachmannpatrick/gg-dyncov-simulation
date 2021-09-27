
sim.dyncov <- function(params, cov.table, names.cov){
  
  # Gamma parameter:
  # alpha > 0
  # r > 0
  # beta > 0
  # process parameter:
  # 0 <= a < 1, set to 0
  # b >=0
  
  require(data.table)
  require(lubridate)
  #require(doFuture)
  #require(foreach)
  
  alpha <- exp(params["alpha"])
  r     <- exp(params["r"])
  b <- params["b"]
  beta <- exp(params["beta"])
  a <- params["a"]
  vec.gamma <- params[6:length(params)]

  
  # bring gammas to same order as columns
  vec.gamma  <- vec.gamma[names.cov]
  
  # extract the names
  cols.gamma  <-  paste0("gamma.",names(vec.gamma))
  
  # Add gammas to tables
  cov.table[,  (cols.gamma)  := as.data.table(t(vec.gamma))]
  
  #get the maximum number of transactions per customer
  cov.table[, max.x:=max(trans.no), by=Id]
  
  
  #At
  cov.table[,  At := exp(rowSums(cov.table[,  .SD, .SDcols = cols.gamma]  * cov.table[,  .SD, .SDcols = names.cov]))]
  
  #x
  cbs <- cov.table[, list(x=max(max.x)), by=Id]
  
  #list of all customers
  customers <- unique(cov.table$Id)
  
  cov.table[, Price:=NaN]
  
  #loop thought the customers
  for(i in 1:length(customers)){
  #foreach(i=1:length(customers), .packages = c("data.table"))%dopar%{
    cust.no <- customers[i]
    
    #print(i)
    #get max x for the customer
    xi <- unique(cov.table[Id==cust.no, max.x])
    
    zeta <- rgamma(1,shape=r, rate=alpha)
    
    sapply(seq(from = 1, to=xi), function(j){
      zetat <- unname(unlist(zeta*cov.table[Id==cust.no & trans.no==j, "At"]))
      cov.table[Id==cust.no & trans.no==j, Price:=ifelse(j==1, rgamma(1,shape=beta, rate=zetat) - b,
                                                          a*cov.table[Id==cust.no & trans.no==(j-1), "Price"] + rgamma(1,shape=beta, rate=zetat) - b)]
    })
    
  }
  
  
  
  
  return(cov.table)
}

