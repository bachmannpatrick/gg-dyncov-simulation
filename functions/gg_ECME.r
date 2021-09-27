


conditionallik<- function(a, cov.table, summarytable, params.noa, names.cov, remove.first.transaction=F){ # alpha, r,b, beta, vec.gamma)
  
  
  require(data.table)
  require(lubridate)
  
  ## Step 1: Calculate B (depends on a!!)
  if("B" %in% colnames(summarytable)){
    summarytable[, B := NULL]
  }
  
  

  
  ## Step 0: Extract everything relevant
  alpha <- exp(params.noa["alpha"])
  r     <- exp(params.noa["r"])
  b <- params.noa["b"] # exp for b is not ideal, since then it cannot be exactly 0
  beta <- exp(params.noa["beta"])
  vec.gamma <- params.noa[5:length(params.noa)]
  # bring gammas to same order as columns
  vec.gamma  <- vec.gamma[names.cov]  
  
  #params.noa["alpha"] <-alpha
  #params.noa["r"] <-r
  #params.noa["b"] <-b
  #params.noa["beta"] <-beta
  #params.noa[5:(length(vec.gamma)+5)] <- vec.gamma
  
  
  ## Important: We already remove the first transaction in this function!
  B <- LL.shift_afixed(params.noa, cov.table = cov.table, names.cov=names.cov, remove.first.transaction = remove.first.transaction, Bcalc = T, a=a)$Bmatrix
  summarytable<-merge(summarytable, B)  
  
  
  if(remove.first.transaction){
    customers.all <- unique(cov.table[, list(Id)])
    cov.table <- cov.table[!(trans.no==1),]
    #adjust the transaction count
    cov.table[, trans.no:=(seq_len(.N)), by=list(Id)]
  }
  
  #get customer ids
  customers=summarytable$Id
  
  
  f0<-function(price){
    val<-c()
    #print(length(vec))
    for (i in 1:length(price)){
      if (i==1){
        val[i]= log(b + price[1]) # (beta-1)*log(b + vec[1])
      }else{
        val[i]<- log(b + price[i] - a*price[i-1]) ##(beta-1)*log(b + vec[i] - a*vec[i-1]) ## (beta-1) is multiplied below!
      }
    }
    return(sum(val))
  }
  
  cov.table[, paramdiff:=f0(Price) , by="Id" ]
  
  summarytable<-merge(summarytable, unique(cov.table[, list(paramdiff), by="Id"]), all.x=T)
  
  #print(summarytable)
  
  #loop thought the customers
  #B1
  #B1 <-(beta-1)*log(b + cov.table[trans.no==1, Price])
  #for(i in 1:length(customers)){
  # summarytable$paramdiff <- unlist( foreach(i=1:length(customers), .packages = c("data.table"))%dopar%{
  #     cust.no <- customers[i]
  #     
  #     #get max x for the customer
  #     xi <- unique(cov.table[Id==cust.no, max.x])
  #     
  #     #B
  #     if(xi>1){
  #       Bi<-sum(sapply(seq(from = 2, to=xi), function(j){
  #         Aj<- cov.table[Id==cust.no & trans.no==j, At]
  #         return(beta*log(Aj) + (beta-1)* log(cov.table[Id==cust.no & trans.no==j, Price] - a * cov.table[Id==cust.no & trans.no==(j-1), Price] + b)) 
  #       }))
  #       #B[i] <-B1[i]+ Bi
  #       B <-B1[i]+ Bi
  #       
  #     } else{
  #       #B[i] <-B1[i]
  #       B <-B1[i]
  #     }
  #     
  # 
  #     
  #     return(B)
  #     #return(list(F1=F1, B=B))
  # })
  
  
  
  #summarytable<-merge(summarytable, cov.table[, c("Id","At")])
  
  ## whole formula would be longer than that (see notes, but if we only solve for a, this is enough!!)
  summarytable[, Lhat:=  (beta-1)*paramdiff - B*vhat , by="Id" ]
  
  
  ## Continue here!!
  
  
  return( -mean(summarytable$Lhat)  ) #sum(summarytable$Lhat)
  
}




ggECME<-function(params, cov.table, names.cov, remove.first.transaction=F){

## parameter translation
## old    ## new
## alpha  ## gamma
## r      ## q
## beta   ## p
  
  
  # Clean!!

 
  ## Step 0: Extract everything relevant
  alpha <- exp(params["alpha"])
  r     <- exp(params["r"])
  b <- params["b"] # exp for b is not ideal, since then it cannot be exactly 0
  beta <- exp(params["beta"])
  #a <- 0 #find a better restriction for a! 0 <= a < 1
  a <- params["a"] #we use the bound of optimx
  vec.gamma <- params[6:length(params)]
  # bring gammas to same order as columns
  vec.gamma  <- vec.gamma[names.cov]  
  
  params.noa <- c(params[1:4], params[6:length(params)])

  summarytable <- data.table(Id=unique(cov.table$Id))
  setkey(summarytable, "Id")

  # extract the names
  cols.gamma  <-  paste0("gamma.",names(vec.gamma))
  # Add gammas to tables
  cov.table[,  (cols.gamma)  := as.data.table(t(vec.gamma))]
  
  cov.table[, max.x:=max(trans.no), by=Id]

  ## Use this in case we remove the first transaciton
  cov.table.withfirsttransaction <- cov.table
  
  if(remove.first.transaction){
    customers.all <- unique(cov.table[, list(Id)])
    cov.table <- cov.table[!(trans.no==1),]
    #adjust the transaction count
    cov.table[, trans.no:=(seq_len(.N)), by=list(Id)]
    
    ## cov.table now has its first transaction removed
  }
  
  #x
  cbs <- cov.table[, list(x=max(max.x)), by=Id]
  #cov.table[,  At := exp(rowSums(cov.table[,  .SD, .SDcols = cols.gamma]  * cov.table[,  .SD, .SDcols = names.cov]))]
  
 
  
  setkey(cov.table, "max.x", "Id")
  setkey(cbs, "x", "Id")
  summarytable<-merge(summarytable,cbs[,x, by="Id"])

  ########## Algorithm ---------------
  # start values for the while loop
  LLnew <- 2
  LLold <- 1
  
  #optimize until the results are close enough
  while((LLnew-LLold)/LLold > 1e-5){
  
    ## Step 1: Calculate B ----
    # Achtung: LL.shift removes the first transaction again
    Res <- LL.shift_afixed(params.noa, cov.table = cov.table.withfirsttransaction, names.cov=names.cov, remove.first.transaction = remove.first.transaction, Bcalc = T,a=a)
 
    #"old" LL value
    LLold <- (-Res$f)
    summarytable<-merge(summarytable, Res$B) 
    
    print(paste("Likelihood old is", LLold))
    
    ## E step ----
    summarytable[, vhat:=  (beta*x + r)/(B+ alpha) ]
    
    #needs to be fixed!
    f0<-function(price){
      val<-c()
      for (i in 1:length(price)){
        if(i==1){
          val[i]=Inf
        }else{
      val[i]<-(price[i]+b)/price[i-1]
        }
      }
      return(min(val))
    }
    
    cov.table[, maxval:=NULL]
    cov.table[, maxval:= f0(Price) , by="Id"]
    up<-min(min(cov.table$maxval)-min(cov.table$maxval)/1000,1)
    
    if (a > up){a=up}
    
    
    if (up <= 0){
      
      a<-0
    }else{
    
    ## Optimization 1 ----
    # optimize parameter a
      # Achtung: conditionallik removes the first transaction again
    results.a<- optimx(a, conditionallik, cov.table = cov.table.withfirsttransaction, summarytable = summarytable, params.noa=params.noa, names.cov=names.cov, remove.first.transaction=remove.first.transaction,
                     #method = ("Nelder-Mead"),
                     method = "L-BFGS-B",
                     lower=c(0),
                     upper=c(up),
                     control=list(kkt=F,
                                  maxit=300000))
    #extract parameter a
    a<-unlist(unname(results.a[1]))
    }
    
    #needs to be fixed!
    f1<-function(price){
      val<-c()
      for (i in 1:length(price)){
        if(i==1){
          val[i]=-Inf
        }else{
          val[i] <- (a*price[i-1]) - price[i]
        }
      }
      return(max(val))
    }
    
    cov.table[, minval:=NULL]
    cov.table[, minval:= f1(Price) , by="Id"]
    low<-max(max(cov.table$minval)+max(cov.table$minval)/1000,0)
    
    
    ## Optimization 2 ----
    # optimize all other parameters
    # Achtung: LL.shift removes the first transaction again
    results<- optimx(params.noa, LL.shift_afixed, cov.table = cov.table.withfirsttransaction, names.cov = names.cov, remove.first.transaction = remove.first.transaction, a=a,
                     #method = ("Nelder-Mead"),
                     method = "L-BFGS-B",
                     lower=c(log(1*10^(-5)),log(1*10^(-5)),max(low,0),log(1*10^(-5)), rep(-5, length(params.noa)-4)),
                     upper=c(log(200),log(200),200,log(200), rep(5, length(params.noa)-4)),
                     control=list(kkt =F,
                                  maxit=300000)) # , trace=6
    
    # ## For a=0, low=0, the above should be the same as this:
    # results<- optimx(params, LL.shift, cov.table = transactions.cov.estimation, names.cov = names.cov, remove.first.transaction = T,
    #                  #method = ("Nelder-Mead"),
    #                  method = "L-BGFS-B",
    #                  lower=c(log(1*10^(-5)),log(1*10^(-5)),0,log(1*10^(-5)), rep(-3, length(params)-4)),
    #                  upper=c(log(200),log(200),200,log(200), rep(3, length(params)-4)),
    #                  control=list(trace =6,
    #                               maxit=300000))
    
    
    
    params.optim <- (unlist(results[1:(length(names.cov)+4)])) ## logparams!!
    LLnew <- (-results$value)
    
    print(paste("Likelihood new is", LLnew))
    
    alpha<-exp(params.optim["alpha"])
    r<-exp(params.optim["r"])
    b<-params.optim["b"]
    beta<-exp(params.optim["beta"])
    vec.gamma <- params.optim[5:length(params.optim)]
    # bring gammas to same order as columns
    vec.gamma  <- vec.gamma[names.cov]
    
    params.noa <- params.optim
    
  }
  
  #add the a value from the first optimization the the result
  params <- c(params.optim[1:4], a=a, params.optim[5:length(params.optim)])
  
  return(params)

}









