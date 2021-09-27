
LL.shift <- function(params, cov.table, names.cov, Bcalc=F, remove.first.transaction = F){
  
  # cov.table = transactions.cov.estimation
  # Bcalc=F
  
  # Gamma parameter:
  # alpha > 0
  # r > 0
  # beta > 0
  # process parameter:
  # 0 <= a < 1
  # b >=0
  
  # cov.table: Contains all nonzero repeaters (zero repeaters are ignored!!!)
  # Id: Id of each customer 
  # date.weekly: Covariate period of their transaction
  # Date: Date of their respective transactions
  # Price: Price of their respective transactions
  
  require(data.table)
  require(lubridate)
  
  alpha <- exp(params["alpha"])
  r     <- exp(params["r"])
  b <- params["b"] # exp for b is not ideal, since then it cannot be exactly 0
  beta <- exp(params["beta"])
  #a <- 1/(1+exp(params["a"])) #find a better restriction for a! 0 <= a < 1
  #a <- params["a"] #we use the bound of optimx
  vec.gamma <- params[5:length(params)]
  
  #fix parameter
  a <- 0
  
  # bring gammas to same order as columns
  vec.gamma  <- vec.gamma[names.cov]
  
  # extract the names
  cols.gamma  <-  paste0("gamma.",names(vec.gamma))
  
  # Add gammas to tables
  cov.table[,  (cols.gamma)  := as.data.table(t(vec.gamma))]
  
  # Remove the first transaction if requested
  if(remove.first.transaction){
    customers.all <- unique(cov.table[, list(Id)])
    cov.table <- cov.table[!(trans.no==1),]
    #adjust the transaction count
    cov.table[, trans.no:=(seq_len(.N)), by=list(Id)]
    }
  
  #get the maximum number of transactions per customer
  cov.table[, max.x:=max(trans.no), by=Id]
  
  #At
  cov.table[,  At := exp(rowSums(cov.table[,  .SD, .SDcols = cols.gamma]  * cov.table[,  .SD, .SDcols = names.cov]))]
  
  #x
  cbs <- cov.table[, list(x=max(max.x)), by=Id]



  # Refactoring starts here ------------------------------------

  #list of all customers
  #customers <- unique(cov.table$Id)
  
  #B1
  #B1 <- cov.table[trans.no==1, At] * (b + cov.table[trans.no==1, Price])
  
  
  # Sort is important as the following calculations depend on the fact
  # that the IDs for transactions = 1 and >1 are in all relevant objects
  # in the same order
  setkey(cov.table, "max.x", "Id")
  setkey(cbs, "x", "Id")
  
  # Split in observations with transactions = 1 and > 1
  B1_maxx_equal_1 <- cov.table[max.x == 1 & trans.no==1, At] * (b + cov.table[max.x == 1 & trans.no==1, Price])
  B1_maxx_equal_2 <- cov.table[max.x == 2 & trans.no==1, At] * (b + cov.table[max.x == 2 & trans.no==1, Price])
  # B1_maxx_larger_2 <- cov.table[max.x > 2 & trans.no==1, At] * (b + cov.table[max.x > 2 & trans.no==1, Price]) # Calculation see below
  xi_maxx_equal_1 <- cbs[x == 1]
  xi_maxx_equal_2 <- cbs[x == 2]
  xi_maxx_larger_2 <- cbs[x > 2]
  customer_number_maxx_equal_1 <- nrow(xi_maxx_equal_1)
  customer_number_maxx_equal_2 <- nrow(xi_maxx_equal_2)
  customer_number_maxx_larger_2 <- nrow(xi_maxx_larger_2)
  
  # ========================================
  # For all customer with single transaction (no loop)
  # ========================================
  
  #B
  # Derivable from previous caculations
  
  # F1
  Aj_maxx_equal_1 <- cov.table[max.x==1, At]
  F1_maxx_equal_1 <- beta * log(Aj_maxx_equal_1) + (beta-1) * log(cov.table[max.x==1, Price] - a * 0 + b)
  
  # RETURN object
  parts_maxx_equal_1 <- data.table(Id=cov.table[max.x==1, Id], F1=F1_maxx_equal_1, B=B1_maxx_equal_1)
  
  
  # parts_maxx_equal_1 <- foreach(i=1:customer_number_maxx_equal_1, .packages = c("data.table"))%dopar%{
  #   
  #   cust.no <- xi_maxx_equal_1[i, Id]
  #   
  #   #xi <- 1
  #   
  #   B <-B1_maxx_equal_1[i]
  #   
  #   #First part of the LL
  #   Aj <- cov.table[Id==cust.no, At]
  #   F1 <- sum(beta * log(Aj) + (beta-1) * log(cov.table[max.x==1, Price] - a * 0 + b))
  # 
  #   return(data.table(Id=cust.no, F1=F1, B=B))
  #   
  #   
  # }
  
  
  # ========================================
  # For all customer with 2 transactions (no loop)
  # ========================================
  
  # B
  Aj_maxx_equal_2_for_B <- cov.table[max.x==2 & trans.no==2, At] # Only for 2nd transaction
  Bi_maxx_equal_2 <- Aj_maxx_equal_2_for_B * (cov.table[max.x==2 & trans.no==2, Price] - a * cov.table[max.x==2 & trans.no==1, Price] + b)
  
  B_maxx_equal_2 <- B1_maxx_equal_2 + Bi_maxx_equal_2
  
  
  # F1
  Aj_maxx_equal_2_for_F1_transaction_1 <- cov.table[max.x==2 & trans.no==1, At] # For 1st transaction
  
  F1_maxx_equal_2_part_1 <- beta * log(Aj_maxx_equal_2_for_F1_transaction_1) + 
    (beta-1) * log(cov.table[max.x==2 & trans.no==1, Price] - a * 0 + b)
  
  F1_maxx_equal_2_part_2 <- beta * log(Aj_maxx_equal_2_for_B) + 
    (beta-1) * log(cov.table[max.x==2 & trans.no==2, Price] - a * cov.table[max.x==2 & trans.no==1, Price] + b)
  
  F1_maxx_equal_2 <- F1_maxx_equal_2_part_1 + F1_maxx_equal_2_part_2
  
  
  # RETURN object
  
  parts_maxx_equal_2 <- data.table(Id=xi_maxx_equal_2[, Id], F1=F1_maxx_equal_2, B=B_maxx_equal_2)
  # cov.table[max.x==2 & trans.no==2, Id]  # Check order of IDs
  
  # ---------------
  
  
  
  # ========================================
  # For all customer with more than 2 transactions (loop)
  # ========================================
  
  
  # ----
  # Get relevant data
  Bi_temp <- cov.table[max.x>2, list(Id, Price, trans.no, max.x, At)]
  
  # CACULATE B-PART ------------------
  
  # Calculate first part of Bi
  #B1_maxx_larger_2_updated <- Bi_temp[trans.no==1, At] * (b + Bi_temp[trans.no==1, Price])
  
  # Lag price variable
  setkey(Bi_temp, Id, trans.no)
  Bi_temp[ , Price_lag:=shift(Price), by=Id]
  
  # Add artifical transaction 0 to make life easier for computation of lag 
  #Bi_temp2 <- rbind(setDT(Bi_temp), Bi_temp[,.SD[.N], Id][, `:=`(trans.no=0, Price=NA, At=NA) ])[order(Id, trans.no)]
  
  # Calculate Bi
  Bi_temp[ , Bi:= At * (Price - a * Price_lag + b)] #trans.no>2
  
  Bi_temp[trans.no==1, Bi:=At * (b + Price)]
  
  B_temp <- Bi_temp[, list(B=sum(Bi)), by=Id]
  
  # OLD:
  #B_temp <- Bi_temp[, list(B=sum(Bi, na.rm=T)), by=Id]
  #B_temp$B_final <- B_temp$B + B1_maxx_larger_2_updated
  
  # Test purposes
  # test <- rbindlist(parts_maxx_larger_2)
  # test2 <- merge(test, B_temp[, list(Id, B_final)], by="Id")
  # all.equal(test2$B, test2$B_final)
  # test2[Id=="10058632040",] # 603 + 156
  # Bi_temp[Id=="10058632040", sum(Bi, na.rm=T)]
  # cov.table[max.x > 2 & trans.no==1]
  # cov.table[max.x > 2 & trans.no==1][order(Id)]
  
  # CACULATE F1-PART ------------------
  
  # Special case for first transaction
  Bi_temp[trans.no==1, F1_temp:=beta * log(At) + (beta-1) * log(Price - a * 0 + b)]
  # For all other transactions (trans.no>1)
  Bi_temp[trans.no>1, F1_temp:=beta * log(At) + (beta-1) * log(Price - a * Price_lag + b)]
  
  F_temp <- Bi_temp[, list(F1=sum(F1_temp)), by=Id]
  
  # For test purposes:
  # test3 <- merge(test, F_temp[, list(Id, F1_final=F1)], by="Id")
  # all.equal(test3$F1, test3$F1_final)
  
  # RETURN object ------------------
  
  parts_maxx_larger_2 <- merge(F_temp, B_temp[, list(Id, B)] , by="Id")


  # # ----
  # 
  # parts_maxx_larger_2 <- foreach(i=1:customer_number_maxx_larger_2, .packages = c("data.table"))%dopar%{
  #   
  #   # i=56
  #   # 760.0342  743.7672
  # 
  #   
  #   cust.no <- xi_maxx_larger_2[i, Id]
  #   
  #   xi <- xi_maxx_larger_2[i, x] # Max transaction per customer
  #   
  #   Bi <- sum(sapply(seq(from = 2, to=xi), 
  #                    function(j){
  #                      # j=1
  #                      Aj<- cov.table[Id==cust.no & trans.no==j, At]
  #                      return(Aj* (cov.table[Id==cust.no & trans.no==j, Price] - a * cov.table[Id==cust.no & trans.no==(j-1), Price] + b)) 
  #                    }))
  #   
  #   B <- B1_maxx_larger_2[i]+ Bi
  #   
  #   
  #   #First part of the LL
  #   F1 <- sum(sapply(seq(from = 1, to=xi), function(j){
  #     
  #     Aj<- cov.table[Id==cust.no & trans.no==j, At]
  #     
  #     return(beta * log(Aj) + (beta-1) * log(cov.table[Id==cust.no & trans.no==j, Price] - a * ifelse(j==1, 0, cov.table[Id==cust.no & trans.no==(j-1), Price]) + b))
  #   }))
  #   
  #   
  #   # RETURN object
  #   return(data.table(Id=cust.no, F1=F1, B=B))
  #   
  #   
  # }
  
  
  # ========================================
  # Combine all 3 parts
  # ========================================
  
  parts <- rbindlist(list(parts_maxx_equal_1, parts_maxx_equal_2, parts_maxx_larger_2))
  setkey(parts, Id)

  # parts_orig <- parts
  # test_parts <- merge(parts_orig, parts, by="Id") 
  # all.equal(test_parts$F1.x, test_parts$F1.y)
  # all.equal(test_parts$B.x, test_parts$B.y)
  # test_parts[test_parts$B.x != test_parts$B.y,]
  # cov.table[Id=="10009126200",]
  
  F1 <- parts$F1
  B <- parts$B
  
  Bmatrix<-parts[, list(Id, B)]
  
  
  #Individual Likelihood
  setkey(cbs, Id)
  LL <- F1 + r*log(alpha) + lgamma(beta*cbs[, x] +r) - (beta*cbs[, x]+r)*log(B + alpha) - lgamma(r) - cbs[, x] * lgamma(beta)
  
  f <- -sum(LL)
  if(is.nan(f)){
    f <- 10^10
  }
  
  
  
  # optimx can only use functions returning scalar arguments
  if (Bcalc==T){
    
    # If first transactions were removed: Add all zero-repeaters again to Bmatrix. They where lost, when removing the first transaction.
    # The value for B in for the missing customers is 0.
    if(remove.first.transaction){
      Bmatrix <- merge(Bmatrix, customers.all, by="Id", all.y=TRUE)
      Bmatrix[is.na(B), B:=0][]
    }
  
    return(list(f=f,Bmatrix=Bmatrix))
  }else{
    return(f)
  }
  
  
}

