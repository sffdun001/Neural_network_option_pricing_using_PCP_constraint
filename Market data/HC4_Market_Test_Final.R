
rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/Market Data")
load('mlpHC4_market_final_2020-09-09.RData')

model = HC4_final

###############################################################################

Option_suffix = c('H8', "M8", "U8", "Z8", 'H9', "M9", "U9")

library(dplyr)

data_main <- vector(mode = "list", length = 7)

for(i in 1:7){ #not including last 2 for test set
  if (i ==6|i ==7){
    quarter_data <- read.csv(paste("AI", Option_suffix[i],
                                   "_processed_data_2019-10-07.csv",
                                   sep = ""), header = T)
  }else{
    quarter_data <- read.csv(paste("AI", Option_suffix[i],
                                   "_processed_data_2019-07-22.csv",
                                   sep = ""), header = T)
  }
  
  
  
  quarter_data.clean  = quarter_data%>%filter(!is.na(CALL_PRICE) | !is.na(PUT_PRICE)) #removing prices that have not year been created
  
  #setting prices to 0 where there were NA values
  na_call_ind = which(is.na(quarter_data.clean$CALL_PRICE))
  na_put_ind = which(is.na(quarter_data.clean$PUT_PRICE))
  
  quarter_data.clean$CALL_PRICE[na_call_ind] = 0
  quarter_data.clean$PUT_PRICE[na_put_ind] = 0
  
  
  quarter_data.clean <- quarter_data.clean%>%
    select(c("CALL_PRICE", "PUT_PRICE", "Strike", "DAY_TILL_EXPIRY",
             "AI1_LAST", "TOP40_LAST", "JIBA3M", "TOP40_EQY_DVD_YLD_12M"))
  names(quarter_data.clean) <- c("CallPrice", "PutPrice", "Strike", "TimeTillExpiry",
                                 "StockPrice","UnderlyingPrice", "InterestRate", "DivYield")
  
  quarter_data.clean$InterestRate = log((1+quarter_data.clean$InterestRate/400)^4)
  
  
  data_main[[i]]= quarter_data.clean
}


options_price_table = vector(mode = 'list',length = 7)
for (i in 1:7){
  data = data_main[[i]]
  strikes_vector = sort(unique(data$Strike))
  price_table = matrix(NA, nrow = 90, ncol = length(strikes_vector))
  for (k in 1:length(strikes_vector)){
    d1 = data[which(data$Strike== strikes_vector[k]),]
    price_table[1:nrow(d1),k] = d1$PutPrice
  }
  tab1 = as.data.frame(price_table)
  colnames(tab1) = strikes_vector
  options_price_table[[i]] = price_table
}

#View( options_price_table[[7]])
Rows_removing = vector(mode = "list", length =  7)

for (i in 1:7){
  data = data_main[[i]]
  Rows_removing[[i]] = matrix(NA, nrow = 2000, ncol = 4)
  l = 1
  for (j in 1:nrow(data)){
    if(data$CallPrice[j] == 0){
      #print(paste(i, ":",j, sep = ""))
      if(sum(data$CallPrice[which((data$Strike == data$Strike[j]) & (data$TimeTillExpiry >=data$TimeTillExpiry[j]))]) == 0 ){
        info1 = c(j, data$StockPrice[j], data$Strike[j], data$TimeTillExpiry[j])
        Rows_removing[[i]][l,] =  info1
        print(paste(i, ":",j, sep = ""))
        l = l+1
      }
    }
    
  }
}

###############################################################################
#head(quarter_data.clean)
#length(na.omit(Rows_removing[[7]][,1]))

#We remove this data point as neven had call price
data_main[[7]] = data_main[[7]][-na.omit(Rows_removing[[7]][,1]),]
data_main[[3]] = data_main[[3]][-na.omit(Rows_removing[[3]][,1]),]


library(dplyr)
#Defining activation functions and derivatives
sig <- function(x)
{
  return(1/(1+ exp(-x)))
}

dsig <- function(x)
{
  return(sig(x)*(1-sig(x)))
}



#Function for fowrard propagation ---------------------------------------------------------------------------------------------------

forwardprop_HC2 <- function(X, S, t, rf, C, P, W_C, W_P, hidden, activation = "sig", L, train_size, output_activation)
{ 
  fl_act <- get(output_activation) #final layer activation
  act <- get(activation)
  A_C = list() # This will be of lenght L+1 and is a list of the activations evaluated at each layer
  A_P = list()
  A_C[[1]] = X
  A_P[[1]] = X
  Z_C = list()
  Z_P = list()
  
  for (l in 1:L){ #forward propogating inputs
    Z_C[[l]] = A_C[[l]]%*%W_C[[l]]
    Z_P[[l]] = A_P[[l]]%*%W_P[[l]]
    A_C[[l+1]] = cbind(1,act(Z_C[[l]]))
    A_P[[l+1]] = cbind(1,act(Z_P[[l]]))
  }
  
  Z_C[[L+1]] = A_C[[L+1]]%*%W_C[[L+1]]
  Z_P[[L+1]] = A_P[[L+1]]%*%W_P[[L+1]]
  chat_star = fl_act(Z_C[[L+1]]) 
  phat_star = fl_act(Z_P[[L+1]])
  
  #applying constraint
  chat = chat_star/2 - (exp(-X[,(t+1)]*rf) - phat_star - X[,S+1])/2
  phat = phat_star/2 + (exp(-X[,(t+1)]*rf) + chat_star - X[,S+1])/2
  
  
  cErr = sum((chat-C)^2) / train_size
  pErr = sum((phat-P)^2) / train_size
  
  cStarErr = sum((chat_star-C)^2) / train_size
  pStarErr = sum((phat_star-P)^2) / train_size
  
  
  return(list(chat =chat, phat = phat, chat_star =chat_star, phat_star = phat_star, A_C= A_C, A_P = A_P, CallError = cErr, PutError = pErr,
              X = X, P=P, C=C, Z_C = Z_C, Z_P = Z_P, pStarErr = pStarErr, cStarErr = cStarErr))
}

Blacks <- function(F0, K, t, r, sd, call = TRUE, DT = 1/365) #So this apparently works
{ 
  d1 <- 1/(sd*sqrt(t*DT)) * (log(F0/K) +(((sd^2)/2)*(t*DT)))
  d2 <- d1 - sd*sqrt(t*DT)
  C <- exp(-r*t*DT)*(F0*pnorm(d1) - K*pnorm(d2))
  if(call){
    return(C) #returning the call price
  }
  else {
    return(C + exp(-r*t*DT)*K - exp(-r*t*DT)*F0) #returning the put price
  }
}


Derivative_HC2 <- function(X, rf, S, C, P, W_C, W_P, hidden, activation, 
                           res, L, output_activation)
{
  
  #X = x.1; S = 1; C = time_seq; P = time_seq; res = test.1
  #W_C = m1$W_C; W_P = m1$W_P; hidden = c(4); activation = "sig"
  #L =1; output_activation = 'sig'
  
  W_C_bar = W_C
  W_P_bar = W_P
  
  for (l in 1:(L+1)){
    W_C[[l]] = W_C[[l]][-1,,drop = F] #drop = F prevents a single column matrix from becoming type = numeric
  }
  
  for (l in 1:(L+1)){
    W_P[[l]] = W_P[[l]][-1,,drop = F] 
  }
  
  d_act = get(paste("d", activation, sep = ""))
  d_fl_act = get(paste("d", output_activation, sep = ""))
  delta_C = W_C; delta_P = W_P; # Set = W some where initialize a list of correct length
  
  
  for (l in 1:L){
    if(l == 1){
      delta_C[[L+1-l]] = (matrix(W_C[[L+2-l]],nrow = nrow(res$Z_C[[L+2-l]]),
                                 ncol = hidden[L], byrow = T))*
        matrix(d_fl_act(res$Z_C[[L+2-l]]), nrow = length(d_fl_act(res$Z_C[[L+2-l]])),
               ncol = hidden[L], byrow = F) 
      
      delta_P[[L+1-l]] = (matrix(W_P[[L+2-l]],nrow = nrow(res$Z_P[[L+2-l]]),
                                 ncol = hidden[L], byrow = T))*
        matrix(d_fl_act(res$Z_P[[L+2-l]]), nrow = length(d_fl_act(res$Z_P[[L+2-l]])),
               ncol = hidden[L], byrow = F) 
    }
    else{
      delta_C[[L+1-l]] = delta_C[[L+2-l]]%*%t(W_C[[L+2-l]])*d_act(res$Z_C[[L+2-l]])
      delta_P[[L+1-l]] = delta_P[[L+2-l]]%*%t(W_P[[L+2-l]])*d_act(res$Z_P[[L+2-l]])
    }
  }
  
  call_star_delta = (d_act(res$Z_C[[1]])*delta_C[[1]])%*%(W_C[[1]][S,])
  #(W_C[[1]][S,]) #(delta_C[[1]])
  put_star_delta = (d_act(res$Z_P[[1]])*delta_P[[1]])%*%(W_P[[1]][S,])
  
  call_derivative = call_star_delta/2 +put_star_delta/2 + 1/2
  put_derivative = call_star_delta/2 +put_star_delta/2 - 1/2
  
  return(list(call_derivative = call_derivative, put_derivative = put_derivative))
}


CallDelta = function(X, Weights, res, hidden, S, activation, 
                     output_activation)
{ #This function finds the derivative of option with respect to the underlying. 
  #X = as.matrix(hutch_data[,c('StockPrice','TimeTillExpiry')])
  #Weights = m1$W ; res = F.Delta; hidden = c(4); S = 1
  #activation = 'sig'; output_activation= 'sig'
  W = Weights
  L = length(hidden)
  W_bar = W
  for (l in 1:(L+1)){
    W[[l]] = W[[l]][-1,,drop = F] #drop = F prevents a single column matrix from becoming type = numeric
  }
  
  d_act = get(paste("d", activation, sep = ""))
  d_fl_act = get(paste("d", output_activation, sep = ""))
  delta = d_fl_act(res$Z[[L+1]])
  #delta_W[[L+1]] = t(res$A[[L+1]])%*%delta[[L+1]]
  for (l in 1:L){
    delta = delta%*%t(W[[L+2-l]])*d_act(res$Z[[L+1-l]])
    #delta_W[[L+1-l]] = t(res$A[[L+1-l]])%*%delta[[L+1-l]]
  }
  
  Delta = delta%*%t(W[[1]][S,,drop = F])
  
  return(list(Delta= Delta))
}



Blacks <- function(F0, K, t, r, sd, call = TRUE, DT = 1/365) #So this apparently works
{ 
  d1 <- 1/(sd*sqrt(t*DT)) * (log(F0/K) +(((sd^2)/2)*(t*DT)))
  d2 <- d1 - sd*sqrt(t*DT)
  C <- exp(-r*t*DT)*(F0*pnorm(d1) - K*pnorm(d2))
  if(call){
    return(C) #returning the call price
  }
  else {
    return(C + exp(-r*t*DT)*K - exp(-r*t*DT)*F0) #returning the put price
  }
}

OptFunction = function(sigma, F0, K, C, P, t, r)
{
  if (C>=1){ans = abs(C- Blacks(F0,  K,  t, r, sd =sigma))}
  else{ans = abs(P- Blacks(F0,  K,  t, r, sd =sigma, call = FALSE))}
  return(ans)
}

#calculating Implied vol



data_main1 = data_main
for(i in 1:7){
  t1 = data_main[[i]]
  for(j in 1:nrow(t1)){
    res = optim(par = 0.173, fn = OptFunction, F0 = t1[j,6], K = t1[j,3],
                C = t1[j,1], P = t1[j,2], t = t1[j,4], r = t1[j,7]) 
    data_main1[[i]]$Implied_vol[j] = res$par
  }
  print(i)
}


vol_means = c()
#Pricing Error
for (i in 1:6){
  vol_means = c(vol_means, mean(data_main1[[i]]$Implied_vol))
}

mean(vol_means)

#17,69291/ 0.1866791
#Hedging Error

MLPMarketHedgeError = vector(mode = 'list', length = 7)
BSMarketHedgeError = vector(mode = 'list', length = 7)
BSImpMarketHedgeError = vector(mode = 'list', length = 7)


for (i in 1:7){
  #i = 1
  print(i)
  Quarter_data = data_main1[[i]]
  k_vec = unique(Quarter_data$Strike)
  k_vec = sort(k_vec, decreasing = FALSE)
  
  N = length(k_vec)
  V_S = matrix(NA, nrow = 90, ncol = N)
  V_C = matrix(NA, nrow = 90, ncol = N)
  V_B = matrix(NA, nrow = 90, ncol = N)
  V = matrix(NA, nrow = 90, ncol = N)
  
  V_bs_S = matrix(NA, nrow = 90, ncol = N)
  #V_bs_C = matrix(NA, nrow = 63, ncol = N)
  V_bs_B = matrix(NA, nrow = 90, ncol = N)
  V_bs = matrix(NA, nrow = 90, ncol = N)
  
  V_bsImp_S = matrix(NA, nrow = 90, ncol = N)
  #V_bs_C = matrix(NA, nrow = 63, ncol = N)
  V_bsImp_B = matrix(NA, nrow = 90, ncol = N)
  V_bsImp = matrix(NA, nrow = 90, ncol = N)
  
  
  
  for (k in 1:N){
    #k=1
    K = k_vec[k]
    
    data = Quarter_data%>%filter(Strike == K)%>%
      arrange(desc(TimeTillExpiry))
    #head(Pathdata)
    #model <-  mlpHC2_black_model###################################
    
    S0 = data$StockPrice[1]
    n = nrow(data)
    
    V_S[n,k] = 0
    V_C[n,k] = -data$CallPrice[1]
    V_B[n,k] = -(V_S[n,k] +V_C[n,k])
    V[n,k]   = V_S[n,k] + V_C[n,k] + V_B[n,k] 
    
    V_bs_S[n,k] = 0 
    V_bs_B[n,k] = -(V_bs_S[n,k]+ V_C[n,k]) 
    V_bs[n,k]   = V_bs_S[n,k] + V_C[n,k] + V_bs_B[n,k] 
    
    V_bsImp_S[n,k] = 0 #S0*BSDelta
    V_bsImp_B[n,k] = -(V_bsImp_S[n,k]+ V_C[n,k]) ### mexp(-0.0687153*TimeToExpiry/365) 
    V_bsImp[n,k]   = V_bsImp_S[n,k]+ V_C[n,k] + V_bsImp_B[n,k] ### *exp(-0.0687153*TimeToExpiry/365) 
    
    data$StockPrice = data$StockPrice / K
    data$TimeTillExpiry = data$TimeTillExpiry / 365
    
    
    d1 = 1/(0.1866791*sqrt(data$TimeTillExpiry[1]))*(log(S0/K) +
                                                       (((0.1866791^2)/2)*(data$TimeTillExpiry[1])))
    
    d1Imp = 1/(data$Implied_vol[1]*sqrt(data$TimeTillExpiry[1]))*(log(S0/K) +
                                                                    (((data$Implied_vol[1]^2)/2)*(data$TimeTillExpiry[1])))
    
    BSImpDelta = pnorm(d1Imp)*exp(-data$InterestRate[1]* data$TimeTillExpiry[1])
    
    #m1 = forwardprop(X =(as.matrix(cbind(1, S0/K, data$TimeTillExpiry[1]))[,,drop = F]), 
    #                 Y= as.matrix(1), 
    #                 W = model$W , hidden = c(4,4), activation = 'sig', 
    #                 output_activation= 'sig', L= 2, train_size = 1)
    m1 = forwardprop_HC2(X = as.matrix(cbind(1,data[1,c('StockPrice','TimeTillExpiry')]))[,,drop = F], 
                         C = as.matrix(1), P  = as.matrix(1),
                         W_C = model$W_C , W_P = model$W_P , hidden = c(4,4), activation = 'sig', 
                         output_activation= 'sig', L= 2, train_size = 1,
                         S =1, t =2, rf = data$InterestRate[1])
    
    BSDelta = pnorm(d1)*exp(-data$InterestRate[1]* data$TimeTillExpiry[1])
    
    
    #MLPDelta = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/365))[,,drop = F]), 
    #                     Weights = model$W, hidden = c(4,4), res = m1, activation = 'sig', 
    #                     output_activation= 'sig', S =1)
    
    MLPDelta = Derivative_HC2(X =as.matrix((data[1,c('StockPrice','TimeTillExpiry')]))[,,drop = F], 
                              C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                              W_C = model$W_C, hidden = c(4,4), res = m1, activation = 'sig', 
                              output_activation= 'sig', S =1, L = 2)
    
    
    for(t in 1:(n-1)){
      #t = 1
      MLPDelta_prev = MLPDelta$call_derivative*(m1$chat>0)
      BSDelta_prev = BSDelta
      BSImpDelta_prev = BSImpDelta
      
      m1 = forwardprop_HC2(X = as.matrix(cbind(1,data[(t+1),c('StockPrice','TimeTillExpiry')]))[,,drop = F], 
                           C = as.matrix(data[,'CallPrice'])[t+1], P  = as.matrix(data[,'PutPrice'])[t+1],
                           W_C = model$W_C , W_P = model$W_P , hidden = c(4,4), activation = 'sig', 
                           output_activation= 'sig', L= 2, train_size = nrow(data),
                           S =1, t =2, rf = data$InterestRate[t+1])
      
      #m1 = forwardprop(X =(as.matrix(cbind(1,(data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
      #                 Y= as.matrix(1), 
      #                 W = model$W , hidden = c(4,4), activation = 'sig', 
      #                 output_activation= 'sig', L= 2, train_size = nrow(data))
      
      #MLPDelta = CallDelta(X = (as.matrix(cbind((data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
      #                     Weights = model$W, hidden = c(4,4), res = m1, activation = 'sig', 
      #                     output_activation= 'sig', S =1)
      
      MLPDelta = Derivative_HC2(X = as.matrix(data[(t+1),c('StockPrice','TimeTillExpiry')])[,,drop = F], 
                                C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                W_C = model$W_C, hidden = c(4,4), res = m1, activation = 'sig', 
                                output_activation= 'sig', S =1, L = 2)
      
      
      
      #r = 0.0687153, sd = 0.173
      V_S[(n-t),k] = 0#data$StockPrice[t+1] * K * MLPDelta$Delta
      V_C[(n-t),k] = -data$CallPrice[t+1]
      V_B[(n-t),k] = V_B[(n-t+1),k]*exp( data$InterestRate[t+1]/365) + MLPDelta_prev * K * (data$StockPrice[t+1] - data$StockPrice[t])
      V[(n-t),k] = V_B[(n-t),k] + V_S[(n-t),k] + V_C[(n-t),k]
      
      d1 =  1/(0.1866791*sqrt(data$TimeTillExpiry[t+1]))*(log(data$StockPrice[t+1]) +
                                                            (((0.1866791^2)/2)*(data$TimeTillExpiry[t+1])))
      
      BSDelta = pnorm(d1)*exp(- data$InterestRate[t+1]*data$TimeTillExpiry[t+1])
      V_bs_S[(n-t),k] =  0 #data$StockPrice[t+1] * K * BSDelta
      V_bs_B[(n-t),k] =  V_bs_B[(n-t+1),k]*exp( data$InterestRate[t+1]/365)  +  BSDelta_prev * K * (data$StockPrice[t+1] - data$StockPrice[t])
      V_bs[(n-t),k] = V_bs_B[(n-t),k] + V_bs_S[(n-t),k] + V_C[(n-t),k]
      
      d1Imp =  1/(data$Implied_vol[t+1]*sqrt(data$TimeTillExpiry[t+1]))*(log(data$StockPrice[t+1]) +
                                                                           (((data$Implied_vol[t+1]^2)/2)*(data$TimeTillExpiry[t+1])))
      
      BSImpDelta = pnorm(d1Imp)*exp(- data$InterestRate[t+1]*data$TimeTillExpiry[t+1])
      V_bsImp_S[(n-t),k] =  0 #data$StockPrice[t+1] * K * BSDelta
      V_bsImp_B[(n-t),k] =  V_bsImp_B[(n-t+1),k]*exp( data$InterestRate[t+1]/365)  +  BSImpDelta_prev * K * (data$StockPrice[t+1] - data$StockPrice[t])
      V_bsImp[(n-t),k] = V_bsImp_B[(n-t),k] + V_bsImp_S[(n-t),k] + V_C[(n-t),k]
      
    }
  }
  
  MLPMarketHedgeError[[i]] = as.data.frame(V)
  colnames(MLPMarketHedgeError[[i]]) = k_vec
  BSMarketHedgeError[[i]]= as.data.frame(V_bs)
  colnames(BSMarketHedgeError[[i]]) = k_vec
  BSImpMarketHedgeError[[i]]= as.data.frame(V_bsImp)
  colnames(BSImpMarketHedgeError[[i]]) = k_vec
  
}


#####
Moneyness_List  = vector(mode = 'list', length = 7)
for (quarter in c(1:7)){
  #quarter = 1
  data = data_main1[[quarter]]
  Strikes = sort(unique(data$Strike))
  a1= c(rep(NA,length(Strikes)))
  Moneyness= c(rep(NA,length(Strikes)))
  for (Strike in 1:length(Strikes)){
    d1 = data[which(data$Strike == Strikes[Strike]),]
    a1[Strike] = d1[which(d1$TimeTillExpiry==max(d1$TimeTillExpiry)),]$TimeTillExpiry
    Moneyness[Strike] = d1$StockPrice[which(d1$TimeTillExpiry ==a1[Strike])]/d1$Strike[which(d1$TimeTillExpiry ==a1[Strike])]
  }
  
  Moneyness_List[[quarter]][[1]] =  which(Moneyness>1.03)
  Moneyness_List[[quarter]][[2]] =  which((Moneyness<1.03)&(Moneyness>0.97))
  Moneyness_List[[quarter]][[3]] =  which((Moneyness<0.97))
}


#################################################
##Pricing Error
#################################################


MLPPriceError = vector(mode = 'list', length = 7)
MLPPriceError_Black = vector(mode = 'list', length = 7)

for(i in 1:7){
  #i =3
  data <-  data_main[[i]]
  k_vec = sort(unique(data$Strike))
  pe = matrix(NA, nrow = length(k_vec), ncol = 3)
  pe_bs = matrix(NA, length(k_vec), ncol = 3)
  
  for (j in 1:length(k_vec)){
    #j = 66
    data1 = data[which(data$Strike==k_vec[j]),]
    data1$StockPrice = data1$StockPrice /data1$Strike
    data1$TimeTillExpiry = data1$TimeTillExpiry / 365
    data1$CallPrice = data1$CallPrice /data1$Strike
    
    m1 = forwardprop_HC2(X = as.matrix(cbind(1,data1[,c('StockPrice','TimeTillExpiry')])), 
                         C = as.matrix(data1[,'CallPrice']), P  = as.matrix(data1[,'PutPrice']),
                         W_C = model$W_C , W_P = model$W_P , hidden = c(4,4), activation = 'sig', 
                         output_activation= 'sig', L= 2, train_size = nrow(data1),
                         S =1, t =2, rf = data1$InterestRate)
    
    #m1 = forwardprop(X = as.matrix(cbind(1,data1[,c('StockPrice','TimeTillExpiry')])), 
    #                 Y= as.matrix(data1[,'CallPrice']), 
    #                 W = model$W , hidden = c(4,4), activation = 'sig', 
    #                 output_activation= 'sig', L= 2, train_size = nrow(data))
    BS = Blacks(F0 = data1$StockPrice*data1$Strike,K = data1$Strike,t =  data1$TimeTillExpiry * 365,
                sd = 0.1866791,r = data1$InterestRate)
    
    BSPrice = BS/data1$Strike
    
    MSE = mean((BSPrice-data1$CallPrice)^2)
    MAE = mean(abs(BSPrice-data1$CallPrice))
    R2 = 1 - sum((BSPrice-data1$CallPrice)^2) / 
      sum((data1$CallPrice - mean(data1$CallPrice))^2)
    pe_bs[j,] = c(MSE,MAE,R2)
    
    MSE = mean((m1$chat*(m1$chat>0)-data1$CallPrice)^2)
    MAE = mean(abs(m1$chat*(m1$chat>0)-data1$CallPrice))
    R2 = 1 - sum((m1$chat*(m1$chat>0)-data1$CallPrice)^2) / 
      sum((data1$CallPrice - mean(data1$CallPrice))^2)
    pe[j,] = c(MSE,MAE,R2)
  }
  MLPPriceError[[i]] = pe
  MLPPriceError_Black[[i]] = pe_bs
}

#Calculating R-squared values
Market_MLP_Price_table = matrix(NA, nrow = 8, ncol = 3)
Market_BS_Price_table = matrix(NA, nrow = 8, ncol = 3)
Market_MLP_Price_table1 = matrix(NA, nrow = 8, ncol = 3)
Market_BS_Price_table1 = matrix(NA, nrow = 8, ncol = 3)

for (i in 1:6){
  for (k in 1:3){
    ind =  Moneyness_List[[i]][[k]]
    Market_MLP_Price_table[i,k] = paste(round(mean(MLPPriceError[[i]][ind,3])*100,2), '%', sep = "")
    Market_BS_Price_table[i,k] = paste(round(mean(MLPPriceError_Black[[i]][ind,3])*100,2), '%', sep = "")
    Market_MLP_Price_table1[i,k] = round(mean(MLPPriceError[[i]][ind,3])*100,2)
    Market_BS_Price_table1[i,k] = round(mean(MLPPriceError_Black[[i]][ind,3])*100,2)
    if(i ==6){
      ind =  Moneyness_List[[7]][[k]]
      Market_MLP_Price_table[8,k] = paste(round(mean(MLPPriceError[[7]][ind,3])*100,2), '%', sep = "")
      Market_BS_Price_table[8,k] = paste(round(mean(MLPPriceError_Black[[7]][ind,3])*100,2), '%', sep = "")
      
      Market_MLP_Price_table[7,k] = paste(round(mean(Market_MLP_Price_table1[c(1:6),k]),2), '%', sep = "")
      Market_BS_Price_table[7,k] = paste(round(mean(Market_BS_Price_table1[c(1:6),k]),2), '%', sep = "")
    }
  }
  
}

t1 = Market_MLP_Price_table
t2 = Market_BS_Price_table

row.names(t1) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')
row.names(t2) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')


write.csv(t1, file = paste('mlpHC4_market_price_error_', Sys.Date(), '.csv', sep = ''))
#write.csv(t2, file = paste('black_market_price_error_', Sys.Date(), '.csv', sep = ''))

#Checking proporttion of insample and out-of-sample
d2 = vector(mode = 'list', length = 7)
for (quarter in c(1:7)){
  data <-  data_main[[quarter]]
  Strikes = sort(unique(data$Strike))
  d2[[quarter]] = matrix(NA, nrow = length(Strikes), ncol = 3)
  for (Strike in 1:length(Strikes)){
    d1 = data[which(data$Strike == Strikes[Strike]),]
    a1 = d1[which(d1$TimeTillExpiry==max(d1$TimeTillExpiry)),]
    d2[[quarter]][Strike,1] = a1$TimeTillExpiry
    d2[[quarter]][Strike,2] = a1$Strike
    d2[[quarter]][Strike,3] = a1$StockPrice/a1$Strike
  }
}

head(d2[[1]])
Moneyness_table = matrix(NA, nrow = 7, ncol = 4)
for (quarter in c(1:7)){
  Moneyness_table[quarter,] =  c( nrow(d2[[quarter]]), 
                                  length(which(d2[[quarter]][,3]>1.03)),
                                  length(which((d2[[quarter]][,3]<1.03)&(d2[[quarter]][,3]>0.97))), 
                                  length(which((d2[[quarter]][,3]<0.97))) 
  )
}
head(data_main1[[1]])
colSums(Moneyness_table)
write.csv(as.data.frame(Moneyness_table), file = "MoneynessTable.csv")

#Identifying columns that have strikes for in-the-money, near-the-money and out-the-money options


#length(Moneyness_List[[4]][[3]])

#################################################
##Hedging Error 
#################################################

Market_MLP_TabV = matrix(NA, nrow = 8, ncol = 4)
Market_MLP_Imp_TabV = matrix(NA, nrow = 8, ncol = 4)
Market_MLP_TabV1 = matrix(NA, nrow = 8, ncol = 4)
Market_MLP_Imp_TabV1 = matrix(NA, nrow = 8, ncol = 4)
for (i in 1:6){
  for (k in 1:3){
    ind =  Moneyness_List[[i]][[k]]
    Market_MLP_TabV[i,k] = round(length(which(abs(MLPMarketHedgeError[[i]][1,ind]) < 
                                                abs(BSMarketHedgeError[[i]][1,ind])))/length(ind), 4)
    Market_MLP_Imp_TabV[i,k] = round(length(which(abs(MLPMarketHedgeError[[i]][1,ind]) < 
                                                    abs(BSImpMarketHedgeError[[i]][1,ind])))/length(ind), 4)
    Market_MLP_TabV1[i,k] = length(which(abs(MLPMarketHedgeError[[i]][1,ind]) < 
                                           abs(BSMarketHedgeError[[i]][1,ind])))/length(ind)
    Market_MLP_Imp_TabV1[i,k] = length(which(abs(BSMarketHedgeError[[i]][1,ind]) < 
                                               abs(BSImpMarketHedgeError[[i]][1,ind])))/length(ind)
    if(i ==6){
      ind =  Moneyness_List[[7]][[k]]
      Market_MLP_TabV[8,k] = round(length(which(abs(MLPMarketHedgeError[[7]][1,ind]) < 
                                                  abs(BSMarketHedgeError[[7]][1,ind])))/length(ind),4)
      Market_MLP_Imp_TabV[8,k] = round(length(which(abs(MLPMarketHedgeError[[7]][1,ind]) < 
                                                      abs(BSImpMarketHedgeError[[7]][1,ind])))/length(ind),4)
      
      
      Market_MLP_TabV[7,k] = round(mean(Market_MLP_TabV1[c(1:6),k]),4)
      Market_MLP_Imp_TabV[7,k] = round(mean(Market_MLP_Imp_TabV1[c(1:6),k]),4)
    }
  }
  Market_MLP_TabV[i,4] = round(length(which(abs(MLPMarketHedgeError[[i]][1,]) < 
                                              abs(BSMarketHedgeError[[i]][1,])))/ncol(MLPMarketHedgeError[[i]]), 4)
  Market_MLP_Imp_TabV[i,4] = round(length(which(abs(MLPMarketHedgeError[[i]][1,]) < 
                                                  abs(BSImpMarketHedgeError[[i]][1,])))/ncol(MLPMarketHedgeError[[i]]), 4)
  if(i ==6){
    Market_MLP_TabV[8,4] = round(length(which(abs(MLPMarketHedgeError[[7]][1,]) < 
                                                abs(BSMarketHedgeError[[7]][1,])))/ncol(MLPMarketHedgeError[[7]]), 4)
    Market_MLP_Imp_TabV[8,4] =  round(length(which(abs(MLPMarketHedgeError[[7]][1,]) < 
                                                     abs(BSImpMarketHedgeError[[7]][1,])))/ncol(MLPMarketHedgeError[[7]]), 4)
    
    
    Market_MLP_TabV[7,4] = round(mean(Market_MLP_TabV[c(1:6),4]),4)
    Market_MLP_Imp_TabV[7,4] = round(mean(Market_MLP_Imp_TabV[c(1:6),4]),4)
  }
  
}

t1 = Market_MLP_TabV
t2 = Market_MLP_Imp_TabV

row.names(t1) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')
row.names(t2) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')
write.csv(t1, file = paste('tab_V_HC4_Market_', Sys.Date(),'.csv', sep = ''))
write.csv(t2, file = paste('tab_V_HC4_Imp_Market_', Sys.Date(),'.csv', sep = ''))

#Tab_VI
Market_MLP_TabVI = matrix(NA, nrow = 8, ncol = 4)
Market_Black_TabVI = matrix(NA, nrow = 8, ncol = 4)
Market_BlackImp_TabVI = matrix(NA, nrow = 8, ncol = 4)
for (i in 1:6){
  for (k in 1:3){
    ind =  Moneyness_List[[i]][[k]]
    Market_MLP_TabVI[i,k] = round(sqrt(mean(as.numeric(MLPMarketHedgeError[[i]][1,ind])^2) + 
                                         var(as.numeric(MLPMarketHedgeError[[i]][1,ind]))), 2)
    Market_Black_TabVI[i,k] = round(sqrt(mean(as.numeric(BSMarketHedgeError[[i]][1,ind])^2) + 
                                           var(as.numeric(BSMarketHedgeError[[i]][1,ind]))), 2)
    Market_BlackImp_TabVI[i,k] =round(sqrt(mean(as.numeric(BSImpMarketHedgeError[[i]][1,ind])^2) + 
                                             var(as.numeric(BSImpMarketHedgeError[[i]][1,ind]))), 2)
    
    if(i ==6){
      ind =  Moneyness_List[[7]][[k]]
      Market_MLP_TabVI[8,k] = round(sqrt(mean(as.numeric(MLPMarketHedgeError[[7]][1,ind])^2) + 
                                           var(as.numeric(MLPMarketHedgeError[[7]][1,ind]))), 2)
      Market_Black_TabVI[8,k] = round(sqrt(mean(as.numeric(BSMarketHedgeError[[7]][1,ind]^2)) +
                                             var(as.numeric(BSMarketHedgeError[[7]][1,ind]))), 2)
      Market_BlackImp_TabVI[8,k] =round(sqrt(mean(as.numeric(BSImpMarketHedgeError[[7]][1,ind])^2) + 
                                               var(as.numeric(BSImpMarketHedgeError[[7]][1,ind]))), 2)
      
      Market_MLP_TabVI[7,k] = mean(Market_MLP_TabVI[c(1:6),k])
      Market_Black_TabVI[7,k] = mean(Market_Black_TabVI[c(1:6),k])
      Market_BlackImp_TabVI[7,k] =mean(Market_BlackImp_TabVI[c(1:6),k])
    }
  }
  
  Market_MLP_TabVI[i,4] = round(sqrt(mean(as.numeric(MLPMarketHedgeError[[i]][1,])^2) + 
                                       var(as.numeric(MLPMarketHedgeError[[i]][1,]))), 2)
  Market_Black_TabVI[i,4] = round(sqrt(mean(as.numeric(BSMarketHedgeError[[i]][1,])^2) + 
                                         var(as.numeric(BSMarketHedgeError[[i]][1,]))), 2)
  Market_BlackImp_TabVI[i,4] = round(sqrt(mean(as.numeric(BSImpMarketHedgeError[[i]][1,])^2) + 
                                            var(as.numeric(BSImpMarketHedgeError[[i]][1,]))), 2)
  
  if(i ==6){
    Market_MLP_TabVI[8,4] = round(sqrt(mean(as.numeric(MLPMarketHedgeError[[7]][1,])^2) + 
                                         var(as.numeric(MLPMarketHedgeError[[7]][1,]))), 2)
    Market_Black_TabVI[8,4] = round(sqrt(mean(as.numeric(BSMarketHedgeError[[7]][1,])^2) + 
                                           var(as.numeric(BSMarketHedgeError[[7]][1,]))), 2)
    Market_BlackImp_TabVI[8,4] = round(sqrt(mean(as.numeric(BSImpMarketHedgeError[[7]][1,])^2) + 
                                              var(as.numeric(BSImpMarketHedgeError[[7]][1,]))), 2)
    
    
    Market_MLP_TabVI[7,4] = mean(Market_MLP_TabVI[c(1:6),4])
    Market_Black_TabVI[7,4] = mean(Market_Black_TabVI[c(1:6),4])
    Market_BlackImp_TabVI[7,4] =mean(Market_BlackImp_TabVI[c(1:6),4])
  }
  
}

t1 = Market_MLP_TabVI
t2 = Market_Black_TabVI
t3 = Market_BlackImp_TabVI

row.names(t1) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')
row.names(t2) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')
row.names(t3) = c('AIH8','AIM8','AIU8','AIZ8','AIH9','AIM9','Train','AIU9/Test')


write.csv(t1, file = paste('tab_VI_HC4_Market_', Sys.Date(),'.csv', sep = ''))
#write.csv(t2, file = paste('tab_VI_Black_Market_', Sys.Date(),'.csv', sep = ''))
#write.csv(t3, file = paste('tab_VI_BlackImp_Market_', Sys.Date(),'.csv', sep = ''))
