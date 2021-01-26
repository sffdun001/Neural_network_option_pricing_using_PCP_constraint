#Duncan Saffy 2019 Dissertation
#Option pricing model


rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/Market Data")

Option_suffix = c('H8', "M8", "U8", "Z8", 'H9', "M9", "U9")

library(dplyr)

data_main <- vector(mode = "list", length = 6)

for(i in 1:6){ #not including last 2 for test set
  if (i ==6 ){
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
  
  data_main[[i]]= quarter_data.clean
}
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
  #fl_act <- get(output_activation) #final layer activation
  #act <- get(activation)
  A_C = list() # This will be of lenght L+1 and is a list of the activations evaluated at each layer
  A_P = list()
  A_C[[1]] = X
  A_P[[1]] = X
  Z_C = list()
  Z_P = list()
  
  for (l in 1:L){ #forward propogating inputs
    Z_C[[l]] = A_C[[l]]%*%W_C[[l]]
    Z_P[[l]] = A_P[[l]]%*%W_P[[l]]
    A_C[[l+1]] = cbind(1,sig(Z_C[[l]]))
    A_P[[l+1]] = cbind(1,sig(Z_P[[l]]))
  }
  
  Z_C[[L+1]] = A_C[[L+1]]%*%W_C[[L+1]]
  Z_P[[L+1]] = A_P[[L+1]]%*%W_P[[L+1]]
  chat_star = sig(Z_C[[L+1]]) 
  phat_star = sig(Z_P[[L+1]])
  
  #applying constraint
  chat = chat_star/2 - (exp(-X[,(t+1)]*rf) - phat_star - exp(-X[,(t+1)]*rf)*X[,S+1])/2
  phat = phat_star/2 + (exp(-X[,(t+1)]*rf) + chat_star - exp(-X[,(t+1)]*rf)*X[,S+1])/2
  
  
  cErr = sum((chat-C)^2) / train_size
  pErr = sum((phat-P)^2) / train_size
  
  cStarErr = sum((chat_star-C)^2) / train_size
  pStarErr = sum((phat_star-P)^2) / train_size
  
  
  return(list(chat =chat, phat = phat, chat_star =chat_star, phat_star = phat_star, A_C= A_C, A_P = A_P, CallError = cErr, PutError = pErr,
              X = X, P=P, C=C, Z_C = Z_C, Z_P = Z_P, pStarErr = pStarErr, cStarErr = cStarErr))
}

#function for backpropagation-----------------------------------------------------------------------------------------------------

backprop_HC2 <- function(X, t, S, C, P, W_C, W_P, hidden, learning_rate, activation, 
                     res, L, lambda, train_size, output_activation)
{#res is the output of the forward propagation
 #for this function only W_bar represent weight matricies that include bias weights and W is without bias weights
    W_C_bar = W_C
    W_P_bar = W_P
    for (l in 1:(L+1)){
      W_C[[l]] = W_C[[l]][-1,,drop = F] #drop = F prevents a single column matrix from becoming type = numeric
    }
    for (l in 1:(L+1)){
      W_P[[l]] = W_P[[l]][-1,,drop = F] 
    }
    #d_act = get(paste("d", activation, sep = ""))
    #d_fl_act = get(paste("d", output_activation, sep = ""))
    delta_C = W_C; delta_P = W_P; # Set = W some where initialize a list of correct length
    delta_W_C = W_C; delta_W_P = W_P
    
    delta_C[[L+1]] = ((res$chat_star - res$C))*dsig(res$Z_C[[L+1]])  #(res$chat - res$C) + (res$phat - res$P) +
    delta_W_C[[L+1]] = t(res$A_C[[L+1]])%*%delta_C[[L+1]]
    
    delta_P[[L+1]] = ((res$phat_star - res$P))*dsig(res$Z_P[[L+1]]) 
    delta_W_P[[L+1]] = t(res$A_P[[L+1]])%*%delta_P[[L+1]]
    
    
    for (l in 1:L){
      delta_C[[L+1-l]] = delta_C[[L+2-l]]%*%t(W_C[[L+2-l]])*dsig(res$Z_C[[L+1-l]])
      delta_W_C[[L+1-l]] = t(res$A_C[[L+1-l]])%*%delta_C[[L+1-l]]
      delta_P[[L+1-l]] = delta_P[[L+2-l]]%*%t(W_P[[L+2-l]])*dsig(res$Z_P[[L+1-l]])
      delta_W_P[[L+1-l]] = t(res$A_P[[L+1-l]])%*%delta_P[[L+1-l]]
    }
    
    for (l in 1:(L+1)){
      W_C_bar[[l]] = W_C_bar[[l]] - learning_rate * delta_W_C[[l]] - learning_rate*lambda * rbind(0,W_C[[l]]) 
      W_P_bar[[l]] = W_P_bar[[l]] - learning_rate * delta_W_P[[l]] - learning_rate*lambda * rbind(0,W_P[[l]]) 
    }
    
    Call_Err   = sum((res$chat-C)^2) / train_size
    Put_Err    = sum((res$phat-P)^2) / train_size
    #Parity_Err = sum(((res$chat + exp(-(X[,(rf+1)]) * X[,(t+1)])- res$phat - X[,(S+1)]))^2) / train_size
    
    return(list(W_C= W_C_bar, W_P = W_P_bar, Call_Err = Call_Err, Put_Err = Put_Err))
}

#Neuralnetwork function -------------------------------------------------------------------------------------------------------------

neural_net_HC2 <- function(X, C, P, rf, t, S , hidden, iter_lim, learning_rate, 
                           activation = 'sig', train_size, lambda = 0, output_activation)
{
  X_bar = cbind(1,X) #adding a column of ones to design matrix
  p = train_size
  L = length(hidden)
  q = ncol(X_bar)-1
  W_C = list()
  for (i in 1:(L+1)){ #initialising weights 
    if (i == L+1){
      k = 1
    }
    else{
      k = hidden[i]
    }
    W_C[[i]] = matrix(runif((q+1)*k,-5,5), q+1, k)
    q = k
  }
  
  q = ncol(X_bar)-1
  W_P = list()
  for (i in 1:(L+1)){ #initialising weights 
    if (i == L+1){
      k = 1
    }
    else{
      k = hidden[i]
    }
    W_P[[i]] = matrix(runif((q+1)*k,-1,1), q+1, k)
    q = k
  }
  
  call_error_tracker = numeric(iter_lim) #This is just something so we can track the error after each epoch
  put_error_tracker = numeric(iter_lim) 
  parity_error_tracker = numeric(iter_lim)
  cntr = 0
  PrevErrC = 0 
  PrevErrP = 0 
  for (m in 1:iter_lim) #or some alternative condition
  { 
    random_vec = sample(c(1:nrow(X_bar)), size = train_size*100000, replace = TRUE)
    for (n in 1:100000){
    samp = random_vec[((n-1)*train_size+1):(n*train_size)]
    X_star = X_bar[samp,]
    C_star = C[samp,]
    P_star = P[samp,]
    rf_star = rf[samp]
    
    forward = forwardprop_HC2(X = X_star, t = t, rf = rf_star, C= C_star, P = P_star, S = S, W_C = W_C, W_P = W_P, hidden = hidden,
                               activation = activation, output_activation= output_activation, L= L, train_size = train_size)
    
    back = backprop_HC2(X = X_star, t = t, S = S, C= C_star, P=P_star, W_C = W_C, W_P = W_P, 
                    hidden = hidden, learning_rate = learning_rate, activation = activation, output_activation= output_activation,
                    L = L, res= forward, lambda = lambda, train_size = train_size)
    W_C = back$W_C
    W_P = back$W_P
   
    }
    
    learning_rate = learning_rate
    Fwd = forwardprop_HC2(X = X_bar, C= C, t = t, rf = rf, P = P, S = S, W_C = W_C, W_P =W_P, hidden = hidden,
                    activation = activation, output_activation= output_activation, L= L, train_size = nrow(X_bar))
    
    #PE = sum(((FC$Yhat + exp(-(X_bar[,(rf+1)]) * X_bar[,(t+1)])- FP$Yhat - X_bar[,(S+1)]))^2) / nrow(X_bar)
    
    parity_error_tracker[m] = Fwd$CallError
    call_error_tracker[m] = Fwd$cStarErr
    put_error_tracker[m] = Fwd$pStarErr
    #parity_error_tracker[m] = PE
    print(c(Fwd$CallError,  Fwd$cStarErr, Fwd$pStarErr))
    
    if((m>10)){
      if(((mean(call_error_tracker[(m-10):(m-5)]) -mean(call_error_tracker[(m-4):m]))<1e-6))
      {break}
    }
     
    if(PrevErrC-Fwd$cStarErr<1e-7) {
        cntr = cntr + 1
    }else{
      cntr = 0}
    
    if(cntr == 5){
      break
    }
    PrevErrC = Fwd$cStarErr
    PrevErrP = Fwd$pStarErr
  }
  
  
  return(list(W_C= W_C, W_P = W_P, Call_Error = call_error_tracker, Put_Error = put_error_tracker, 
              parity_error_tracker = parity_error_tracker, Chat = Fwd$chat, Phat = Fwd$phat, C = C, P = P))
}




set.seed(1995)

mlpHC3_market_model <- vector(mode = "list", length = 32)
mlpHC3_market_cv <- vector(mode = "list", length = 32)
periods = c(1,2,3,4,5,6)


hidden_layers = list()
hidden_layers[[1]] = 4;hidden_layers[[2]] = 8;hidden_layers[[3]] = c(4,4);hidden_layers[[4]] = c(8,8);
learning_rates = c(0.1,0.01)
lambdas = c(0.001,0.0001,0.00001,0)

#creating model names
model_names = c()
model_no = 1
for(hidden in 1:4){
  for( lr in learning_rates){
    for(l in lambdas){
      model_names = c(model_names, paste(model_no,"-", "Hidden:", hidden, ",","LearingRate:",lr, ",","Lambda:",l,sep = ""))
      model_no = model_no + 1
    }
  }
}

length(model_names)
CV_res = c()
model_no = 1

StartTime = Sys.time()
for(h in 1:4){
  for(lr in learning_rates){
    for(l in lambdas){
      for (cv_fold in 1:length(data_main)){
        
        train = c()
        trainperiod = periods[-cv_fold]
        for(i in c(1:6)[-cv_fold]){
          train = rbind(train,data_main[[i]])
        }
        test= data_main[[cv_fold]]
        
        
        market_data = train[,c('StockPrice', 'TimeTillExpiry','CallPrice','PutPrice', 'InterestRate')]
        market_data[,c('StockPrice','CallPrice', 'PutPrice')] = train[,c('StockPrice','CallPrice', 'PutPrice')]/ train$Strike
        market_data[, 'TimeTillExpiry'] =  train[,'TimeTillExpiry'] / 365
        market_data[, 'InterestRate'] = log((1+train[,"InterestRate"]/400)^4)
        
        
        m1 <- neural_net_HC2(X = as.matrix(market_data[,c('StockPrice','TimeTillExpiry')]), 
                             C = as.matrix(market_data[, 'CallPrice']), P= as.matrix(market_data[, 'PutPrice']),
                             S = 1, t = 2, hidden = hidden_layers[[h]], iter_lim  = 50, learning_rate = lr, 
                             rf = as.matrix(market_data[, 'InterestRate']),
                             train_size = 100, lambda = l, activation = 'sig', output_activation = 'sig')
        
        
        MSE = mean((m1$Chat-market_data[, 'CallPrice'])^2)
        MAE = mean(abs(m1$Chat-market_data[, 'CallPrice']))
        R2 = 1 - sum((m1$Chat-(market_data[, 'CallPrice']))^2) /
          sum((market_data[, 'CallPrice'] - mean(market_data[, 'CallPrice']))^2)
        print(paste('Train', cv_fold ,c(MSE,MAE,R2), sep = " "))
        
        market_test = test[,c('StockPrice', 'TimeTillExpiry','CallPrice', 'PutPrice')]
        market_test[,c('StockPrice','CallPrice', 'PutPrice')] = test[,c('StockPrice','CallPrice', 'PutPrice')]/ test$Strike
        market_test[, 'TimeTillExpiry'] =  test[,'TimeTillExpiry'] / 365
        market_test[, 'InterestRate'] = log((1+test[,"InterestRate"]/400)^4)
        
        forward = forwardprop_HC2(X = as.matrix(cbind(1,market_test[,c('StockPrice', 'TimeTillExpiry')])), t = 2, 
                                  rf = market_test[, 'InterestRate'], C= market_test[, 'CallPrice'], P = market_test[, 'PutPrice'],
                                  S = 1, W_C = m1$W_C, W_P = m1$W_P, hidden = hidden_layers[[h]],
                                  activation = 'sig', output_activation= 'sig', L= length(hidden_layers[[h]]), train_size = nrow(market_test))
      
        
        CALL_CV_MSE = mean((forward$chat-market_test[, 'CallPrice'])^2)
        CALL_CV_MAE = mean(abs(forward$chat-market_test[, 'CallPrice']))
        CALL_CV_R2 = 1 - sum((forward$chat-(market_test[, 'CallPrice']))^2) /
          sum((market_test[, 'CallPrice'] - mean(market_test[, 'CallPrice']))^2)
        
        PUT_CV_MSE = mean((forward$phat-market_test[, 'PutPrice'])^2)
        PUT_CV_MAE = mean(abs(forward$phat-market_test[, 'PutPrice']))
        PUT_CV_R2 = 1 - sum((forward$phat-(market_test[, 'PutPrice']))^2) /
          sum((market_test[, 'PutPrice'] - mean(market_test[, 'PutPrice']))^2)
      
        
        model_res = c(CALL_CV_MSE, CALL_CV_MAE, CALL_CV_R2, PUT_CV_MSE, PUT_CV_MAE, PUT_CV_R2)
        CV_res= rbind(CV_res, model_res)
        
        mlpHC3_market_model[[model_no]][[i]] = m1
        
      }
      
      CV_Err = colMeans(CV_res)
      CV_res  = rbind(CV_res,CV_Err)
      CV_res = as.data.frame(CV_res)
      names(CV_res) = c("CALL_CV_MSE", "CALL_CV_MAE", "CALL_CV_R2", "PUT_CV_MSE", "PUT_CV_MAE", "PUT_CV_R2")
      mlpHC3_market_cv[[model_no]] = CV_res
      print(paste("#########################################MODEL NUMBER:", model_no, sep = " "))
      model_no = model_no + 1
      print(paste("CV Error:", CV_Err, sep = " "))
      CV_res = c()
    }
  }
}


EndTime = Sys.time()

RunTime = StartTime-EndTime

#setNames(mlpHC3_market_model, model_names)
#setNames(mlpHC3_market_cv, model_names)

save(mlpHC3_market_model, file = paste('mlpHC3_market_001_', Sys.Date(), ".RData", sep = ''))
save(mlpHC3_market_cv, file = paste('mlpHC3_market_cv_001_', Sys.Date(), ".RData", sep = ''))








