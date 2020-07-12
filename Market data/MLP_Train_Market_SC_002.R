#Duncan Saffy 2019 Dissertation
#Option pricing model

rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/Market Data")

Option_suffix = c('H8', "M8", "U8", "Z8", 'H9', "M9", "U9")

library(dplyr)

data_main <- vector(mode = "list", length = 5)

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

forwardprop_soft <- function(X,Y, W, hidden, activation = "sig", L, train_size, output_activation)
{ #This function performs the forward propogation
  #X is the input data (each row is one observation)
  #Y is the output data for which to train the model to
  #W is a list of weights - each element in the list is a set of weights for a layer in the NN. NOTE IT INCLUDES WEIGHTS FOR ONES
  #res1#W = matrix(runif(16), 4,4); Y = 1; X = runif(4);activation = "sig"; hidden = 1
  #fl_act <- get(output_activation) #final layer activation
  #act <- get(activation)
  A = list() # This will be of lenght L+1 and is a list of the activations evaluated at each layer
  A[[1]] = X
  Z = list()
  for (l in 1:L){
    Z[[l]] = A[[l]]%*%W[[l]]
    A[[l+1]] = cbind(1,sig(Z[[l]]))
  }
  Z[[L+1]] = A[[L+1]]%*%W[[L+1]]
  Yhat = sig(Z[[L+1]]) #linear activation in final layer
  Err = sum((Yhat-Y)^2) / train_size
  
  return(list(Yhat = Yhat, A = A, Err = Err, X = X, Y = Y, Z = Z))
}

#function for backpropagation-----------------------------------------------------------------------------------------------------

backprop_soft <- function(X, rf, t, S, C, P, W_C, W_P, hidden, iter_lim, learning_rate, activation, 
                     res_call, res_put , L, lambda, train_size, output_activation)
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
    delta_C[[L+1]] = (res_call$Yhat - res_call$Y)*dsig(res_call$Z[[L+1]]) + 
                     (res_call$Yhat + exp(-rf * X[,(t+1)])- res_put$Yhat - exp(- rf * X[,(t+1)])*X[,(S+1)])*dsig(res_call$Z[[L+1]])# derivative wr.r.t weights connection L -> output layer
    delta_W_C[[L+1]] = t(res_call$A[[L+1]])%*%delta_C[[L+1]]
    delta_P[[L+1]] = (res_put$Yhat - res_put$Y)*dsig(res_put$Z[[L+1]]) - 
                     (res_call$Yhat + exp(- rf * X[,(t+1)])- res_put$Yhat - exp(- rf * X[,(t+1)])*X[,(S+1)])*dsig(res_put$Z[[L+1]])
    delta_W_P[[L+1]] = t(res_put$A[[L+1]])%*%delta_P[[L+1]]
    for (l in 1:L){
      delta_C[[L+1-l]] = delta_C[[L+2-l]]%*%t(W_C[[L+2-l]])*dsig(res_call$Z[[L+1-l]])
      delta_W_C[[L+1-l]] = t(res_call$A[[L+1-l]])%*%delta_C[[L+1-l]]
      delta_P[[L+1-l]] = delta_P[[L+2-l]]%*%t(W_P[[L+2-l]])*dsig(res_put$Z[[L+1-l]])
      delta_W_P[[L+1-l]] = t(res_put$A[[L+1-l]])%*%delta_P[[L+1-l]]
    }
    
    for (l in 1:(L+1)){
      W_C_bar[[l]] = W_C_bar[[l]] - learning_rate * delta_W_C[[l]] - learning_rate*lambda * rbind(0,W_C[[l]]) 
      W_P_bar[[l]] = W_P_bar[[l]] - learning_rate * delta_W_P[[l]] - learning_rate*lambda * rbind(0,W_P[[l]]) 
    }
    
    Call_Err   = sum((res_call$Yhat-C)^2) / train_size
    Put_Err    = sum((res_put$Yhat-P)^2) / train_size
    Parity_Err = sum(((res_call$Yhat + exp(-rf * X[,(t+1)])- res_put$Yhat - exp(-rf * X[,(t+1)])*X[,(S+1)]))^2) / train_size
    
    return(list(W_C= W_C_bar, W_P = W_P_bar, Call_Err = Call_Err, Put_Err = Put_Err, Parity_Err = Parity_Err))
}



#Neuralnetwork function -------------------------------------------------------------------------------------------------------------

neural_net_soft <- function(X, C, P, rf, t, S , hidden, iter_lim, learning_rate, activation = 'sig', train_size, lambda = 0, output_activation)
{ # This function uses stoc_NN to update weights
  # epochs is the number of iterations the of backpropagation we will do
  # learning_rate 
  # sample_size is the number of X values selected for trainging at each epoch
  #rf, t and S are the indexes for the risk free, time and stock price in the input matrix X.
  #activation = 'sig'
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
    W_C[[i]] = matrix(runif((q+1)*k,-1,1), q+1, k)
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
    #check_put = numeric(1000)
    #check_call = numeric(1000)
    #check_parity = numeric(1000)
    for (n in 1:100000){
    samp = random_vec[((n-1)*train_size+1):(n*train_size)]
    X_star = X_bar[samp,]
    C_star = C[samp,]
    P_star = P[samp,]
    rf_star = rf[samp]
    forward_call = forwardprop_soft(X = X_star, Y= C_star, W = W_C, hidden = hidden,
                               activation = activation, output_activation= output_activation, L= L, train_size = train_size)
    forward_put = forwardprop_soft(X = X_star, Y= P_star, W = W_P, hidden = hidden, 
                              activation = activation, output_activation= output_activation, L= L, train_size = train_size)
    back = backprop_soft(X = X_star,  rf = rf_star, t = t, S = S, C= C_star, P=P_star, W_C = W_C, W_P = W_P, 
                    hidden = hidden, learning_rate = learning_rate, activation = activation, output_activation= output_activation,
                    L = L, res_call  = forward_call, res_put = forward_put, lambda = lambda, train_size = train_size)
    W_C = back$W_C
    W_P = back$W_P
    #check_put[n] = back$Put_Err 
    #check_call[n] = back$Call_Err 
    #check_parity[n] = back$Parity_Err 
    
    }
    
    FC = forwardprop_soft(X = X_bar, Y= C, W = W_C, hidden = hidden,
                    activation = activation, output_activation= output_activation, L= L, train_size = nrow(X_bar))
    FP = forwardprop_soft(X = X_bar, Y= P, W = W_P, hidden = hidden, 
                    activation = activation, output_activation= output_activation, L= L, train_size = nrow(X_bar))
    PE = sum(((FC$Yhat + exp(- rf * X_bar[,(t+1)])- FP$Yhat - X_bar[,(S+1)]*exp(- rf * X_bar[,(t+1)])))^2) / nrow(X_bar)
    
    call_error_tracker[m] = FC$Err
    put_error_tracker[m] = FP$Err
    parity_error_tracker[m] = PE
    
    learning_rate = learning_rate #* 0.9
    
    
    print(c(FC$Err, FP$Err, PE))
    
    if((m>10)){
      if(((mean(call_error_tracker[(m-10):(m-5)]) -mean(call_error_tracker[(m-4):m]))<1e-6))
      {break}
    }
    
    if(((PrevErrC - FC$Err) < 1e-7)&((PrevErrP - FP$Err) < 1e-7)){
      cntr = cntr + 1
    }else{
      cntr = 0}
    
    if(cntr == 5){
      break
    }
    PrevErrC = FC$Err
    PrevErrP = FP$Err
  }
  
  
  return(list(W_C= W_C, W_P = W_P, Call_Error = call_error_tracker, Put_Error = put_error_tracker, Parity_Error = parity_error_tracker, 
              Chat = FC$Yhat, Phat = FP$Yhat, C = C, P = P))
}




#Running NN-------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------




#


set.seed(1995)

mlpSC_market_model <- vector(mode = "list", length = 32)
mlpSC_market_cv <- vector(mode = "list", length = 32)
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
        market_data[,c('StockPrice','CallPrice')] = train[,c('StockPrice','CallPrice', 'PutPrice')]/ train$Strike
        market_data[, 'TimeTillExpiry'] =  train[,'TimeTillExpiry'] / 365
        market_data[, 'InterestRate'] = log((1+train[,"InterestRate"]/400)^4)
        
        m1 <- neural_net_soft(X = as.matrix(market_data[,c('StockPrice','TimeTillExpiry')]), 
                        C = as.matrix(market_data[, 'CallPrice']), P= as.matrix(market_data[, 'PutPrice']),
                        S = 1, t = 2, hidden = hidden_layers[[h]], iter_lim  = 50, learning_rate = lr, 
                        rf = as.matrix(market_data[, 'InterestRate']),
                        train_size = 100, lambda = l, activation = 'sig', output_activation = 'sig')
        
        
        MSE = mean((m1$Chat-market_data[, 'CallPrice'])^2)
        MAE = mean(abs(m1$Chat-market_data[, 'CallPrice']))
        R2 = 1 - sum((m1$Chat-(market_data[, 'CallPrice']))^2) /
          sum((market_data[, 'CallPrice'] - mean(market_data[, 'CallPrice']))^2)
        print(paste('Train', cv_fold ,c(MSE,MAE,R2), sep = " "))
        
        market_test = test[,c('StockPrice', 'TimeTillExpiry','CallPrice', 'PutPrice', 'InterestRate')]
        market_test[,c('StockPrice','CallPrice')] = test[,c('StockPrice','CallPrice')]/ test$Strike
        market_test[, 'TimeTillExpiry'] =  test[,'TimeTillExpiry'] / 365
        #market_test[, 'InterestRate'] = test[,"InterestRate"]
        
        t_call = forwardprop_soft(X = as.matrix(cbind(1,market_test[,c('StockPrice','TimeTillExpiry')])), 
                                  Y= as.matrix(market_test[,'CallPrice']), W = m1$W_C, hidden = hidden_layers[[h]],
                                  activation = 'sig', output_activation= 'sig', 
                                  L= length(hidden_layers[[h]]), train_size = nrow(market_test))
        t_put  = forwardprop_soft(X = as.matrix(cbind(1,market_test[,c('StockPrice','TimeTillExpiry')])), 
                                  Y= as.matrix(market_test[,'PutPrice']), W = m1$W_P, hidden = hidden_layers[[h]],
                                  activation = 'sig', output_activation= 'sig', 
                                  L= length(hidden_layers[[h]]), train_size = nrow(market_test))
        
        CALL_CV_MSE = mean((t_call$Yhat-market_test[, 'CallPrice'])^2)
        CALL_CV_MAE = mean(abs(t_call$Yhat-market_test[, 'CallPrice']))
        CALL_CV_R2 = 1 - sum((t_call$Yhat-(market_test[, 'CallPrice']))^2) /
          sum((market_test[, 'CallPrice'] - mean(market_test[, 'CallPrice']))^2)
        
        PUT_CV_MSE = mean((t_put$Yhat-market_test[, 'PutPrice'])^2)
        PUT_CV_MAE = mean(abs(t_put$Yhat-market_test[, 'PutPrice']))
        PUT_CV_R2 = 1 - sum((t_put$Yhat-(market_test[, 'PutPrice']))^2) /
          sum((market_test[, 'PutPrice'] - mean(market_test[, 'PutPrice']))^2)
        
        
        model_res = c(CALL_CV_MSE, CALL_CV_MAE, CALL_CV_R2, PUT_CV_MSE, PUT_CV_MAE, PUT_CV_R2)
        CV_res= rbind(CV_res, model_res)
        
        mlpSC_market_model[[model_no]][[i]] = m1
        
      }
      
      
        
      
      CV_Err = colMeans(CV_res)
      CV_res  = rbind(CV_res,CV_Err)
      CV_res = as.data.frame(CV_res)
      names(CV_res) = c("CALL_CV_MSE", "CALL_CV_MAE", "CALL_CV_R2", "PUT_CV_MSE", "PUT_CV_MAE", "PUT_CV_R2")
      mlpSC_market_cv[[model_no]] = CV_res
      print(paste("#########################################MODEL NUMBER:", model_no, sep = " "))
      model_no = model_no + 1
      print(paste("CV Error:", CV_Err, sep = " "))
      CV_res = c()
    }
  }
}
EndTime = Sys.time()

RunTime = StartTime-EndTime

#setNames(mlpSC_market_model, model_names)
#setNames(mlpSC_market_cv, model_names)


save(mlpSC_market_model, file = paste('mlpSC_market_001_', Sys.Date(), ".RData", sep = ''))
save(mlpSC_market_cv, file = paste('mlpSC_market_cv_001_', Sys.Date(), ".RData", sep = ''))
