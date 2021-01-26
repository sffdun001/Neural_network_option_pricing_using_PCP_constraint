#Duncan Saffy 2019 Dissertation
#Option pricing model

#Version 003 is a neural netork that prices puts and calls and enforces put call parity as well as includes cross validation code
#Version 004 is the same as 003 but now we try use a softplus and RelU function in the output layer to avoid negative prices.
#Version 008 has a BS data simulator that test the soft contraint

rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/HestonTest")
load('Heston_train_final_2020-03-28.RData')
#Defining activation functions and derivatives
sig <- function(x)
{
 return(1/(1+ exp(-x)))
}

dsig <- function(x)
{
  return(sig(x)*(1-sig(x)))
}

dtanh <- function(x)
{
  return(1-(tanh(x))^2)
}

softplus <- function(x)
{
  return(log(1+exp(x)))
}


dsoftplus <- function(x)
{
  return(1/(1+exp(-x)))
}

drelu <- function(x)
{
  return(1*(x>=0))
}

relu <- function(x)
{
  return(x*(x>=0))
}
#Function for fowrard propagation ---------------------------------------------------------------------------------------------------

forwardprop_soft <- function(X,Y, W, hidden, activation = "sig", L, train_size, output_activation)
{ #This function performs the forward propogation
  #X is the input data (each row is one observation)
  #Y is the output data for which to train the model to
  #W is a list of weights - each element in the list is a set of weights for a layer in the NN. NOTE IT INCLUDES WEIGHTS FOR ONES
  #res1#W = matrix(runif(16), 4,4); Y = 1; X = runif(4);activation = "sig"; hidden = 1
  fl_act <- get(output_activation) #final layer activation
  act <- get(activation)
  A = list() # This will be of lenght L+1 and is a list of the activations evaluated at each layer
  A[[1]] = X
  Z = list()
  for (l in 1:L){
    Z[[l]] = A[[l]]%*%W[[l]]
    A[[l+1]] = cbind(1,act(Z[[l]]))
  }
  Z[[L+1]] = A[[L+1]]%*%W[[L+1]]
  Yhat = fl_act(Z[[L+1]]) #linear activation in final layer
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
    d_act = get(paste("d", activation, sep = ""))
    d_fl_act = get(paste("d", output_activation, sep = ""))
    delta_C = W_C; delta_P = W_P; # Set = W some where initialize a list of correct length
    delta_W_C = W_C; delta_W_P = W_P
    delta_C[[L+1]] = (res_call$Yhat - res_call$Y)*d_fl_act(res_call$Z[[L+1]]) + 
                     (res_call$Yhat + exp(-rf * X[,(t+1)])- res_put$Yhat - exp(-rf * X[,(t+1)])*X[,(S+1)])*d_fl_act(res_call$Z[[L+1]])# derivative wr.r.t weights connection L -> output layer
    delta_W_C[[L+1]] = t(res_call$A[[L+1]])%*%delta_C[[L+1]]
    delta_P[[L+1]] = (res_put$Yhat - res_put$Y)*d_fl_act(res_put$Z[[L+1]]) - 
                     (res_call$Yhat + exp(- rf * X[,(t+1)])- res_put$Yhat - exp(-rf * X[,(t+1)])*X[,(S+1)])*d_fl_act(res_put$Z[[L+1]])
    delta_W_P[[L+1]] = t(res_put$A[[L+1]])%*%delta_P[[L+1]]
    for (l in 1:L){
      delta_C[[L+1-l]] = delta_C[[L+2-l]]%*%t(W_C[[L+2-l]])*d_act(res_call$Z[[L+1-l]])
      delta_W_C[[L+1-l]] = t(res_call$A[[L+1-l]])%*%delta_C[[L+1-l]]
      delta_P[[L+1-l]] = delta_P[[L+2-l]]%*%t(W_P[[L+2-l]])*d_act(res_put$Z[[L+1-l]])
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
    forward_call = forwardprop_soft(X = X_star, Y= C_star, W = W_C, hidden = hidden,
                               activation = activation, output_activation= output_activation, L= L, train_size = train_size)
    forward_put = forwardprop_soft(X = X_star, Y= P_star, W = W_P, hidden = hidden, 
                              activation = activation, output_activation= output_activation, L= L, train_size = train_size)
    back = backprop_soft(X = X_star,  rf = rf, t = t, S = S, C= C_star, P=P_star, W_C = W_C, W_P = W_P, 
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
    PE = sum(((FC$Yhat + exp(- rf * X_bar[,(t+1)])- FP$Yhat - exp(- rf * X_bar[,(t+1)])*X_bar[,(S+1)]))^2) / nrow(X_bar)
    
    call_error_tracker[m] = FC$Err
    put_error_tracker[m] = FP$Err
    parity_error_tracker[m] = PE
    
    learning_rate = learning_rate #* 0.9
    
    print(c(FC$Err, FP$Err, PE))
  }
  
  
  return(list(W_C= W_C, W_P = W_P, Call_Error = call_error_tracker, Put_Error = put_error_tracker, Parity_Error = parity_error_tracker, 
              Chat = FC$Yhat, Phat = FP$Yhat, C = C, P = P))
}


set.seed(1995)

mlpSC_hest_model <- vector(mode = "list", length = 10)

for (i in 1:10){
  Data = df_hest_train[[i]]
  hest_data = Data[,c('StockPrice', 'TimeTillExpiry','CallPrice', 'PutPrice')]
  hest_data[,c('StockPrice','CallPrice', 'PutPrice')] = Data[,c('StockPrice','CallPrice', 'PutPrice')]/ Data$Strike
  hest_data[, 'TimeTillExpiry'] =  Data[,'TimeTillExpiry'] / 365
  
  m1 <- neural_net_soft(X = as.matrix(hest_data[,c('StockPrice','TimeTillExpiry')]), 
                       C = as.matrix(hest_data[, 'CallPrice']), P= as.matrix(hest_data[, 'PutPrice']),
                       S = 1, t = 2, hidden = c(4), iter_lim  = 50, learning_rate = 0.1, rf = 0.0687153,
                       train_size = 100, lambda = 0, activation = 'sig', output_activation = 'sig')
  
  MSE = mean((m1$Chat-hest_data[, 'CallPrice'])^2)
  MAE = mean(abs(m1$Chat-hest_data[, 'CallPrice']))
  R2 = 1 - sum((m1$Chat-(hest_data[, 'CallPrice']))^2) / 
    sum((hest_data[, 'CallPrice'] - mean(hest_data[, 'CallPrice']))^2)
  print(c(MSE,MAE,R2))
  mlpSC_hest_model[[i]] = m1
}


save(mlpSC_hest_model, file = paste('mlpSC_hest_001_', Sys.Date(), ".RData", sep = ''))
