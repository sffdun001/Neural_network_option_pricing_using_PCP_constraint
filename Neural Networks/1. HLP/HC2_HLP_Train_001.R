#Duncan Saffy 2019 Dissertation
#Option pricing model


rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/Technical Paper/Hard Constraint")
load('Hutch_train_001_2019-10-31.RData')
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

#function for backpropagation-----------------------------------------------------------------------------------------------------

backprop_HC2 <- function(X, rf, t, S, C, P, W_C, W_P, hidden, learning_rate, activation, 
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
    d_act = get(paste("d", activation, sep = ""))
    d_fl_act = get(paste("d", output_activation, sep = ""))
    delta_C = W_C; delta_P = W_P; # Set = W some where initialize a list of correct length
    delta_W_C = W_C; delta_W_P = W_P
    
    delta_C[[L+1]] = ((res$chat - res$C) + (res$phat - res$P))*d_fl_act(res$Z_C[[L+1]])  #(res$chat - res$C) + (res$phat - res$P) +
    delta_W_C[[L+1]] = t(res$A_C[[L+1]])%*%delta_C[[L+1]]
    
    delta_P[[L+1]] = ((res$chat - res$C) + (res$phat - res$P))*d_fl_act(res$Z_P[[L+1]]) 
    delta_W_P[[L+1]] = t(res$A_P[[L+1]])%*%delta_P[[L+1]]
    
    
    for (l in 1:L){
      delta_C[[L+1-l]] = delta_C[[L+2-l]]%*%t(W_C[[L+2-l]])*d_act(res$Z_C[[L+1-l]])
      delta_W_C[[L+1-l]] = t(res$A_C[[L+1-l]])%*%delta_C[[L+1-l]]
      delta_P[[L+1-l]] = delta_P[[L+2-l]]%*%t(W_P[[L+2-l]])*d_act(res$Z_P[[L+1-l]])
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

  for (m in 1:iter_lim) #or some alternative condition
  { 
    random_vec = sample(c(1:nrow(X_bar)), size = train_size*100000, replace = TRUE)
    for (n in 1:100000){
    samp = random_vec[((n-1)*train_size+1):(n*train_size)]
    X_star = X_bar[samp,]
    C_star = C[samp,]
    P_star = P[samp,]
    
    forward = forwardprop_HC2(X = X_star, t = t, rf = rf, C= C_star, P = P_star, S = S, W_C = W_C, W_P = W_P, hidden = hidden,
                               activation = activation, output_activation= output_activation, L= L, train_size = train_size)
    
    back = backprop_HC2(X = X_star, t = t, S = S, C= C_star, P=P_star, W_C = W_C, W_P = W_P, 
                    hidden = hidden, learning_rate = learning_rate, activation = activation, output_activation= output_activation,
                    L = L, res= forward, lambda = lambda, train_size = train_size)
    W_C = back$W_C
    W_P = back$W_P
   
    }
    
    learning_rate = learning_rate
    Fwd = forwardprop_HC2(X = X_bar, C= C, t = t, r = rf, P = P, S = S, W_C = W_C, W_P =W_P, hidden = hidden,
                    activation = activation, output_activation= output_activation, L= L, train_size = nrow(X_bar))
    
    #PE = sum(((FC$Yhat + exp(-(X_bar[,(rf+1)]) * X_bar[,(t+1)])- FP$Yhat - X_bar[,(S+1)]))^2) / nrow(X_bar)
    
    parity_error_tracker[m] = Fwd$CallError
    call_error_tracker[m] = Fwd$cStarErr
    put_error_tracker[m] = Fwd$pStarErr
    #parity_error_tracker[m] = PE
    print(c(Fwd$CallError,  Fwd$cStarErr, Fwd$pStarErr))
  }
  
  
  return(list(W_C= W_C, W_P = W_P, Call_Error = call_error_tracker, Put_Error = put_error_tracker, 
              parity_error_tracker = parity_error_tracker, Chat = Fwd$chat, Phat = Fwd$phat, C = C, P = P))
}


set.seed(1995)


mlpHC2_sig_model2 <- vector(mode = "list", length = 10)

for (i in 1:10){
  #i =2
  Data = TRAIN_DATA[[i]]
  hutch_data = Data[,c('StockPrice', 'TimeTillExpiry','CallPrice', 'PutPrice')]
  hutch_data[,c('StockPrice','CallPrice', 'PutPrice')] = Data[,c('StockPrice','CallPrice', 'PutPrice')]/ Data$Strike
  hutch_data[, 'TimeTillExpiry'] =  Data[,'TimeTillExpiry'] / 252
  
  m1 <- neural_net_HC2(X = as.matrix(hutch_data[,c('StockPrice','TimeTillExpiry')]), 
                      C = as.matrix(hutch_data[, 'CallPrice']), P= as.matrix(hutch_data[, 'PutPrice']),
                      S = 1, t = 2, hidden = c(4), iter_lim  = 30, learning_rate = 0.1, rf = 0.04681,
                      train_size = 100, lambda = 0, activation = 'sig', output_activation = 'sig')
  
  MSE = mean((m1$Chat-hutch_data[, 'CallPrice'])^2)
  MAE = mean(abs(m1$Chat-hutch_data[, 'CallPrice']))
  R2 = 1 - sum((m1$Chat-(hutch_data[, 'CallPrice']))^2) / 
    sum((hutch_data[, 'CallPrice'] - mean(hutch_data[, 'CallPrice']))^2)
  print(c(MSE,MAE,R2))
  mlpHC2_sig_model2[[i]] = m1
}

mlpHC4_sig_model2 = mlpHC2_sig_model2
save(mlpHC4_sig_model2, file = paste('mlpHC4_model_002_', Sys.Date(), ".RData", sep = ''))
mlpHC2_sig_model[[]]$Call_Error
mlpHC2_sig_model[[1]]$Put_Error


#load('mlpHC2_model_001_2019-11-04.RData')
# Verifying PC parity
# plotting portfolio values


price_seq = seq(0.8, 1.2, length = 100)
time_seq = rep(1/252, 100)

x.1 = cbind(1,price_seq,time_seq)
test.1 = forwardprop_HC2(X = x.1, S = 1, rf = 0.04861, t= 2, C = time_seq, P = time_seq, 
                         W_C = m1$W_C, W_P = m1$W_P, hidden = c(4), activation = "sig", 
                         L =1, train_size =100, output_activation = 'sig')

x = price_seq
y1 = test.1$chat + 1
y2 = test.1$phat +price_seq
plot(x ,y1, type = 'l')
lines(x,y2, type = 'l', col = 'red')


min(test.1$phat)


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


delta_HC = Derivative_HC2(X = x.1, S = 1, C = time_seq, P = time_seq, res = test.1, 
                          W_C = m1$W_C, W_P = m1$W_P, hidden = c(4), activation = "sig", 
                          L =1, output_activation = 'sig')
delta_HC$call_derivative







