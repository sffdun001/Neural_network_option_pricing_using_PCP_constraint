#Duncan Saffy 2019 Dissertation
#Option pricing model

#Version 003 is a neural netork that prices puts and calls and enforces put call parity as well as includes cross validation code
#Version 004 is the same as 003 but now we try use a softplus and RelU function in the output layer to avoid negative prices.
#version 005 only prices puts or calls. Not both.
#Version 007 only prices puts or calls and 



rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/Technical Paper/Only Puts and Calls")
load('Hutch_train_001_2019-10-29.RData')
#Load in training data

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

identity <- function(x)
{
  return(x)
}

didentity <- function(x)
{
  return(1)
}

#Function for fowrard propagation ---------------------------------------------------------------------------------------------------

forwardprop <- function(X,Y, W, hidden, activation = "sig", L, 
                        train_size, output_activation)
{ #This function performs the forward propogation
  #X is the input data (each row is one observation), the first column is a column of 1's
  #Y is the output data for which to train the model to
  #W is a list of weights - each element in the list is a set of weights for a layer in the NN. NOTE IT INCLUDES WEIGHTS FOR ONES
  #L is the number of hidden layers
  #activation  is the activation in the hidden layers
  #output activation is the activation in the final layer
  #sample_size is the batch size of training
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
  Yhat = fl_act(Z[[L+1]]) 
  Err = sum((Yhat-Y)^2) / train_size
  
  return(list(Yhat = Yhat, A = A, Err = Err, X = X, Y = Y, Z = Z))
}

#function for backpropagation-----------------------------------------------------------------------------------------------------

backprop<- function(X, Y, W, hidden, learning_rate, activation, 
                    res , L, lambda, output_activation)
{#res is the output of the forward propagation
  #for this function only W_bar represent weight matricies that include bias weights and W is without bias weights
  W_bar = W
  for (l in 1:(L+1)){
    W[[l]] = W[[l]][-1,,drop = F] #drop = F prevents a single column matrix from becoming type = numeric
  }
  
  d_act = get(paste("d", activation, sep = ""))
  d_fl_act = get(paste("d", output_activation, sep = ""))
  delta = W # Set = W some where initialize a list of correct length
  delta_W = W
  delta[[L+1]] = (res$Yhat - res$Y)*d_fl_act(res$Z[[L+1]])
  delta_W[[L+1]] = t(res$A[[L+1]])%*%delta[[L+1]]
  for (l in 1:L){
    delta[[L+1-l]] = delta[[L+2-l]]%*%t(W[[L+2-l]])*d_act(res$Z[[L+1-l]])
    delta_W[[L+1-l]] = t(res$A[[L+1-l]])%*%delta[[L+1-l]]
  }
  #updating Weights
  for (l in 1:(L+1)){
    W_bar[[l]] = W_bar[[l]] - learning_rate * delta_W[[l]] - learning_rate*lambda * rbind(0,W[[l]])
  }
  
  #Err   = sum((res$Yhat-Y)^2) / batch_size
  
  return(list(W= W_bar))
}


#Neuralnetwork function -------------------------------------------------------------------------------------------------------------

neural_net_single <- function(X, Y, hidden, iter, learning_rate, 
                              activation = 'sig', train_size, lambda = 0, 
                              output_activation)
{ # This function uses stoc_NN to update weights
  # iter is the number of iterations the of backpropagation we will do
  # learning_rate 
  # sample_size is the number of X values selected for trainging at each epoch
  #rf, t and S are the indexes for the risk free, time and stock price in the input matrix X.
  X_bar = cbind(1,X) #adding a column of ones to design matrix
  p = train_size
  L = length(hidden)
  q = ncol(X_bar)-1
  W = list()
  for (i in 1:(L+1)){ #initialising weights 
    if (i == L+1){
      k = 1
    }
    else{
      k = hidden[i]
    }
    W[[i]] = matrix(runif((q+1)*k,-1,1), nrow = q+1, ncol = k)
    q = k
  }
  
  error_tracker = numeric(iter) #This is just something so we can track the error after each epoch
  
  for (m in 1:iter) #or some alternative condition
  { 
    random_vec = sample(c(1:nrow(X_bar)), size = train_size*100000, replace = TRUE)
    for (n in 1:100000){
      samp = random_vec[((n-1)*train_size+1):(n*train_size)]
      X_star = X_bar[samp,]
      Y_star = Y[samp,]
      forward = forwardprop(X = X_star, Y= Y_star, W = W, hidden = hidden,
                            activation = activation, output_activation= output_activation, L= L, train_size = train_size)
      back = backprop(X = X_star, Y= Y_star, W = W, hidden = hidden, learning_rate = learning_rate, 
                      activation = activation, output_activation= output_activation, L = L, 
                      res= forward, lambda = lambda)
      
      W = back$W
    }
    FY = forwardprop(X = X_bar, Y= Y, W = W, hidden = hidden, activation = activation, 
                     output_activation= output_activation, L= L, train_size = nrow(X_bar))
    error_tracker[m] = FY$Err
    print(FY$Err)
    learning_rate = learning_rate #*0.95
  }
  
  
  return(list(W= W, Error = error_tracker, Yhat = FY$Yhat, Y = Y))
}


#Training neural network

set.seed(1995)

mlp_sig_model <- vector(mode = "list", length = 10)

for (i in 1:10){
  #i = 1
  Data = TRAIN_DATA[[i]]
  hutch_data = Data[,c('StockPrice', 'TimeTillExpiry','CallPrice')]
  hutch_data[,c('StockPrice','CallPrice')] = Data[,c('StockPrice','CallPrice')]/ Data$Strike
  hutch_data[, 'TimeTillExpiry'] =  Data[,'TimeTillExpiry'] / 252
  
  m1 <- neural_net_single(X = as.matrix(hutch_data[,c('StockPrice','TimeTillExpiry')]), 
             Y= as.matrix(hutch_data[, 'CallPrice']), 
             hidden = c(4), iter = 50, learning_rate = 0.1, 
             train_size = 100, lambda = 0, activation = 'sig', output_activation = 'sig')
  
  MSE = mean((m1$Yhat-hutch_data[, 'CallPrice'])^2)
  MAE = mean(abs(m1$Yhat-hutch_data[, 'CallPrice']))
  R2 = 1 - sum((m1$Yhat-(hutch_data[, 'CallPrice']))^2) / 
    sum((hutch_data[, 'CallPrice'] - mean(hutch_data[, 'CallPrice']))^2)
  print(c(MSE,MAE,R2))
  mlp_sig_model[[i]] = m1
}

save(mlp_sig_model, file = paste('mlp_model_001_', Sys.Date(), ".RData", sep = ''))
  


CallDelta = function(X, Weights, res, hidden, S, activation, 
                     output_activation)
{ #This function finds the derivative of option with respect to the underlying. 
  X = as.matrix(hutch_data[,c('StockPrice','TimeTillExpiry')])
  Weights = m1$W ; res = F.Delta; hidden = c(4); S = 1
  activation = 'sig'; output_activation= 'sig'
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

F.Delta = forwardprop(X = as.matrix(cbind(1,hutch_data[,c('StockPrice','TimeTillExpiry')])), 
                      Y= as.matrix(hutch_data[,'CallPrice']), 
                      W = m1$W, hidden = c(4), activation = 'sig', 
                      output_activation= 'sig', L= 1, train_size = nrow(hutch_data))
t.1 = CallDelta(X = as.matrix(hutch_data[,c('StockPrice','TimeTillExpiry')]),
                Weights = m1$W, res = F.Delta, hidden = c(4), S = 1,
                activation = 'sig', 
                output_activation= 'sig')



