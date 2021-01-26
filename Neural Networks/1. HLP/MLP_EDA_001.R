#Duncan Saffy 2019 Dissertation
#Option pricing model

#Version 003 is a neural netork that prices puts and calls and enforces put call parity as well as includes cross validation code
#Version 004 is the same as 003 but now we try use a softplus and RelU function in the output layer to avoid negative prices.
#version 005 only prices puts or calls. Not both.
#Version 007 only prices puts or calls and 



rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/Technical Paper/Only Puts and Calls")
load('mlp_model_001_2019-11-05.RData')
load('Hutch_test_001_2019-10-31.RData')
#Load in training data

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

BlackScholes <- function(S, K, t, r, sd, call = TRUE, DT = 1/252) #So this apparently works
{ 
  d1 <- 1/(sd*sqrt(t*DT)) * (log(S/K) +((r + (sd^2)/2)*(t*DT)))
  d2 <- d1 - sd*sqrt(t*DT)
  C <- S*pnorm(d1) - K*exp(-r*t*DT)*pnorm(d2)
  if(call){
    return(C) #returning the call price
  }
  else {
    return(C + K*exp(-r*t*DT) - S) #returning the put price
  }
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


##############################################################################
#Plotting Surfaces
##############################################################################
set.seed(1995)

load('Hutch_train_001_2019-10-31.RData')


data <-  TRAIN_DATA[[1]]%>%mutate(std.StockPrice = StockPrice / Strike)%>%
  arrange(desc(std.StockPrice))%>%arrange(desc(TimeTillExpiry))




n <- 100
time_seq <- 1:max(data$TimeTillExpiry) / 252
price_seq <- seq(from = 0.79, to = 1.5 ,length = n)

predicted_price= matrix(NA, nrow= n, ncol = max(data$TimeTillExpiry))
predicted_derivative= matrix(NA, nrow= n, ncol = max(data$TimeTillExpiry))
bs = matrix(NA, nrow= n, ncol = max(data$TimeTillExpiry))
bs_derivative = matrix(NA, nrow= n, ncol = max(data$TimeTillExpiry))

X = matrix(NA, nrow = n*max(data$TimeTillExpiry), ncol = 2)

time_new = c()
price_new = c()
predicted_price_new = c()
predicted_derivative_new = c()

min_price = 0.95
max_price = 1.05

model <- mlp_sig_model[[3]]

for (i in max(data$TimeTillExpiry):1){
  #i = 168
  #Setting input
  min_price = max(0.8, min_price - 0.0025)
  #max_price = max(max_price, data$std.StockPrice[which(data$TimeTillExpiry == i)])
  max_price = min(max_price +0.004, 140)
  
  X[((i-1)*n+1):(n*i),] = cbind(price_seq, i/252)
  
  #Calcating model value and derivative
  output = forwardprop(X = as.matrix(cbind(1,X[((i-1)*n+1):(n*i),]))[,,drop = F], 
                       Y = as.matrix(X[((i-1)*n+1):(n*i),1]) , #doesnt matter what this vector is..
                       W = model$W, hidden = c(4), train_size = 100, output_activation = 'sig',
                       activation = 'sig', L = 1)
  
  #forwardprop <- function(X,Y, W, hidden, activation = "sig", L, 
  #                        train_size, output_activation)
  
  output_derivative = CallDelta(X = as.matrix(cbind(1,X[((i-1)*n+1):(n*i),]))[,,drop = F], S = 1,
                                      Weights = model$W, hidden = c(4), output_activation = 'sig',
                                      activation = 'sig', res = output)
  # CallDelta = function(X, Weights, res, hidden, S, activation, 
  #                     output_activation)
  
  #Calculating BS derivative
  d1 = 1/(0.2*sqrt(X[((i-1)*n+1):(n*i),2]))*(log(X[((i-1)*n+1):(n*i),1]) +
                                               ((0.04681 + (0.2^2)/2)*(X[((i-1)*n+1):(n*i),2])))
  d2 <- d1 - 0.2*sqrt(X[((i-1)*n+1):(n*i),2])
  bs[((i-1)*n+1):(n*i)] =  X[((i-1)*n+1):(n*i),1]*pnorm(d1) - 
    exp(-0.04681*X[((i-1)*n+1):(n*i),2])*pnorm(d2)
  bs_derivative[((i-1)*n+1):(n*i)] =  pnorm(d1)
  
  
  #Storing values
  ind = which((X[((i-1)*n+1):(n*i),1] > min_price) & (X[((i-1)*n+1):(n*i),1] < max_price))

  predicted_price[ind,i] = output$Yhat[ind]
  predicted_derivative[ind,i] = output_derivative$Delta[ind]
  time_new = c(time_new, time_seq)
  price_new = c(price_new, rep(price_seq[i], n))
  predicted_price_new = c(predicted_price_new,
                          predicted_price[ind,i])
  predicted_derivative_new = c(predicted_derivative_new, 
                               predicted_derivative[ind,i])
}


row_remove = seq(from = 2, to = nrow(predicted_price) , by = 2)
col_remove = seq(from = 2, to = ncol(predicted_price) , by = 2)  

price_seq = price_seq[-row_remove]
time_seq = time_seq[-col_remove]
predicted_price =predicted_price[-row_remove, -col_remove]
bs = bs[-row_remove, -col_remove]
predicted_derivative =predicted_derivative[-row_remove, -col_remove]
bs_derivative = bs_derivative[-row_remove, -col_remove]

col_remove2 = seq(from = 2, to = ncol(predicted_price) , by = 2) 

time_seq = time_seq[-col_remove2]
predicted_price =predicted_price[, -col_remove2]
bs = bs[, -col_remove2]
predicted_derivative =predicted_derivative[, -col_remove2]
bs_derivative = bs_derivative[, -col_remove2]


pdf('mlp_call_001.pdf')  
persp(x = price_seq, y =time_seq, z = predicted_price, phi = 25, theta = 320,
      xlab = "S/K", ylab = "T-t", zlab = '', ticktype = 'detailed', d = 4, expand = 0.5,
      zlim = c(0, 0.5)
)
mtext(expression(frac(hat(c), K)),side=2,las=1,line=1, at = c(0,0))
dev.off()


pdf('mlp_delta_001.pdf')  
persp(x = price_seq, y = time_seq, z = predicted_derivative, phi = 25, theta = 320,
      xlab = "Moneyness", ylab = "T-t", zlab = '', ticktype = 'detailed', expand = 0.5, las = 0,
      d =4, zlim = c(-0.3,1.3)
)
mtext(expression(Delta),side=2,las=2,line=0.5, at = c(0,0), cex = 1)
dev.off()


pdf('mlp_calldiff_001.pdf')  
persp(x = price_seq, y =time_seq, z = predicted_price  - bs, phi = 25, theta = 320,
      main = NULL, ticktype = 'detailed', d = 4, expand = 0.5,
      xlab = "Moneyness", ylab = "T-t", zlab = '', zlim = c(-0.05, 0.02)
)
mtext(expression(paste(frac(hat(c), K), "  Error", sep = " ")),side=2,las=0,line=0, at = c(0,0))
dev.off()



pdf('mlp_deltadiff_001.pdf')  
persp(x = price_seq, y = time_seq, z = predicted_derivative - bs_derivative, phi = 25, theta = 320,
      xlab = "Moneyness", ylab = "T-t", zlab = '', d =4, ticktype = 'detailed', expand = 0.5, shade = 0.01
)
mtext(expression(paste(Delta, "  Error", sep = " ")),side=2,las = 0,line= 0.5, at = c(2,0))
dev.off()







