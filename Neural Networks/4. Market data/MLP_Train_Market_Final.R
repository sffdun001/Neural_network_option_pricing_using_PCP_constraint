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

#This function performs the forward propogation
#X is the input data (each row is one observation), the first column is a column of 1's
#Y is the output data for which to train the model to
#W is a list of weights - each element in the list is a set of weights for a layer in the NN. NOTE IT INCLUDES WEIGHTS FOR ONES
#L is the number of hidden layers
#activation  is the activation in the hidden layers
#output activation is the activation in the final layer
#sample_size is the batch size of training

forwardprop <- function(X,Y, W, hidden, activation = "sig", L,
                        train_size, output_activation)
{ 
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
  Yhat = sig(Z[[L+1]])
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

  #d_act = get(paste("d", activation, sep = ""))
  #d_fl_act = get(paste("d", output_activation, sep = ""))
  #delta = W # Set = W some where initialize a list of correct length
  delta_W = W
  #delta[[L+1]] = (res$Yhat - res$Y)*dsig(res$Z[[L+1]])
  delta = (res$Yhat - res$Y)*dsig(res$Z[[L+1]])
  delta_W[[L+1]] = t(res$A[[L+1]])%*%delta
  for (l in 1:L){
    #delta[[L+1-l]] = delta[[L+2-l]]%*%t(W[[L+2-l]])*dsig(res$Z[[L+1-l]])
    delta = delta%*%t(W[[L+2-l]])*dsig(res$Z[[L+1-l]])
    delta_W[[L+1-l]] = t(res$A[[L+1-l]])%*%delta
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
  cntr = 0
  PrevErr = 0 
  for (m in 1:iter) #or some alternative condition
  {
    random_vec = sample(c(1:nrow(X_bar)), size = train_size*1e5, replace = TRUE)
    for (n in 1:1e5){
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
    
    
    PrevErr = FY$Err
  }


  return(list(W= W, Error = error_tracker, Yhat = FY$Yhat, Y = Y))
}


#Training neural network

set.seed(2020)

train = data_main[[1]]
for(i in 2:6){
  train = rbind(train, data_main[[i]])
}



market_data = train[,c('StockPrice', 'TimeTillExpiry','CallPrice')]
market_data[,c('StockPrice','CallPrice')] = train[,c('StockPrice','CallPrice')]/ train$Strike
market_data[, 'TimeTillExpiry'] =  train[,'TimeTillExpiry'] / 365


m1 <- neural_net_single(X = as.matrix(market_data[,c('StockPrice','TimeTillExpiry')]),
                        Y= as.matrix(market_data[, 'CallPrice']),
                        hidden = c(4,4), iter = 50, learning_rate = 0.1,
                        train_size = 100, lambda = 0.00001, activation = 'sig', output_activation = 'sig')

mlp_final = m1
save(mlp_final, file = paste('mlp_market_final_', Sys.Date(), ".RData", sep = ''))



MSE = mean((m1$Yhat-market_data[, 'CallPrice'])^2)
MAE = mean(abs(m1$Yhat-market_data[, 'CallPrice']))
R2 = 1 - sum((m1$Yhat-(market_data[, 'CallPrice']))^2) /
  sum((market_data[, 'CallPrice'] - mean(market_data[, 'CallPrice']))^2)
print(paste('Train', cv_fold ,c(MSE,MAE,R2), sep = " "))

market_test = test[,c('StockPrice', 'TimeTillExpiry','CallPrice')]
market_test[,c('StockPrice','CallPrice')] = test[,c('StockPrice','CallPrice')]/ test$Strike
market_test[, 'TimeTillExpiry'] =  test[,'TimeTillExpiry'] / 365


t1 = forwardprop(X = as.matrix(cbind(1,market_test[,c('StockPrice','TimeTillExpiry')])),
                 Y= as.matrix(market_test[,'CallPrice']),
                 W = m1$W , hidden = hidden_layers[[h]], activation = 'sig',
                 output_activation= 'sig', L= length(hidden_layers[[h]]), train_size = nrow(market_test))










