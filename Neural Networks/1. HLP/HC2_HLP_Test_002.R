#Duncan Saffy 2019 Dissertation
#Option pricing model


rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/MLP HC Hutch")
load('Hutch_test_001_2019-10-31.RData')
load('mlpHC4_model_001_2019-11-18.RData') 

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
    
    delta_C[[L+1]] = ((res$chat_star - res$C))*d_fl_act(res$Z_C[[L+1]])  #(res$chat - res$C) + (res$phat - res$P) +
    delta_W_C[[L+1]] = t(res$A_C[[L+1]])%*%delta_C[[L+1]]
    
    delta_P[[L+1]] = ((res$phat_star - res$P))*d_fl_act(res$Z_P[[L+1]]) 
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


#delta_HC = Derivative_HC2(X = x.1, S = 1, C = time_seq, P = time_seq, res = test.1, 
#                          W_C = m1$W_C, W_P = m1$W_P, hidden = c(4), activation = "sig", 
#                          L =1, output_activation = 'sig')
#delta_HC$call_derivative
library(dplyr)

set.seed(1995)

MLP_HC4_PriceError = vector(mode = 'list', length = 10)

for (m in 1:10){
  #m = 1
  model <-  mlpHC4_sig_model[[m]]
  pe = matrix(NA, nrow = 500, ncol = 4)
  for(i in 1:500){
    #i = 1
    data <-  TEST_DATA[[i]]%>%filter(TerminalDate == 63 & Strike == 50)%>%
      arrange(Day)
    data$StockPrice = data$StockPrice /data$Strike
    data$TimeTillExpiry = data$TimeTillExpiry / 252
    
    m1 = forwardprop_HC2(X = as.matrix(cbind(1,data[,c('StockPrice','TimeTillExpiry')])), 
                     C = as.matrix(data[,'CallPrice']), P  = as.matrix(data[,'PutPrice']),
                     W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                     output_activation= 'sig', L= 1, train_size = nrow(data),
                     S =1, t =2, rf = 0.04861)
    
    MSE = mean((m1$chat*(m1$chat>0)-data$CallPrice /data$Strike)^2)
    MAE = mean(abs(m1$chat*(m1$chat>0)-data$CallPrice /data$Strike))
    R2 = 1 - sum((m1$chat*(m1$chat>0)-data$CallPrice  /data$Strike)^2) / 
      sum((data$CallPrice /data$Strike - mean(data$CallPrice /data$Strike))^2)
    R2_call_star = 1 - sum((m1$chat_star*(m1$chat>0) - data$CallPrice  /data$Strike)^2) / 
      sum((data$CallPrice /data$Strike - mean(data$CallPrice /data$Strike))^2)
    pe[i,] = c(MSE,MAE,R2, R2_call_star)
  }
  MLP_HC4_PriceError[[m]] = pe
}


r2mean = c()
r2min = c()
r2max = c()
r2mean_call_star = c()
r2min_call_star = c()
r2max_call_star = c()
for (i in 1:10){
  r2mean = c(r2mean, mean(MLP_HC4_PriceError[[i]][,3]))
  r2min = c(r2min, min(MLP_HC4_PriceError[[i]][,3]))
  r2max = c(r2max, max(MLP_HC4_PriceError[[i]][,3]))
  r2mean_call_star = c(r2mean_call_star, mean(MLP_HC4_PriceError[[i]][,4]))
  r2min_call_star = c(r2min_call_star, min(MLP_HC4_PriceError[[i]][,4]))
  r2max_call_star = c(r2max_call_star, max(MLP_HC4_PriceError[[i]][,4]))
}

max(r2max)
mean(r2mean)
min(r2min)
r2max_call_star
r2mean_call_star
r2min_call_star

price_error = as.data.frame(matrix(NA, nrow = 4, ncol = 1))
colnames(price_error) = "HC2"
price_error[1,1] = paste(round(mean(r2mean)*100,2), '%',sep = "")
price_error[2,1] = paste(round(sd(r2mean)*100,2), '%',sep = "")
price_error[3,1] = paste(round(max(r2max)*100,2), '%',sep = "")
price_error[4,1] = paste(round(min(r2min)*100,2), '%',sep = "")
#min(r2min)
write.csv(price_error, file = paste('mlpHC4_price_error_table_001_', Sys.Date(), '.csv', sep = '')) 


save(MLP_HC4_PriceError, file = paste('mlp_HC4_price_error_001_', Sys.Date(), '.csv', sep = ''))



##############################################################################
#Hedging Error Calculations
##############################################################################

MLP_HC4_HedgeError = vector(mode = 'list', length = 5)
MLP_HC4_HedgeError_call_star = vector(mode = 'list', length = 5)

BSHedgeError = vector(mode = 'list', length = 5)


#Calculating 3 month hedging error.
TimeToExpiry = 63
for (K in c(40, 45, 50, 55, 60)){
  #K = 60
  for (m in 1:10){
   # K=50
   #m=1
    model <-  mlpHC4_sig_model[[m]]
    
    S0 =50
    N = 500
    V_S = matrix(NA, nrow = 253, ncol = N)
    V_C = matrix(NA, nrow = 253, ncol = N)
    V_B = matrix(NA, nrow = 253, ncol = N)
    V = matrix(NA, nrow = 253, ncol = N)
    
    V_call_star_S = matrix(NA, nrow = 253, ncol = N)
    V_call_star_B = matrix(NA, nrow = 253, ncol = N)
    V_call_star = matrix(NA, nrow = 253, ncol = N)
    
    V_bs_S = matrix(NA, nrow = 253, ncol = N)
    V_bs_B = matrix(NA, nrow = 253, ncol = N)
    V_bs = matrix(NA, nrow = 253, ncol = N)
    
    m1_call_star = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/252))[,,drop = F]), 
                     Y= as.matrix(1), 
                     W = model$W_C , hidden = c(4), activation = 'sig', 
                     output_activation= 'sig', L= 1, train_size = 1)
    
    
    MLP_HC_Delta_call_star = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/252))[,,drop = F]), 
                         Weights = model$W_C, hidden = c(4), res = m1_call_star, activation = 'sig', 
                         output_activation= 'sig', S =1)
    

    m1 = forwardprop_HC2(X = (as.matrix(cbind(1, (S0/K), TimeToExpiry/252))[,,drop = F]), 
                         C = as.matrix(1), P  = as.matrix(1),
                         W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                         output_activation= 'sig', L= 1, train_size = 1,
                         S =1, t =2, rf = 0.04861)
    
  
    MLP_HC_Delta = Derivative_HC2(X = (as.matrix(cbind( (S0/K), TimeToExpiry/252))[,,drop = F]), 
                                  C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                  W_C = model$W_C, hidden = c(4), res = m1, activation = 'sig', 
                                  output_activation= 'sig', S =1, L = 1)
    

    
    V_S[1,] = S0 * MLP_HC_Delta$call_derivative *(m1$chat>0)
    V_C[1,] = -BlackScholes(S = rep(S0, N), K = K, t = TimeToExpiry, r = 0.04681, sd = 0.2) 
    V_B[1,] = -(V_S[1,] + V_C[1,]) 
    V[1,]   = V_S[1,] + V_C[1,] + V_B[1,]
    
    V_call_star_S[1,] = S0 * MLP_HC_Delta_call_star$Delta*(m1$chat>0)
    V_call_star_B[1,] = -(V_call_star_S[1,] + V_C[1,]) 
    V_call_star[1,]   = V_call_star_S[1,] + V_C[1,] + V_call_star_B[1,]
    
    d1 = 1/(0.2*sqrt(TimeToExpiry/252))*(log(S0/K) +
                                 ((0.04681 + (0.2^2)/2)*(TimeToExpiry/252)))
    BSDelta = pnorm(d1)
    V_bs_S[1,] = S0*BSDelta*(m1$chat>0)
    V_bs_B[1,] = -(V_bs_S[1,] + V_C[1,]) 
    V_bs[1,]   = V_bs_S[1,] + V_C[1,] + V_bs_B[1,]
    
    for (i in 1:500){
      #i =1
      data <-  TEST_DATA[[i]]%>%filter(TerminalDate == TimeToExpiry & Strike == K)%>%
        arrange(Day)
      data$StockPrice = data$StockPrice / K 
      data$TimeTillExpiry = data$TimeTillExpiry / 252
      
      d1 = 1/(0.2*sqrt(TimeToExpiry/252))*(log(S0/K) +
                                   ((0.04681 + (0.2^2)/2)*(TimeToExpiry/252)))

      m1_call_star = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/252))[,,drop = F]), 
                                 Y= as.matrix(1), 
                                 W = model$W_C , hidden = c(4), activation = 'sig', 
                                 output_activation= 'sig', L= 1, train_size = 1)
      
      
      MLP_HC_Delta_call_star = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/252))[,,drop = F]), 
                                         Weights = model$W_C, hidden = c(4), res = m1_call_star, activation = 'sig', 
                                         output_activation= 'sig', S =1)
      
      
      m1 = forwardprop_HC2(X = (as.matrix(cbind(1, S0/K, TimeToExpiry/252))[,,drop = F]), 
                           C = as.matrix(1), P  = as.matrix(1),
                           W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                           output_activation= 'sig', L= 1, train_size = 1,
                           S =1, t =2, rf = 0.04861)
      
      MLP_HC_Delta = Derivative_HC2(X = (as.matrix(cbind(S0/K, TimeToExpiry/252))[,,drop = F]), 
                                    C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                    W_C = model$W_C, hidden = c(4), res = m1, activation = 'sig', 
                                    output_activation= 'sig', S =1, L = 1)
      
      BSDelta = pnorm(d1)
     
      for(t in 1:(TimeToExpiry-1)){
        MLP_HC_Delta_prev = MLP_HC_Delta$call_derivative * K*(m1$chat>0)
        MLP_HC_Delta_call_star_prev = MLP_HC_Delta_call_star$Delta * K
        BSDelta_prev = BSDelta
  
        
        m1 = forwardprop_HC2(X =(as.matrix(cbind(1,(data$StockPrice[t]), data$TimeTillExpiry[t]))[,,drop = F]), 
                             C = as.matrix(1), P  = as.matrix(1),
                             W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                             output_activation= 'sig', L= 1, train_size = nrow(data),
                             S =1, t =2, rf = 0.04861)
        
        MLP_HC_Delta = Derivative_HC2(X = (as.matrix(cbind((data$StockPrice[t]), data$TimeTillExpiry[t]))[,,drop = F]),
                                      C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                      W_C = model$W_C, hidden = c(4), res = m1, activation = 'sig', 
                                      output_activation= 'sig', S =1, L = 1)
        
        m1_call_star = forwardprop(X =(as.matrix(cbind(1,(data$StockPrice[t]), data$TimeTillExpiry[t]))[,,drop = F]), 
                                   Y= as.matrix(1), 
                                   W = model$W_C , hidden = c(4), activation = 'sig', 
                                   output_activation= 'sig', L= 1, train_size = nrow(data))
        
        
        MLP_HC_Delta_call_star = CallDelta(X = (as.matrix(cbind((data$StockPrice[t]), data$TimeTillExpiry[t]))[,,drop = F]),
                                           Weights = model$W_C, hidden = c(4), res = m1_call_star, activation = 'sig', 
                                           output_activation= 'sig', S =1)
        

        #Hedging with HC
        V_S[(t+1),i] = data$StockPrice[t] * MLP_HC_Delta$call_derivative  * K*(m1$chat>0)
        V_C[(t+1),i] = -BlackScholes(S = data$StockPrice[t]*K, K = K, t = data$TimeTillExpiry[t]*252, r = 0.04681, sd = 0.2) 
        V_B[(t+1),i] = exp(0.04681/252)* V_B[t,i] - data$StockPrice[t] * (MLP_HC_Delta$call_derivative*(m1$chat>0)  *K - MLP_HC_Delta_prev)
        V[(t+1),i] = V_S[(t+1),i] + V_B[(t+1),i] + V_C[(t+1),i]
        
        #Hedging without HC
        V_call_star_S[(t+1),i] = data$StockPrice[t] * MLP_HC_Delta_call_star$Delta * K
        V_call_star_B[(t+1),i] = exp(0.04681/252)* V_call_star_B[t,i] - data$StockPrice[t] *
          (MLP_HC_Delta_call_star$Delta  *K - MLP_HC_Delta_call_star_prev)
        V_call_star[(t+1),i]   = V_call_star_S[(t+1),i] + V_C[(t+1),i] + V_call_star_B[(t+1),i]
        
        
        
        d1 = 1/(0.2*sqrt(data$TimeTillExpiry[t]))*(log(data$StockPrice[t]) +
                                                     ((0.04681 + (0.2^2)/2)*(data$TimeTillExpiry[t])))
        
        BSDelta = pnorm(d1)
        V_bs_S[(t+1),i] = data$StockPrice[t]* BSDelta * K
        V_bs_B[(t+1),i] = exp(0.04681/252)* V_bs_B[t,i] - data$StockPrice[t]*K*(BSDelta - BSDelta_prev)
        V_bs[(t+1),i] = V_bs_S[(t+1),i] + V_bs_B[(t+1),i] + V_C[(t+1),i]
        
      }
    }
    k = K/5 - 7
    MLP_HC4_HedgeError[[k]][[m]] = V
    MLP_HC4_HedgeError_call_star[[k]][[m]] = V_call_star
    BSHedgeError[[k]] = V_bs
    print(m)
    
  }
  
}


MLP_HC4_HedgeError_ = MLP_HC4_HedgeError
save(MLP_HC4_HedgeError_63, file = paste('mlpHC4_model_001_63_', Sys.Date(), ".RData", sep = ''))
save(BSHedgeError, file = paste('BSHedgeError_63_', Sys.Date(), ".RData", sep = ''))



mean(abs(MLP_HC4_HedgeError[[3]][[1]][TimeToExpiry,]))
mean(abs(BSHedgeError[[3]][TimeToExpiry,]))
length(which(abs(MLP_HC4_HedgeError[[3]][[1]][TimeToExpiry,]) <abs(BSHedgeError[[3]][TimeToExpiry,])))/500
######TAB III

tab_III = matrix(NA, nrow = 11, ncol = 5)
tab_III[11,] = abs(V_bs[TimeToExpiry,c(100,200,300,400,500)])

for (i in 1:10){
  tab_III[i,] = abs(MLP_HC4_HedgeError[[3]][[i]][TimeToExpiry,c(100,200,300,400,500)])
}
#View(round(tab_III, 4))
write.csv(round(tab_III,4), file = paste('tab_III_MLP_HC4_63_', Sys.Date(),'.csv', sep = ''))

######TAB IV

tab_IV = c()

for (i in 1:10){
  one = length(which(abs(MLP_HC4_HedgeError[[3]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[3]][TimeToExpiry,])))/500
  two = paste("(", round(((1- one)* one / sqrt(500)),3), ")", sep = "")
  tab_IV = c(tab_IV, one, two)
}


#View(tab_IV)
write.csv(tab_IV, file = paste('tab_IV_MLP_HC4_63_', Sys.Date(),'.csv', sep = ''))

######TAB V

tab_V = as.data.frame(matrix(NA, nrow = 4, ncol = 5))
names(tab_V)  = c('K = 40','K = 45','K = 50','K = 55','K = 60')
row.names(tab_V) = c('Mean' , '(SE)', 'Minimum', 'Maximum')

for (k in 1:5){
  HedgeErrFrac = numeric(10)
  for (i in 1:10){
    HedgeErrFrac[i]= length(which(abs(MLP_HC4_HedgeError[[k]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[k]][TimeToExpiry,])))/500
  }
  tab_V[1,k] = round(mean(HedgeErrFrac),4)
  sd = tab_V[1,k]*(1-tab_V[1,k])/sqrt(500)
  tab_V[2,k] = paste("(", round(sd, 3), ")", sep = "")
  tab_V[3,k] = round(min(HedgeErrFrac),4)
  tab_V[4,k] = round(max(HedgeErrFrac),4)
}


write.csv(tab_V, file = paste('tab_V_MLP_HC4_63_', TimeToExpiry ,"_", Sys.Date(),'.csv', sep = ''))


#Tab VI
tab_VI = as.data.frame(matrix(NA, nrow = 5, ncol = 5))
names(tab_VI)  = c('K = 40','K = 45','K = 50','K = 55','K = 60')
row.names(tab_VI) = c('Mean' , '(SE)', 'Minimum', 'Maximum')

for (k in 1:5){
  PreditionErr = numeric(10)
  PreditionErrBS= sqrt(mean(BSHedgeError[[k]][TimeToExpiry,]^2) + var(BSHedgeError[[k]][TimeToExpiry,]))
  for (i in 1:10){
    PreditionErr[i]= sqrt(mean(MLP_HC4_HedgeError[[k]][[i]][TimeToExpiry,]^2) + var(MLP_HC4_HedgeError[[k]][[i]][TimeToExpiry,]))
  }
  tab_VI[1,k] = round(mean(PreditionErr),4)
  sd = sd(PreditionErr)
  tab_VI[2,k] = paste("(", round(sd, 3), ")", sep = "")
  tab_VI[3,k] = round(min(PreditionErr),4)
  tab_VI[4,k] = round(max(PreditionErr),4)
  tab_VI[5,k] = round(PreditionErrBS,4)
  
}

#View(tab_VI)
write.csv(tab_VI, file = paste('tab_VI_HC4_63_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))


