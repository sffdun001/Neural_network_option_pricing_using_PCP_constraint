#Duncan Saffy 2019 Dissertation
#Option pricing model

#002 accounts for the model applying a 


rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/MLP Hutch")
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
#Pricing Error Calculations
##############################################################################

MLPPriceError = vector(mode = 'list', length = 10)

for (m in 1:10){
  model <-  mlp_sig_model[[m]]
  pe = matrix(NA, nrow = 500, ncol = 3)
  for(i in 1:500){
    data <-  TEST_DATA[[i]]%>%filter(TerminalDate == 63 & Strike == 50)%>%
      arrange(Day)
    data$StockPrice = data$StockPrice /data$Strike
    data$TimeTillExpiry = data$TimeTillExpiry / 252
    
    m1 = forwardprop(X = as.matrix(cbind(1,data[,c('StockPrice','TimeTillExpiry')])), 
                          Y= as.matrix(data[,'CallPrice']), 
                          W = model$W , hidden = c(4), activation = 'sig', 
                          output_activation= 'sig', L= 1, train_size = nrow(data))
    
    MSE = mean((m1$Yhat-data$CallPrice /data$Strike)^2)
    MAE = mean(abs(m1$Yat-data$CallPrice /data$Strike))
    R2 = 1 - sum((m1$Yhat-data$CallPrice  /data$Strike)^2) / 
      sum((data$CallPrice /data$Strike - mean(data$CallPrice /data$Strike))^2)
    pe[i,] = c(MSE,MAE,R2)
  }
  MLPPriceError[[m]] = pe
}

r2mean = c()
r2min = c()
r2max = c()
for (i in 1:10){
  r2mean = c(r2mean, mean(MLPPriceError[[i]][,3]))
  r2min = c(r2min, min(MLPPriceError[[i]][,3]))
  r2max = c(r2max, max(MLPPriceError[[i]][,3]))
}

price_error = as.data.frame(matrix(NA, nrow = 4, ncol = 1))
colnames(price_error) = "MLP"
price_error[1,1] = paste(round(mean(r2mean)*100,2), '%',sep = "")
price_error[2,1] = paste(round(sd(r2mean)*100,2), '%',sep = "")
price_error[3,1] = paste(round(max(r2max)*100,2), '%',sep = "")
price_error[4,1] = paste(round(min(r2min)*100,2), '%',sep = "")
#min(r2min)
write.csv(price_error, file = paste('mlp_price_error_table_001_', Sys.Date(), '.csv', sep = '')) 

save(MLPPriceError, file = paste('mlp_price_error_001_', Sys.Date(), '.csv', sep = ''))



##############################################################################
#Hedging Error Calculations
##############################################################################

BSHedgeError = vector(mode = 'list', length = 5)
MLPHedgeError = vector(mode = 'list', length = 5)



#Calculating 3 month hedging error.
TimeToExpiry = 126
for (K in c(40, 45, 50, 55, 60)){
  #K = 50
  for (m in 1:10){
   #m =1
   #K = 50
    model <-  mlp_sig_model[[m]]
    
    S0 =50
    N = 500
    V_S = matrix(NA, nrow = 253, ncol = N)
    V_C = matrix(NA, nrow = 253, ncol = N)
    V_B = matrix(NA, nrow = 253, ncol = N)
    V = matrix(NA, nrow = 253, ncol = N)
    
    V_bs_S = matrix(NA, nrow = 253, ncol = N)
    #V_bs_C = matrix(NA, nrow = 63, ncol = N)
    V_bs_B = matrix(NA, nrow = 253, ncol = N)
    V_bs = matrix(NA, nrow = 253, ncol = N)

    
    
    m1 = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/252))[,,drop = F]), 
                     Y= as.matrix(1), 
                     W = model$W , hidden = c(4), activation = 'sig', 
                     output_activation= 'sig', L= 1, train_size = 1)
    
    
    MLPDelta = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/252))[,,drop = F]), 
                         Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                         output_activation= 'sig', S =1)
    
    if (MLPDelta$Delta <0){print(MLPDelta$Delta)}
    
    V_S[1,] = S0*MLPDelta$Delta
    V_C[1,] = -BlackScholes(S = rep(S0, N), K = K, t = TimeToExpiry, r = 0.04681, sd = 0.2) 
    V_B[1,] = -(V_S[1,] + V_C[1,]) 
    V[1,]   = V_S[1,] + V_C[1,] + V_B[1,]
    
    d1 = 1/(0.2*sqrt(TimeToExpiry/252))*(log(S0/K) +
                                 ((0.04681 + (0.2^2)/2)*(TimeToExpiry/252)))
    BSDelta = pnorm(d1)
    V_bs_S[1,] = S0*BSDelta
    V_bs_B[1,] = -(V_bs_S[1,] + V_C[1,]) 
    V_bs[1,]   = V_bs_S[1,] + V_C[1,] + V_bs_B[1,]
    
    for (i in 1:500){
      data <-  TEST_DATA[[i]]%>%filter(TerminalDate == TimeToExpiry & Strike == K)%>%
        arrange(Day)
      data$StockPrice = data$StockPrice / K 
      data$TimeTillExpiry = data$TimeTillExpiry / 252
      
      d1 = 1/(0.2*sqrt(TimeToExpiry/252))*(log(S0/K) +
                                   ((0.04681 + (0.2^2)/2)*(TimeToExpiry/252)))
      m1 = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/252))[,,drop = F]), 
                       Y= as.matrix(1), 
                       W = model$W , hidden = c(4), activation = 'sig', 
                       output_activation= 'sig', L= 1, train_size = 1)
      BSDelta = pnorm(d1)
      MLPDelta = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/252))[,,drop = F]), 
                           Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                           output_activation= 'sig', S =1)
      
      for(t in 1:(TimeToExpiry-1)){
        MLPDelta_prev = MLPDelta$Delta *K
        BSDelta_prev = BSDelta
        
        m1 = forwardprop(X =(as.matrix(cbind(1,(data$StockPrice[t]), data$TimeTillExpiry[t]))[,,drop = F]), 
                         Y= as.matrix(1), 
                         W = model$W , hidden = c(4), activation = 'sig', 
                         output_activation= 'sig', L= 1, train_size = nrow(data))
        
        MLPDelta = CallDelta(X = (as.matrix(cbind((data$StockPrice[t]), data$TimeTillExpiry[t]))[,,drop = F]), 
                             Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                             output_activation= 'sig', S =1)
        if (MLPDelta$Delta <0){print(MLPDelta$Delta)}
        
        V_S[(t+1),i] = data$StockPrice[t] * MLPDelta$Delta * K
        V_C[(t+1),i] = -BlackScholes(S = data$StockPrice[t]*K, K = K, t = data$TimeTillExpiry[t]*252, r = 0.04681, sd = 0.2) 
        V_B[(t+1),i] = exp(0.04681/252)* V_B[t,i] - data$StockPrice[t] * (MLPDelta$Delta *K+ - MLPDelta_prev)
        V[(t+1),i] = V_S[(t+1),i] + V_B[(t+1),i] + V_C[(t+1),i]
        
        d1 = 1/(0.2*sqrt(data$TimeTillExpiry[t]))*(log(data$StockPrice[t]) +
                                                         ((0.04681 + (0.2^2)/2)*(data$TimeTillExpiry[t])))
        
        BSDelta = pnorm(d1)
        V_bs_S[(t+1),i] = data$StockPrice[t]* BSDelta * K
        V_bs_B[(t+1),i] = exp(0.04681/252)* V_bs_B[t,i] - data$StockPrice[t]*K*(BSDelta - BSDelta_prev)
        V_bs[(t+1),i] = V_bs_S[(t+1),i] + V_bs_B[(t+1),i] + V_C[(t+1),i]
        
      }
    }
    k = K/5 - 7
    MLPHedgeError[[k]][[m]] = V
    BSHedgeError[[k]]= V_bs
    print(m)
    
  }
  
}



mean(abs(MLPHedgeError[[3]][[1]][63,]))
mean(abs(BSHedgeError[[3]][63,]))
length(which(abs(MLPHedgeError[[3]][[1]][63,]) <abs(BSHedgeError[[3]][63,])))/500


hist(MLPHedgeError[[3]][[1]][63,])
hist(BSHedgeError[[3]][63,])

hist(V[63,], breaks = 15, xlim = c(-1,1), 
     xlab = "Terminal Hedging Error", 
     main = "Distribution of B-S Network Hedging Error") #, breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
hist(V_bs[63,], breaks = 15, xlim = c(-1,1),
     xlab = "Terminal Hedging Error", 
     main = "Distribution of B-S Hedging Error")



?hist

######TAB III

tab_III = matrix(NA, nrow = 11, ncol = 5)
tab_III[11,] = abs(V_bs[TimeToExpiry,c(100,200,300,400,500)])

for (i in 1:10){
  tab_III[i,] = abs(MLPHedgeError[[3]][[i]][TimeToExpiry,c(100,200,300,400,500)])
}
View(round(tab_III, 4))
write.csv(round(tab_III,4), file = paste('tab_III_MLP_',TimeToExpiry,'_', Sys.Date(),'.csv', sep = ''))

######TAB IV

tab_IV = c()

for (i in 1:10){
  one = length(which(abs(MLPHedgeError[[3]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[3]][TimeToExpiry,])))/500
  two = paste("(", round(((1- one)* one / sqrt(500)),3), ")", sep = "")
  tab_IV = c(tab_IV, one, two)
}


View(tab_IV)
write.csv(tab_IV, file = paste('tab_IV_MLP_',TimeToExpiry, "_" Sys.Date(),'.csv', sep = ''))

######TAB V

tab_V = as.data.frame(matrix(NA, nrow = 4, ncol = 5))
names(tab_V)  = c('K = 40','K = 45','K = 50','K = 55','K = 60')
row.names(tab_V) = c('Mean' , '(SE)', 'Minimum', 'Maximum')

for (k in 1:5){
  HedgeErrFrac = numeric(10)
  for (i in 1:10){
    HedgeErrFrac[i]= length(which(abs(MLPHedgeError[[k]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[k]][TimeToExpiry,])))/500
  }
  tab_V[1,k] = round(mean(HedgeErrFrac),4)
  sd = tab_V[1,k]*(1-tab_V[1,k])/sqrt(500)
  tab_V[2,k] = paste("(", round(sd, 3), ")", sep = "")
  tab_V[3,k] = round(min(HedgeErrFrac),4)
  tab_V[4,k] = round(max(HedgeErrFrac),4)
}

View(tab_V)
write.csv(tab_V, file = paste('tab_V_MLP_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))


######TAB VI

tab_VI = as.data.frame(matrix(NA, nrow = 4, ncol = 5))
names(tab_VI)  = c('K = 40','K = 45','K = 50','K = 55','K = 60')
row.names(tab_VI) = c('Mean' , '(SE)', 'Minimum', 'Maximum')

for (k in 1:5){
  PreditionErr = numeric(10)
  #PreditionErrBS= sqrt(mean(BSHedgeError[[k]][TimeToExpiry,]^2) + var(BSHedgeError[[k]][TimeToExpiry,]))
  for (i in 1:10){
    PreditionErr[i]= sqrt(mean(MLPHedgeError[[k]][[i]][TimeToExpiry,]^2) + var(MLPHedgeError[[k]][[i]][TimeToExpiry,]))
  }
  tab_VI[1,k] = round(mean(PreditionErr),4)
  sd = sd(PreditionErr)
  tab_VI[2,k] = paste("(", round(sd, 3), ")", sep = "")
  tab_VI[3,k] = round(min(PreditionErr),4)
  tab_VI[4,k] = round(max(PreditionErr),4)
  
  
}

#View(tab_VI)
write.csv(tab_VI, file = paste('tab_VI_MLP_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))


