#Duncan Saffy 2019 Dissertation
#Option pricing model

#Hedging on Futures done correctly 


rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/BlackTest")
load('mlp_black_001_2020-03-11.RData')
load('Black_test_final_2020-03-28.RData')
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

OptFunction = function(sigma, F0, C, t, r)
{
  ans = abs(C- Blacks(F0,  K,  t, r, sd =sigma))
  return(ans)
}


##############################################################################
#Pricing Error Calculations
##############################################################################

for (K in c(45000, 47500, 50000, 52500, 55000)){
  #K = 45000
  MLPPriceError = vector(mode = 'list', length = 10)
  
  for (m in 1:10){
    #m = 1
    model <-  mlp_black_model[[m]]
    pe = matrix(NA, nrow = 500, ncol = 3)
    for(i in 1:500){
     # i =1
      data1 <-  df_black_test[[i]][which(df_black_test[[i]]$day%in%c(1:91)),]
      data = data1[which(data1$Strike == K),]
      data$StockPrice = data$StockPrice /data$Strike
      data$TimeTillExpiry = data$TimeTillExpiry / 365
      data$CallPrice = data$CallPrice /data$Strike
      
      m1 = forwardprop(X = as.matrix(cbind(1,data[,c('StockPrice','TimeTillExpiry')])), 
                       Y= as.matrix(data[,'CallPrice']), 
                       W = model$W , hidden = c(4), activation = 'sig', 
                       output_activation= 'sig', L= 1, train_size = nrow(data))
      
      MSE = mean((m1$Yhat-data$CallPrice)^2)
      MAE = mean(abs(m1$Yhat-data$CallPrice))
      R2 = 1 - sum((m1$Yhat-data$CallPrice)^2) / 
        sum((data$CallPrice - mean(data$CallPrice))^2)
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
  write.csv(price_error, file = paste('mlp_price_error_table_001_',K,'_', Sys.Date(), '.csv', sep = '')) 
  
  save(MLPPriceError, file = paste('mlp_price_error_001_',K,'_', Sys.Date(), '.csv', sep = ''))

}





##############################################################################
#Hedging Error Calculations
##############################################################################

BSHedgeError = vector(mode = 'list', length = 5)
MLPHedgeError = vector(mode = 'list', length = 5)



#Calculating 3 month hedging error.
TimeToExpiry = 91
for (K in c(45000, 47500, 50000, 52500, 55000)){
  for (m in 1:10){
  
    model <-  mlp_black_model[[m]]
    data <-  df_black_test[[1]]%>%filter(day < 92 & Strike == K)%>%
      arrange(day)
    
    S0 =50000
    N = 500
    V_S = matrix(NA, nrow = 253, ncol = N)
    V_C = matrix(NA, nrow = 253, ncol = N)
    V_B = matrix(NA, nrow = 253, ncol = N)
    V = matrix(NA, nrow = 253, ncol = N)
    
    V_bs_S = matrix(NA, nrow = 253, ncol = N)
    #V_bs_C = matrix(NA, nrow = 63, ncol = N)
    V_bs_B = matrix(NA, nrow = 253, ncol = N)
    V_bs = matrix(NA, nrow = 253, ncol = N)

    
    
    m1 = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/365))[,,drop = F]), 
                     Y= as.matrix(1), 
                     W = model$W , hidden = c(4), activation = 'sig', 
                     output_activation= 'sig', L= 1, train_size = 1)
    
    
    MLPDelta = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/365))[,,drop = F]), 
                         Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                         output_activation= 'sig', S =1)
    
    if (MLPDelta$Delta <0){print(MLPDelta$Delta)}
    
    V_S[1,] = 0#S0*MLPDelta$Delta
    V_C[1,] = -data$CallPrice[1]
    V_B[1,] = -(V_S[1,] +V_C[1,])####*exp(-0.0687153*TimeToExpiry/365)
    V[1,]   = V_S[1,] + V_C[1,] + V_B[1,] #### *exp(-0.0687153*TimeToExpiry/365)
    
    d1 = (1/(0.173*sqrt(TimeToExpiry/365)))*(log(S0/K)+(0.173^2)*0.5*TimeToExpiry/365)
    BSDelta = pnorm(d1)*exp(-0.0687153*TimeToExpiry/365)
    V_bs_S[1,] = 0 #S0*BSDelta
    V_bs_B[1,] = -(V_bs_S[1,]+ V_C[1,]) ### mexp(-0.0687153*TimeToExpiry/365) 
    V_bs[1,]   = V_bs_S[1,]+ V_C[1,] + V_bs_B[1,] ### *exp(-0.0687153*TimeToExpiry/365) 
    
    for (i in 1:500){
      #i = 1
      data <-  df_black_test[[i]]%>%filter(day < 92 & Strike == K)%>%
        arrange(day)
      data$StockPrice = data$StockPrice / K 
      data$TimeTillExpiry = data$TimeTillExpiry / 365
      
      d1 = 1/(0.173*sqrt(TimeToExpiry/365))*(log(S0/K) +
                                   (((0.173^2)/2)*(TimeToExpiry/365)))
      m1 = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/365))[,,drop = F]), 
                       Y= as.matrix(1), 
                       W = model$W , hidden = c(4), activation = 'sig', 
                       output_activation= 'sig', L= 1, train_size = 1)
      BSDelta = pnorm(d1)*exp(-0.0687153*TimeToExpiry/365)
      MLPDelta = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/365))[,,drop = F]), 
                           Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                           output_activation= 'sig', S =1)
      
      for(t in 1:(TimeToExpiry-1)){
        #t =1
        MLPDelta_prev = MLPDelta$Delta
        BSDelta_prev = BSDelta
        
        m1 = forwardprop(X =(as.matrix(cbind(1,(data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
                         Y= as.matrix(1), 
                         W = model$W , hidden = c(4), activation = 'sig', 
                         output_activation= 'sig', L= 1, train_size = nrow(data))
        
        MLPDelta = CallDelta(X = (as.matrix(cbind((data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
                             Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                             output_activation= 'sig', S =1)
      
        #r = 0.0687153, sd = 0.173
        V_S[(t+1),i] = 0#data$StockPrice[t+1] * K * MLPDelta$Delta
        V_C[(t+1),i] = -data$CallPrice[t+1]
        V_B[(t+1),i] = V_B[t,i]*exp(0.0687153/365) + MLPDelta_prev * K * (data$StockPrice[t+1] - data$StockPrice[t])
        V[(t+1),i] = V_B[(t+1),i] + V_S[(t+1),i] + V_C[(t+1),i]
        
        d1 =  1/(0.173*sqrt(data$TimeTillExpiry[t+1]))*(log(data$StockPrice[t+1]) +
                                                  (((0.173^2)/2)*(data$TimeTillExpiry[t+1])))
        
        BSDelta = pnorm(d1)*exp(-0.0687153*data$TimeTillExpiry[t+1])
        V_bs_S[(t+1),i] =  0 #data$StockPrice[t+1] * K * BSDelta
        V_bs_B[(t+1),i] =  V_bs_B[t,i]*exp(0.0687153/365)  +  BSDelta_prev * K * (data$StockPrice[t+1] - data$StockPrice[t])
        V_bs[(t+1),i] = V_bs_B[(t+1),i] + V_bs_S[(t+1),i] + V_C[(t+1),i]
        
      }
    }
    k = K/2500 - 17
    MLPHedgeError[[k]][[m]] = V
    BSHedgeError[[k]]= V_bs
    print(m)
    
  }
  
}

mean(MLPHedgeError[[3]][[1]][TimeToExpiry,])
mean(abs(MLPHedgeError[[3]][[1]][TimeToExpiry,]))
mean(abs(BSHedgeError[[3]][TimeToExpiry,]))
length(which(abs(MLPHedgeError[[3]][[1]][TimeToExpiry,]) <abs(BSHedgeError[[3]][TimeToExpiry,])))/500

length(which(abs(BSHedgeError[[3]][90,])>500))
hist(V[90,], breaks = 15,
     xlab = "Terminal Hedging Error", 
     main = "Distribution of Black Network Hedging Error") #, breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
hist(V_bs[90,],breaks = 15,
     xlab = "Terminal Hedging Error", 
     main = "Distribution of Black Hedging Error")


plot(x = c(1:90), y = V_bs_B[c(1:90),1], type = 'l')
for (j in 2:10){
  lines(x = c(1:90), y = V_bs_B[c(1:90),j], type = 'l')
}

hist(MLPHedgeError[[3]][[1]][63,])
hist(BSHedgeError[[3]][63,])

######TAB III

tab_III = matrix(NA, nrow = 11, ncol = 5)
tab_III[11,] = abs(V_bs[TimeToExpiry,c(100,200,300,400,500)])

for (i in 1:10){
  tab_III[i,] = abs(MLPHedgeError[[3]][[i]][TimeToExpiry,c(100,200,300,400,500)])
}
#View(round(tab_III, 4))
write.csv(round(tab_III,4), file = paste('tab_III_MLP_Market', Sys.Date(),'.csv', sep = ''))

######TAB IV

tab_IV = c()

for (i in 1:10){
  one = length(which(abs(MLPHedgeError[[3]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[3]][TimeToExpiry,])))/500
  two = paste("(", round(((1- one)* one / sqrt(500)),3), ")", sep = "")
  tab_IV = c(tab_IV, one, two)
}


#View(tab_IV)
write.csv(tab_IV, file = paste('tab_IV_MLP_Market', Sys.Date(),'.csv', sep = ''))

######TAB V

tab_V = as.data.frame(matrix(NA, nrow = 4, ncol = 5))
names(tab_V)  = c('K = 40','K = 45','K = 50','K = 55','K = 60')
row.names(tab_V) = c('Mean' , '(SE)', 'Minimum', 'Maximum')

for (k in 1:5){
  HedgeErrFrac = numeric(10)
  for (i in 1:10){
    HedgeErrFrac[i]= length(which(abs(MLPHedgeError[[k]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[k]][TimeToExpiry,])))/500
  }
  tab_V[1,k] = mean(HedgeErrFrac)
  sd = tab_V[1,k]*(1-tab_V[1,k])/sqrt(500)
  tab_V[2,k] = paste("(", round(sd, 3), ")", sep = "")
  tab_V[3,k] = min(HedgeErrFrac)
  tab_V[4,k] = max(HedgeErrFrac)
}

#View(tab_V)
write.csv(tab_V, file = paste('tab_V_MLP_Market', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))


######TAB VI

tab_VI = as.data.frame(matrix(NA, nrow = 4, ncol = 5))
names(tab_VI)  = c('K = 40','K = 45','K = 50','K = 55','K = 60')
row.names(tab_VI) = c('Mean' , '(SE)', 'Minimum', 'Maximum')

for (k in 1:5){
  PreditionErr = numeric(10)
  #PreditionErrBS= sqrt(mean(BSHedgeError[[k]][TimeToExpiry,]^2) + var(BSHedgeError[[k]][TimeToExpiry,]))
  for (i in 1:10){
    PreditionErr[i]= sqrt(mean(MLPHedgeError[[k]][[i]][TimeToExpiry,]^2) + var(MLPHedgeError[[k]][[i]][TimeToExpiry,]))
    PreditionErr_Black= sqrt(mean(BSHedgeError[[k]][TimeToExpiry,]^2) + var(BSHedgeError[[k]][TimeToExpiry,]))
  }
  tab_VI[1,k] = round(mean(PreditionErr),4)
  sd = sd(PreditionErr)
  tab_VI[2,k] = paste("(", round(sd, 3), ")", sep = "")
  tab_VI[3,k] = round(min(PreditionErr),4)
  tab_VI[4,k] = round(max(PreditionErr),4)
  tab_VI[5,k] = round(PreditionErr_Black,4)
  
}

#View(tab_VI)
write.csv(tab_VI, file = paste('tab_VI_MLP_Black_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))

