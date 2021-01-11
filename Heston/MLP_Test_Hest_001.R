#Duncan Saffy 2019 Dissertation
#Option pricing model

#Version 003 is a neural netork that prices puts and calls and enforces put call parity as well as includes cross validation code
#Version 004 is the same as 003 but now we try use a softplus and RelU function in the output layer to avoid negative prices.
#version 005 only prices puts or calls. Not both.
#Version 007 only prices puts or calls and 



rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/HestonTest")
load('mlp_hest_001_2020-03-12.RData')
load('Heston_test_final_2020-03-28.RData')
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

Blacks <- function(F0, K, t, r, sd, call = TRUE, DT = 1/252) #So this apparently works
{ 
  d1 <- 1/(sd*sqrt(t*DT)) * (log(F0/K) +(((sd^2)/2)*(t*DT)))
  d2 <- d1 - sd*sqrt(t*DT)
  C <- exp(-r*t*DT)*(F0*pnorm(d1) - K*pnorm(d2))
  if(call){
    return(C) #returning the call price
  }
  else {
    return(C + K*exp(-r*t) - F0*exp(-r*t)) #returning the put price
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

cf_heston <- function(u, lnS, r, d, t, V0, theta, kappa, omega, rho)
{
  alfa = -0.5*(u*u + u*1i)
  beta = kappa-rho*omega*u*1i
  omega2 = omega*omega
  
  
  D = sqrt(beta*beta-4*alfa*(omega2/2))
  
  bD = beta - D
  eDt = exp(-D*t)
  
  G = bD / (beta  + D) 
  B = (bD / omega2) * ((1 - eDt) / (1 - G*eDt))
  psi = (G * eDt - 1)/ (G - 1) 
  A = ((kappa*theta) / omega2) * (bD*t - 2*log(psi))
  
  y = exp(A + B*V0 + 1i*u*(lnS)) # + (r-d)*t   ----- - theta*t*u^2/2
  return(y)
}


cf_b <- function(u, lnS, r, d, t, sig)
{
  u1 = lnS - ((sig^2) / 2) * t
  u2 = -0.5 * (sig^2)*t *(u^2)
  phi = exp(1i*u*u1)*exp(u2)
  return(phi)
}


psi <- function(model, v, alpha, lnS, r, d, t, ...)
{
  cf = get(model)
  cf_result = cf(u = (v - (alpha + 1) * 1i), lnS, r, d, t, ...) 
  ret = cf_result/ 
    (alpha^2 + alpha - v^2 + 1i * (2*alpha +1) * v )
  return(ret)
}


call_price_cm_vectorised <- function(model, S ,K ,r ,d, t, ...)
{
  lnS = log(S)
  lnK = log(K)
  
  alpha = 0.5
  
  FFT_N = 2^18
  FFT_eta = 0.05
  
  FFT_lambda = (2*pi) / (FFT_N * FFT_eta)
  FFT_b = (FFT_N *FFT_lambda) / 2
  
  uvec = 1:FFT_N
  ku = -FFT_b +FFT_lambda*(uvec - 1)
  jvec = 1:FFT_N
  vj = (uvec - 1) * FFT_eta
  
  tmp = exp(-r*t) * psi(model, v = vj, alpha = alpha, lnS, r, d, t, ...) * exp(1i * vj * (FFT_b)) * FFT_eta
  tmp = (tmp / 3) * (3 + (-1)^jvec - ((jvec -1)==0))
  cpvec = Re(exp(-alpha * ku) * fft(tmp) / pi)
  
  indexOfStrike = floor((lnK + FFT_b)/FFT_lambda + 1)
  
  xp = cbind(ku[indexOfStrike], ku[indexOfStrike + 1])
  yp = cbind(cpvec[indexOfStrike], cpvec[indexOfStrike + 1])
  
  deltavec = Re(exp(-alpha*ku)*fft((1i*vj + (alpha + 1)) * tmp)/S/pi)
  ydelta = cbind(deltavec[indexOfStrike], deltavec[indexOfStrike + 1])
  delta_fft = Re(interpl_vec(x = xp, y = ydelta, xi = lnK))
  
  call_price = Re(interpl_vec(x = xp, y = yp, xi = lnK))
  
  return(list(call_price = call_price, y = cpvec, x = ku, delta = delta_fft))
}

interpl_vec <- function(x, y, xi)
{
  x0 = x[,1]; x1 = x[,2]
  y0 = y[,1]
  y1 = y[,2]
  
  m = (y1 - y0)/(x1 -x0)
  
  b = y0 - m*x0
  
  yi = m * xi + b
  return(yi)
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

mean(r2mean)
max(r2max)
min(r2min)

save(MLPPriceError, file = paste('mlp_price_error_001_', Sys.Date(), '.csv', sep = ''))



##############################################################################
#Hedging Error Calculations
##############################################################################


params = read.csv('param_90.csv')
model_parameters = c(exp(params[1,2]), exp(params[2,2]), exp(params[3,2]), tanh(params[4,2]),
                   exp(params[5,2]) )

kappa = model_parameters[1]
theta = model_parameters[2]
omega = model_parameters[3]
rho = model_parameters[4]
V0 = model_parameters[5]

BSHedgeError = vector(mode = 'list', length = 5)
MLPHedgeError = vector(mode = 'list', length = 5)



#Calculating 3 month hedging error.
TimeToExpiry = 90
for (K in c(45000, 47500, 50000, 52500, 55000)){
  #K = 50
  for (m in 1:10){
    K = 50000
    m = 1
    model <-  mlp_hest_model[[m]]
    
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
    
    Ct = call_price_cm_vectorised(model = 'cf_heston', S = S0 ,K = K ,r = 0.0687153 ,d = 0, 
                                  t = TimeToExpiry / 365, kappa = kappa, V0 = V0, theta = theta,
                                  rho = rho, omega = omega)
    
    
    
    V_S[1,] = S0*MLPDelta$Delta
    V_C[1,] = -Ct$call_price
    V_B[1,] = -(V_S[1,] + V_C[1,]) 
    V[1,]   = V_S[1,] + V_C[1,] + V_B[1,]
    

    BSDelta = Ct$delta
    V_bs_S[1,] = S0*BSDelta
    V_bs_B[1,] = -(V_bs_S[1,] + V_C[1,]) 
    V_bs[1,]   = V_bs_S[1,] + V_C[1,] + V_bs_B[1,]
    
    for (i in 1:500){
      data <-  df_hest_test[[i]]%>%filter(day < 92 & Strike == K)%>%
        arrange(day)
      data$StockPrice = data$StockPrice / K 
      data$TimeTillExpiry = data$TimeTillExpiry / 365
      
      
      m1 = forwardprop(X =(as.matrix(cbind(1, S0/K, TimeToExpiry/365))[,,drop = F]), 
                       Y= as.matrix(1), 
                       W = model$W , hidden = c(4), activation = 'sig', 
                       output_activation= 'sig', L= 1, train_size = 1)
      
      Ct = call_price_cm_vectorised(model = 'cf_heston', S = data$StockPrice[1] ,K = K ,r = 0.0687153 ,d = 0, 
                                     t = TimeToExpiry / 365, kappa = kappa, V0 = V0, theta = theta,
                                     rho = rho, omega = omega)
       
      
      MLPDelta = CallDelta(X = (as.matrix(cbind(S0/K, TimeToExpiry/365))[,,drop = F]), 
                           Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                           output_activation= 'sig', S =1)
      
      BSDelta =  Ct$delta
      
      for(t in 1:(TimeToExpiry-1)){
        MLPDelta_prev = MLPDelta$Delta *K
        BSDelta_prev = BSDelta
        
        m1 = forwardprop(X =(as.matrix(cbind(1,(data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
                         Y= as.matrix(1), 
                         W = model$W , hidden = c(4), activation = 'sig', 
                         output_activation= 'sig', L= 1, train_size = nrow(data))
        
        MLPDelta = CallDelta(X = (as.matrix(cbind((data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
                             Weights = model$W, hidden = c(4), res = m1, activation = 'sig', 
                             output_activation= 'sig', S =1)
        
        Ct = call_price_cm_vectorised(model = 'cf_heston', S = data$StockPrice[t+1] ,K = K ,r = 0.0687153 ,d = 0, 
                                      t = data$TimeTillExpiry[t+1]/ 365, kappa = kappa, V0 = V0, theta = theta,
                                      rho = rho, omega = omega)
        
        V_S[(t+1),i] = data$StockPrice[t+1] * MLPDelta$Delta * K
        V_C[(t+1),i] = -Ct$call_price
        V_B[(t+1),i] = exp(0.0687153/365)* V_B[t,i] - data$StockPrice[t+1] * (MLPDelta$Delta *K - MLPDelta_prev)
        V[(t+1),i] = V_S[(t+1),i] + V_B[(t+1),i] + V_C[(t+1),i]
        
        
        BSDelta = Ct$delta
        
        V_bs_S[(t+1),i] = data$StockPrice[t+1]* BSDelta * K
        V_bs_B[(t+1),i] = exp(0.0687153/365)* V_bs_B[t,i] - data$StockPrice[t]*K*(BSDelta - BSDelta_prev)
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


hist(V[90,])
hist(V_bs[1:90,1:170])
######TAB III

tab_III = matrix(NA, nrow = 11, ncol = 5)
tab_III[11,] = abs(V_bs[TimeToExpiry,c(100,200,300,400,500)])

for (i in 1:10){
  tab_III[i,] = abs(MLPHedgeError[[3]][[i]][TimeToExpiry,c(100,200,300,400,500)])
}
View(round(tab_III, 4))
write.csv(round(tab_III,4), file = paste('tab_III_MLP_', Sys.Date(),'.csv', sep = ''))

######TAB IV

tab_IV = c()

for (i in 1:10){
  one = length(which(abs(MLPHedgeError[[3]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[3]][TimeToExpiry,])))/500
  two = paste("(", round(((1- one)* one / sqrt(500)),3), ")", sep = "")
  tab_IV = c(tab_IV, one, two)
}


View(tab_IV)
write.csv(tab_IV, file = paste('tab_IV_MLP_', Sys.Date(),'.csv', sep = ''))

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

View(tab_V)
write.csv(tab_V, file = paste('tab_V_MLP_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))




