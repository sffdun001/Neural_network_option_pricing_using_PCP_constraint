#Duncan Saffy 2020 Dissertation

rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/BlackTest")
load('mlpHC4_black_001_2020-03-29.RData')
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

for (K in c(45000, 47500, 50000, 52500, 55000)){
  MLPPriceError = vector(mode = 'list', length = 10)
  
  for (m in 1:10){
    model <-  mlpHC2_black_model[[m]]
    pe = matrix(NA, nrow = 500, ncol = 3)
    for(i in 1:500){
      data1 <-  df_black_test[[i]][which(df_black_test[[i]]$day%in%c(1:91)),]
      data = data1[which(data1$Strike == K),]
      data$StockPrice = data$StockPrice /data$Strike
      data$TimeTillExpiry = data$TimeTillExpiry / 365
      data$CallPrice = data$CallPrice /data$Strike
      
      
      m1 = forwardprop_HC2(X = as.matrix(cbind(1,data[,c('StockPrice','TimeTillExpiry')])), 
                           C = as.matrix(data[,'CallPrice']), P  = as.matrix(data[,'PutPrice']),
                           W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                           output_activation= 'sig', L= 1, train_size = nrow(data),
                           S =1, t =2, rf = 0.0687153)
      
      MSE = mean((m1$chat*(m1$chat>0)-data$CallPrice)^2)
      MAE = mean(abs(m1$chat*(m1$chat>0)-data$CallPrice))
      R2 = 1 - sum((m1$chat*(m1$chat>0)-data$CallPrice)^2) / 
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
  colnames(price_error) = "MLP HC4"
  price_error[1,1] = paste(round(mean(r2mean)*100,2), '%',sep = "")
  price_error[2,1] = paste(round(sd(r2mean)*100,2), '%',sep = "")
  price_error[3,1] = paste(round(max(r2max)*100,2), '%',sep = "")
  price_error[4,1] = paste(round(min(r2min)*100,2), '%',sep = "")
  #min(r2min)
  write.csv(price_error, file = paste('mlpHC4_price_error_table_001_',K,'_', Sys.Date(), '.csv', sep = '')) 
  
  save(MLPPriceError, file = paste('mlpHC4_price_error_001_',K,'_', Sys.Date(), '.csv', sep = ''))

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
  
    model <-  mlpHC2_black_model[[m]]
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


    
    m1 = forwardprop_HC2(X = as.matrix(cbind(1, S0/K, TimeToExpiry/365)), 
                         C = as.matrix(1), P  = as.matrix(1),
                         W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                         output_activation= 'sig', L= 1, train_size = nrow(data),
                         S =1, t =2, rf = 0.0687153)
    
    
    
    MLPDelta = Derivative_HC2(X = (as.matrix(cbind((S0/K), TimeToExpiry/365))[,,drop = F]), 
                                  C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                  W_C = model$W_C, hidden = c(4), res = m1, activation = 'sig', 
                                  output_activation= 'sig', S =1, L = 1)
    
    if (MLPDelta$call_derivative <0){print(MLPDelta$call_derivative)}
    
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
     
      
      m1 = forwardprop_HC2(X = as.matrix(cbind(1, S0/K, TimeToExpiry/365))[,,drop = F], 
                           C = as.matrix(1), P  = as.matrix(1),
                           W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                           output_activation= 'sig', L= 1, train_size = nrow(data),
                           S =1, t =2, rf = 0.0687153)
      
      BSDelta = pnorm(d1)*exp(-0.0687153*TimeToExpiry/365)
     
      
      MLPDelta =  Derivative_HC2(X = (as.matrix(cbind((S0/K), TimeToExpiry/365))[,,drop = F]), 
                                 C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                 W_C = model$W_C, hidden = c(4), res = m1, activation = 'sig', 
                                 output_activation= 'sig', S =1, L = 1)
      
      for(t in 1:(TimeToExpiry-1)){
        #t =1
        MLPDelta_prev = MLPDelta$call_derivative*(m1$chat>0)
        BSDelta_prev = BSDelta
        
     
        m1 = forwardprop_HC2(X = as.matrix(cbind(1, (data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F], 
                             C = as.matrix(1), P  = as.matrix(1),
                             W_C = model$W_C , W_P = model$W_P , hidden = c(4), activation = 'sig', 
                             output_activation= 'sig', L= 1, train_size = nrow(data),
                             S =1, t =2, rf = 0.0687153)
        
        MLPDelta =  Derivative_HC2(X = (as.matrix(cbind((data$StockPrice[t+1]), data$TimeTillExpiry[t+1]))[,,drop = F]), 
                                   C = as.matrix(1), P  = as.matrix(1), W_P = model$W_P,
                                   W_C = model$W_C, hidden = c(4), res = m1, activation = 'sig', 
                                   output_activation= 'sig', S =1, L = 1)
        
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

mean(MLPHedgeError[[3]][[1]][90,])
mean(abs(MLPHedgeError[[3]][[1]][90,]))
mean(abs(BSHedgeError[[3]][90,]))
length(which(abs(MLPHedgeError[[3]][[1]][90,]) <abs(BSHedgeError[[3]][90,])))/500

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
write.csv(round(tab_III,4), file = paste('tab_III_black_HC4_', Sys.Date(),'.csv', sep = ''))

######TAB IV

tab_IV = c()

for (i in 1:10){
  one = length(which(abs(MLPHedgeError[[3]][[i]][TimeToExpiry,]) < abs(BSHedgeError[[3]][TimeToExpiry,])))/500
  two = paste("(", round(((1- one)* one / sqrt(500)),3), ")", sep = "")
  tab_IV = c(tab_IV, one, two)
}


#View(tab_IV)
write.csv(tab_IV, file = paste('tab_IV_black_HC4_', Sys.Date(),'.csv', sep = ''))

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
write.csv(tab_V, file = paste('tab_V_black_HC4_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))


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
write.csv(tab_VI, file = paste('tab_VI_HC4_Black_', TimeToExpiry, "_", Sys.Date(),'.csv', sep = ''))


