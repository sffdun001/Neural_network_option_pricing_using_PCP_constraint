# Heston Simulation
# Duncan Saffy Dissertation
rm(list = ls(all = TRUE))
setwd("C:/Users/dunca/Desktop/Dissertation/GitHub/Data Simulation")
#setwd("~/UCT Academic/Masters/Dissertation/GitHub/BlackTest")

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


Blacks <- function(F0, K, t, r, sd, call = TRUE, DT = 1/252) #So this apparently works
{ 
  d1 <- 1/(sd*sqrt(t*DT)) * (log(F0/K) +(((sd^2)/2)*(t*DT)))
  d2 <- d1 - sd*sqrt(t*DT)
  C <- exp(-r*t*DT)*(F0*pnorm(d1) - K*pnorm(d2))
  if(call){
    return(C) #returning the call price
  }
  else {
    return(C + K - F0) #returning the put price
  }
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


params = read.csv('param_90.csv')


Heston_simulation <- function(Ft0, model_parameters, r, d = 0)
{
  #d = 0
  #r = 0.07
  #model_parameters = initial_params
  #Ft0 = 50000
  
  kappa = model_parameters[1]
  theta = model_parameters[2]
  omega = model_parameters[3]
  rho = model_parameters[4]
  V0 = model_parameters[5]
  
  N <- 91*5
  z1 <-  rnorm(N)
  z2 <- rnorm(N)
  
  w2 <- z1
  w1 <- rho*z1 +  sqrt(1-rho^2)*z2
  
  V = numeric(N); Ft = numeric(N); lnFt = numeric(N);
  Ct = numeric(N); Pt = numeric(N)
  V[1]= V0
  V1 = V0
  lnFt[1] = log(Ft0)
  lnFt1 = log(Ft0)
  Ft[1] = Ft0
  quarter_data = list()
  DT = 1/365

  quarter_data= matrix(NA,nrow = 91*201*5, ncol = 7)
  for (quarter in 1:5){
    if (quarter == 1){
      nearest_strike =(round(exp(lnFt1)/50, 0)*50)
      days = (2:91)
      } else{
      nearest_strike = (round(Ft[(quarter-1)*91]/50, 0)*50)
      days = (quarter-1)*91 + (1:91)
    }
    
    K = seq(from = (nearest_strike - (50*100)), to = (nearest_strike + (50*100)), by = 50)
    
    if(quarter == 1){
      if(sum(days == c(2:91))==90){
        Ct = call_price_cm_vectorised(model = 'cf_heston', S = Ft[1] ,K = K ,r = r ,d = 0, 
                                      t = 90 / 365, kappa = kappa, V0 = V0, theta = theta,
                                      rho = rho, omega = omega)
        Pt = Ct$call_price + K*exp(90*r/365) - Ft[1]*exp((90)*(r-d)/365)
        quarter_data[(1:201),] = cbind(Ft[1], 90, K, V[1], Ct$call_price, Pt, 1)
      }
    } 
    
    for(day in days){
      V[day] = V1 + kappa*(theta -max(0,V[day-1]))*DT + omega*max(0,V[day-1])*w2[day]*sqrt(DT)
      lnFt[day] = lnFt1 - (1/2) * max(0, V[day]) * DT + sqrt(max(0, V[day]) * DT)*w1[day]
      Ft[day] = exp(lnFt[day])
      Ct = call_price_cm_vectorised(model = 'cf_heston', S = Ft[day] ,K = K ,r = r ,d = 0, 
                                         t = (quarter*91 - day) / 365, kappa = kappa, V0 = V0, theta = theta,
                                         rho = rho, omega = omega)
      Pt = Ct$call_price + K*exp((quarter*91 - day)*r/365) - Ft[day]*exp((quarter*91 - day)*(r-d)/365)
      
      quarter_data[(day-1)*201 +(1:201),] = cbind(Ft[day], (quarter*91- day), K, V[day], Ct$call_price, Pt, day)
      V1 = V[day]
      lnFt1 = lnFt[day]
    }
    
    
    Ft[day] = exp((r-d)*0.25)*Ft[quarter*91]
    lnFt[day] = log(Ft[quarter*91])
    lnFt1 = lnFt[day] 
  }
  
  quarter_data = as.data.frame(quarter_data)
  names(quarter_data) = c("StockPrice", 'TimeTillExpiry','Strike', 'Volatility', 'CallPrice', 'PutPrice', 'day')
  
  return(quarter_data)
}


Black_simulation <- function(Ft0, model_parameters, r, d = 0)
{
  #d = 0
  #r = 0.07
  #model_parameters = initial_params
  sig = model_parameters
 
  N <- 91*5
  z1 <-  rnorm(N)

  Ft = numeric(N); lnFt = numeric(N);
  Ct = numeric(N); Pt = numeric(N)
  lnFt[1] = log(Ft0)
  lnFt1 = log(Ft0)
  Ft[1] = Ft0
  quarter_data = list()
  DT = 1/365
  
  quarter_data= matrix(NA,nrow = 91*201*5, ncol = 7)
  for (quarter in 1:5){
    #quarter = 1
    
    if (quarter == 1){
      nearest_strike =(round(exp(lnFt1)/50, 0)*50)
      days = (2:91)
    } else{
      nearest_strike = (round(Ft[(quarter-1)*91]/50, 0)*50)
      days = (quarter-1)*91 + (1:91)
    }
    
    K = seq(from = (nearest_strike - (50*100)), to = (nearest_strike + (50*100)), by = 50)
    
    if(quarter == 1){
      if(sum(days == c(2:91))==90){
        Ct = call_price_cm_vectorised(model = 'cf_b', S = Ft[1] ,K = K ,r = r ,d = 0, 
                                      t = 90 / 365, sig = sig)
        Pt = Ct$call_price + K*exp(-90*r/365) - Ft[1]*exp((90)*(-r)/365)
        quarter_data[(1:201),] = cbind(Ft[1], 90, K, sig^2, Ct$call_price, Pt, 1)
      }
    } 
   
    for(day in days){
      lnFt[day] = lnFt1 - (1/2) * sig^2 * DT + sig*sqrt(DT)*z1[day]
      Ft[day] = exp(lnFt[day])
      Ct = call_price_cm_vectorised(model = 'cf_b', S = Ft[day] ,K = K ,r = r ,d = 0, 
                                    t = (quarter*91 - day) / 365, sig = sig)
      Pt = Ct$call_price + K*exp(-(quarter*91 - day)*r/365) - Ft[day]*exp(-(quarter*91 - day)*(r)/365)
      
      quarter_data[(day-1)*201 +(1:201),] = cbind(Ft[day], (quarter*91- day), K, sig^2, Ct$call_price, Pt, day)
      lnFt1 = lnFt[day]
    }
    
    
    Ft[day] = exp((r-d)*0.25)*Ft[quarter*91]
    lnFt[day] = log(Ft[quarter*91])
    lnFt1 = lnFt[day] 
  }
  
  quarter_data = as.data.frame(quarter_data)
  names(quarter_data) = c("StockPrice", 'TimeTillExpiry','Strike', 'Volatility', 'CallPrice', 'PutPrice', 'day')
  
  return(quarter_data)
}

initial_params = c( exp(params[1,2]), exp(params[2,2]), exp(params[3,2]), tanh(params[4,2]),
                  exp(params[5,2]) )


Option_suffix = c('H8', "M8", "U8", "Z8", 'H9', "M9", "U9")

#
library(dplyr)

data <- c()

for(i in 1:5){
  quarter_data <- read.csv(paste("AI", Option_suffix[i],
                                 "_processed_data_2019-07-22.csv", 
                                 sep = ""), header = T)
  
  
  quarter_data.clean  = quarter_data%>%filter(!is.na(CALL_PRICE) | !is.na(PUT_PRICE)) #removing prices that have not year been created
  
  #setting prices to 0 where there were NA values
  na_call_ind = which(is.na(quarter_data.clean$CALL_PRICE))
  na_put_ind = which(is.na(quarter_data.clean$PUT_PRICE))
  
  quarter_data.clean$CALL_PRICE[na_call_ind] = 0
  quarter_data.clean$PUT_PRICE[na_put_ind] = 0
  
  
  quarter_data.clean <- quarter_data.clean%>%
    select(c("CALL_PRICE", "PUT_PRICE", "Strike", "DAY_TILL_EXPIRY", 
             "AI1_LAST", "JIBA3M", "TOP40_EQY_DVD_YLD_12M"))
  
  
  data = rbind(data, cbind(quarter_data.clean, i))
  #changing so only the last month is selected
}

data = data%>%arrange(DAY_TILL_EXPIRY,AI1_LAST)


sig = 0.173 #the implied volatility of at the money strike of AIM8 options 90 days to expiry, found in HestonMOdel.R file
d = log(1+mean(data$TOP40_EQY_DVD_YLD_12M[which(data$i ==2)])/100)
r = log(1+mean(data$JIBA3M[which(data$i ==2)]/4)/100)*4 #0.0687153

#training data
set.seed(1995)
df_hest_train = list()
df_black_train = list()

for (i in 1:10) { 
  #i =1
  df_hest_train[[i]] = Heston_simulation(Ft0 = 50000, model_parameters = initial_params, r = r, d = d)
  df_black_train[[i]] = Black_simulation(Ft0 = 50000, model_parameters = sig, r = r, d= d)
  print(i)
}

save(df_hest_train, file = paste('Heston_train_001_',Sys.Date(), ".RData", sep = ''))
save(df_black_train, file = paste('Black_train_001_',Sys.Date(), ".RData", sep = ''))



plot(df_hest_train[[1]][(seq(from = 0, to = 201*90*5, by = 201)+1),1], 
     type = 'l', ylim = c(35000,80000), xlab = 'Time',ylab = 'Price')
for(i in 2:10){
  lines(df_hest_train[[i]][(seq(from = 0, to = 201*90*5, by = 201)+1),1])
}
grid()

plot(df_black_train[[1]][(seq(from = 0, to = 201*90*5, by = 201)+1),1], 
     type = 'l', ylim = c(35000,80000), xlab = 'Time',ylab = 'Price')
for(i in 2:10){
  lines(df_black_train[[i]][(seq(from = 0, to = 201*90*5, by = 201)+1),1])
}
grid()




set.seed(2020)
df_hest_test = list()
df_black_test = list()

for (i in 1:500) { 
  df_hest_test[[i]] = Heston_simulation(Ft0 = 50000, model_parameters = initial_params, r = r, d = d)
  df_black_test[[i]] = Black_simulation(Ft0 = 50000, model_parameters = sig, r = r, d= d)
  print(i)
}

save(df_hest_test, file = paste('Heston_test_001_',Sys.Date(), ".RData", sep = ''))
save(df_black_test, file = paste('Black_test_001_',Sys.Date(), ".RData", sep = ''))








