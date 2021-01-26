#Duncan Saffy Dissertation 
#Converting code from Finacial Modelling to R for Heston model

rm(list = ls(all = TRUE))


#Function code
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

cf_bs <- function(u, lnS, r, d, t, sig)
{
  u1 = lnS + (r -d - (sig^2) / 2) * t
  u2 = -0.5 * (sig^2)*t *(u^2)
  phi = exp(1i*u*u1)*exp(u2)
  return(phi)
}


cf_b <- function(u, lnS, r, d, t, sig)
{
  u1 = lnS - ((sig^2) / 2) * t
  u2 = -0.5 * (sig^2)*t *(u^2)
  phi = exp(1i*u*u1)*exp(u2)
  return(phi)
}

interpl <- function(x, y, xi)
{
  x0 = x[1]; x1 = x[2]
  y0 = y[1]; y1 = y[2]
  
  m = (y1 - y0)/(x1 -x0)
  
  b = y0 - m*x0
  
  yi = m * xi + b
  return(yi)
}

psi <- function(model, v, alpha, lnS, r, d, t, ...)
{
  cf = get(model)
  cf_result = cf(u = (v - (alpha + 1) * 1i), lnS, r, d, t, ...) 
  ret = cf_result/ 
    (alpha^2 + alpha - v^2 + 1i * (2*alpha +1) * v )
  return(ret)
}

call_price_cm <- function(model, S ,K ,r ,d, t, ...)
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
  
  xp = c(ku[indexOfStrike], ku[indexOfStrike + 1])
  yp = c(cpvec[indexOfStrike], cpvec[indexOfStrike + 1])
  
  
  deltavec = Re(exp(-alpha*ku)*fft((1i*vj + (alpha + 1)) * tmp)/S/pi)
  delta_fft = Re(interpl(x = xp, y = deltavec[c(indexOfStrike, (indexOfStrike+1))], xi = lnK))
  
  call_price = Re(interpl(x = xp, y = yp, xi = lnK))
  
  return(list(call_price = call_price, y = cpvec, x = ku, delta = delta_fft))
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

densityrecovery <- function(model, a, N, ...)
{
 cf = get(model)
 b = a/N
 u = ((0:(N-1)) - N/2)*b
 h2 = ((-1)^(0:(N-1))) * cf(u, ...)
 g = fft(h2)
 y = Re((-1)^(0:(N-1))* g ) * (b / (2 * pi))
 return(y)
}




#parameters for testing CM and charfunctions
S = 100; r = 0.02;d = 0;sig = 0.2;t = 1;v = (1:(2^18)-1)*0.05
u = v; lnS = log(100);t = 1;V0 = 0.04;theta = 0.04;kappa = 0.1
omega = 0.000001; rho = 0; sig = 0.2; alpha = 0.5

resblack = call_price_cm_vectorised(model = 'cf_b', S = 100, K = c(90,100, 110) ,t = 1 ,r = 0.02 ,d = 0, sig = 0.2)
  
reshest = call_price_cm_vectorised(model = 'cf_heston', S = 100, K = c(90,100, 110) ,t = 1 ,r = 0.02 ,d = 0, 
                               rho = 0, theta = 0.04, omega = 0.000001, kappa = 0.1, V0 = 0.04)

#Checking each model is producing the same results
rbind(resblack$call_price, reshest$call_price)
rbind(resblack$delta,reshest$delta)

Blacks(F0 = 100, K =c(90,100, 110), t = 252, r = 0.02, sd = 0.2)
d1 = 1/(0.2*sqrt(1)) * (log(100/c(90,100, 110)) +(((0.2^2)/2)*(1)))
pnorm(d1)*exp(-0.02)

#############################################################################
#Reading and cleaning real data in Real Data for calibration
#############################################################################
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
head(data)
#View(data)

implied_vol <- function(sig, Callval, S, K, t, r){
  ret  = (Callval - Blacks(F0 = S, K = K, t = t, r = r, sd = sig, call = TRUE, DT = 1))^2
  return(ret)
}
#############################################################################
#Calibrating for 90 days
#############################################################################

heston_call_pricer <- function(x, y, model = 'cf_heston', S ,K ,r ,d, t)
{ # this is the function to input into the calibration function

  kappa = exp(x[1])
  theta = exp(x[2])
  V0 = exp(x[5])
  omega = exp(x[3])
  rho = tanh(x[4])
  
  yhat = call_price_cm_vectorised(model = 'cf_heston', S = S ,K = K, r = r,
                                    d = d ,t = t, V0 = V0,  theta = theta, 
                                    kappa = kappa, omega = omega, rho = rho)
  
  diff = sum(((yhat$call_price - y)/y)^2)
  
  return(diff)
}


ind = which(data$i ==2 & data$DAY_TILL_EXPIRY ==90)
y = data$CALL_PRICE[ind]
S = data$AI1_LAST[ind][1]
K = data$Strike[ind]
t = data$DAY_TILL_EXPIRY[ind][1]/365
r = log((1+(data$JIBA3M[ind][1]/100)/4)^4)
d = log(1+data$TOP40_EQY_DVD_YLD_12M[ind][1]/100)


atm_strike = which((abs(S/K-1)) == min(abs(S/K-1)))



implied_vol(sig = 0.2, Callval = y[atm_strike], S = S, K = K[atm_strike], t = t, r = 0.07)


ImpVol_90 = nlm(f = implied_vol, p = 0.2, Callval = y[atm_strike],  S= S, K = K[atm_strike], r = r, t = t)
ImpVol_90$estimate  #0.1729683
y_black_90 = Blacks(F0 = S, K = K, t = t, r = r, sd = ImpVol_90$estimate, call = TRUE, DT = 1)

#0.9479766 -3.2480291 -0.5784336 -0.6955646 -3.2188758
x0 = c(0.9479766,-3.2480291,-0.5784336, -0.6955646, -3.2480291)
#x0 = c(log(0.1),log(0.04),log(0.2), 0, log(0.04))
opt_results_nlm_90 = nlm(f = heston_call_pricer, p = x0, 
                  y = y, S =S ,K =K ,t = t,r = r ,d = d, hessian = T)

x1 = (opt_results_nlm_90$estimate)
opt_results_nlm_90$minimum

y_check_90 = call_price_cm_vectorised(model = 'cf_heston', S = S, K = K ,t = t ,r = r ,d = 0, 
                                   rho = tanh(x1[4]), theta =  exp(x1[2]), 
                                   omega =  exp(x1[3]), kappa = exp(x1[1]), V0 = exp(x1[5]))




#sum((y_check_90$call_price - y)^2/length(y))

#cbind(y_check_90$call_price, y)
#plot((y_check_90$call_price - y)/y)

plot(x = K, y = y, type = 'l', ylab = 'Price', xlab = 'Strike',
     xlim = c(46000, 57000), ylim =  c(0, 4200),lwd = 2)
lines(x = K, y = y_check_90$call_price, col = 'blue',
      lwd = 2, type = 'l', lty = 2)
lines(x = K, y = y_black_90, col = 'cyan',lwd = 2, lty = 2)
legend(x = 53500, y = 4000,legend = c('True Price', 'Heston','Black'), 
       col = c('black','blue', 'cyan'),lty = c(1,2,2), bty = 'n', lwd = 2)
grid()


param_name = c("kappa", "theta", "omega", "rho", "V0")
param_90 = data.frame(opt_results_nlm_90$estimate)
row.names(param_90) <- param_name
write.csv(x = param_90, file = 'param_90.csv')
#############################################################################
#Calibrating at 59 days
#############################################################################
#kappa, theta, omega, rho,V0
#other parameters

#ind_rem = seq(from = 2, to  = 193, by = 2)
ind = which(data$i ==2 & data$DAY_TILL_EXPIRY ==59)
##ind = ind[-ind_rem]
#ind = sample(ind, size = length(ind), replace = F)
#?base::sample

y = data$CALL_PRICE[ind]
S = data$AI1_LAST[ind][1]
K = data$Strike[ind]
t = data$DAY_TILL_EXPIRY[ind][1]/365
r = log((1+(data$JIBA3M[ind][1]/100)/4)^4)
d = log(1+data$TOP40_EQY_DVD_YLD_12M[ind][1]/100)


atm_strike = which((abs(S/K-1)) == min(abs(S/K-1)))


implied_vol(sig = 0.2, Callval = y[atm_strike], S = S, K = K[atm_strike], t = t, r = 0.07)


ImpVol_59 = nlm(f = implied_vol, p = 0.2, Callval = y[atm_strike],  S= S, K = K[atm_strike], r = r, t = t)
ImpVol_59$estimate
y_black_59 = Blacks(F0 = S, K = K, t = t, r = r, sd = ImpVol_59$estimate, call = TRUE, DT = 1)

#0.9479766 -3.2480291 -0.5784336 -0.6955646 -3.2188758
x0 = c(0.9479766, -3.2480291, -0.5784336, -0.6955646, -3.2188758)
#x0 = x1
opt_results_nlm_59 = nlm(f = heston_call_pricer, p = x0, 
                       y = y, S =S ,K =K ,t = t,r = r ,d = d, hessian = T)

x1 = (opt_results_nlm_59$estimate)
#solve(opt_results_nlm_59$hessian)
opt_results_nlm_59$minimum

y_check_59 = call_price_cm_vectorised(model = 'cf_heston', S = S, K = K ,t = t ,r = r ,d = 0, 
                                   rho = tanh(x1[4]), theta =  exp(x1[2]), 
                                   omega =  exp(x1[3]), kappa = exp(x1[1]), V0 = exp(x1[5]))




#sum((y_check_59$call_price - y)^2/length(y))

#cbind(y_check_59$call_price, y)
#plot((y_check_59$call_price - y)/y)

plot(x = K, y = y, type = 'l', ylab = 'Price', xlab = 'Strike',
     xlim = c(46000, 57000), ylim =  c(0, 5000),lwd = 2)
lines(x = K, y = y_check_59$call_price, col = 'blue',
      lwd = 2, type = 'l', lty = 2)
lines(x = K, y = y_black_59, col = 'cyan',lwd = 2, lty = 2)
legend(x = 53500, y = 5000,legend = c('True Price', 'Heston','Black'), 
       col = c('black','blue', 'cyan'),lty = c(1,2,2), bty = 'n', lwd = 2)
grid()


param_name = c("kappa", "theta", "omega", "rho", "V0")
param_59 = data.frame(opt_results_nlm_59$estimate)
row.names(param_59) <- param_name
write.csv(x = param_59, file = 'param_59.csv')

#############################################################################
#Calibrating at 30 days
#############################################################################


#kappa, theta, omega, rho,V0

#ind_rem = seq(from = 2, to  = 193, by = 2)
ind = which(data$i == 2 & data$DAY_TILL_EXPIRY ==30)
##ind = ind[-ind_rem]
#ind = sample(ind, size = length(ind), replace = F)
#?base::sample

y = data$CALL_PRICE[ind]
S = data$AI1_LAST[ind][1]
K = data$Strike[ind]
t = data$DAY_TILL_EXPIRY[ind][1]/365
r = log((1+(data$JIBA3M[ind][1]/100)/4)^4)
d = log(1+data$TOP40_EQY_DVD_YLD_12M[ind][1]/100)


atm_strike = which((abs(S/K-1)) == min(abs(S/K-1)))

implied_vol(sig = 0.2, Callval = y[atm_strike], S = S, K = K[atm_strike], t = t, r = 0.07)


ImpVol_30 = nlm(f = implied_vol, p = 0.2, Callval = y[atm_strike],  S= S, K = K[atm_strike], r = r, t = t)
#ImpVol_30$estimate
y_black_30 = Blacks(F0 = S, K = K, t = t, r = r, sd = ImpVol_30$estimate, call = TRUE, DT = 1)

#0.9479766 -3.2480291 -0.5784336 -0.6955646 -3.2188758
x0 = c(0.9479766, -3.2480291, -0.5784336, -0.6955646, -3.2188758)
# c(log(0.1),log(0.04),log(0.2), 0, log(0.04))
opt_results_nlm_30 = nlm(f = heston_call_pricer, p = x0, 
                       y = y, S =S ,K =K ,t = t,r = r ,d = d, hessian = T)

x1 = (opt_results_nlm_30$estimate)
#solve(opt_results_nlm_30$hessian)
#opt_results_nlm_30$minimum

y_check_30 = call_price_cm_vectorised(model = 'cf_heston', S = S, K = K ,t = t ,r = r ,d = 0, 
                                   rho = tanh(x1[4]), theta =  exp(x1[2]), 
                                   omega =  exp(x1[3]), kappa = exp(x1[1]), V0 = exp(x1[5]))




#sum((y_check_30$call_price - y)^2/length(y))

#cbind(y_check_30$call_price, y)
#plot((y_check_30$call_price - y)/y)

plot(x = K, y = y, type = 'l', ylab = 'Price', xlab = 'Strike',
     xlim = c(46000, 57000), ylim =  c(0, 5300),lwd = 2)
lines(x = K, y = y_check_30$call_price, col = 'blue',
      lwd = 2, type = 'l', lty = 2)
lines(x = K, y = y_black_30, col = 'cyan',lwd = 2, lty = 2)
legend(x = 53500, y = 5000,legend = c('True Price', 'Heston','Black'), 
       col = c('black','blue', 'cyan'),lty = c(1,2,2), bty = 'n', lwd = 2)
grid()

param_name = c("kappa", "theta", "omega", "rho", "V0")
param_30 = data.frame(opt_results_nlm_30$estimate)
row.names(param_30) <- param_name
write.csv(x = param_30, file = 'param_30.csv')

#############################################################################
#Simulation
#############################################################################





############################################################
#Simulation
############################################################


param_nlm = (opt_results_nlm4$estimate)

kappa = exp(param_nlm[1])
V0 = exp(param_nlm[2])
theta = exp(param_nlm[3])
omega = exp(param_nlm[4])
rho = tanh(param_nlm[5])

c(kappa, V0, theta, omega, rho)


N <- 500
z1 <-  rnorm(N)
z2 <- rnorm(N)
rho <- 0.1

w2 <- z1
w1 <- rho*z1 +  sqrt(1-rho^2)*z2

V = numeric(N); Ft = numeric(N); lnFt = numeric(N)
V[1]= V0
lnFt[1] = log(Ft0)
Ft[1] = Ft0

for (day in 2:N){
  V[day] = V[day-1] + kappa*(theta -max(0,V[day-1]))*DT + nu*max(0,V[day-1])*w2[day]*sqrt(DT)
  lnFt[day] = lnFt[day-1] - (1/2) * max(0, V[day]) * DT + sqrt(max(0, V[day]))*w1[day]*sqrt(DT)
  Ft[day] = Ft[day-1]exp(lnFt[day])
}




##################
#Plotting Densities
####################



densityrecovery <- function(model, a, N, ...)
{
  cf = get(model)
  b = a/N
  u = ((0:(N-1)) - N/2)*b
  h2 = ((-1)^(0:(N-1))) * cf(u, ...)
  g = fft(h2)
  y = Re((-1)^(0:(N-1))* g ) * (b / (2 * pi))
  return(y)
}


y.plot = densityrecovery(model = 'cf_heston', a = 2*pi*96,N = 2^8, lnS = log(S) ,t = t ,r = r ,d = 0, 
                         rho = tanh(x2[5]), theta =  sig(x2[3]), omega =  sig(x2[4]), kappa = sig(x2[1]), V0 =  sig(x2[2]))

plot(y.plot, type = 'l')


pdf1 = c()  
Nt = 1500
for(i in 1:Nt){
  x = i
  a=-5; b=5;N=Nt 
  k <- seq (1,N -1)
  model <- get( 'cf_heston' )
  Fk <- function (k){
    2/(b-a) * Re ( model (k*pi/(b-a),
                          lnS = log(S) ,t = t ,r = r ,d = 0,
                          rho = tanh(x2[5]), theta =  sig(x2[3]), 
                          omega =  sig(x2[4]), kappa = sig(x2[1]), 
                          V0 =  sig(x2[2])) * exp (-1i*k*a*pi/(b-a)))
                          
  }
  
  Fk.k <- vector ("double", N -1)
  Fk.k <- Fk(k)
  Fk.o <- Fk (0)
  trans <- vector ("double", N)
  trans <- 0.5*Fk.o +sum (Fk.k* cos(k*pi*(x-a)/(b-a)))
  pdf1 = c(pdf1,trans)
}


plot(pdf1)
lines(pdf)









