#This code gives us plots and figures for the ALSI data extracted from bloomberg

rm(list = ls(all.names = TRUE)) 
setwd("~/UCT Academic/Masters/Dissertation/Technical Paper/Data")
data = read.csv("AIH8_processed_data_2019-07-22.csv", header = T)

Black<- function(Ft, K, t, sd, r, call = TRUE) #So this apparently works
{
  d1 <- 1/(sd*sqrt(t/365))*(log(Ft/K) +((sd^2)/2)*(t/365))
  d2 <- d1 - sd*sqrt(t/365)
  C <- (Ft*pnorm(d1) - K*pnorm(d2))*exp(-r*t / 365)
  if(call){
    return(C) #returning the call price
  }
  else {
    return(C + K*exp(-r*t/365) - Ft*exp(-r*t/365)) #returning the put price
  }
}

# Plotting data
#Subsetting data and setting NA values to 0
subset_index = which(AIH8$Strike == 52000)
data_subset = AIH8[subset_index,]
data_subset$Date = as.Date(data_subset$Date)
na_index = which(is.na(data_subset$CALL_PRICE))
data_subset$CALL_PRICE[na_index] = 0
blacks_price = Black(Ft = data_subset$AI1_LAST, K= data_subset$Strike, 
                     t = data_subset$DAY_TILL_EXPIRY, r = log((1+data_subset$JIBA3M/400)^4), 
                     sd = data_subset$SAVIT40/100)


par(mar = c(5,5,2,5))
plot(x = data_subset$Date, y = data_subset$CALL_PRICE, main = "",
     type = 'l',yaxs="i", ylim = c(0, 3500), col = 'blue',
     xlab = 'Date', ylab = 'Call Price') 
lines(x =  data_subset$Date, y = blacks_price, col = 'lightblue3')
grid(nx = 8, ny = 7, col = "lightgray",
     lwd = par("lwd"), equilogs = TRUE)
par(new = T)
plot(data_subset$Date, y = data_subset$AI1_LAST, pch=16, axes=F, xlab=NA,
     ylab=NA, type = 'l', col = 'black')
axis(side = 4)
mtext(side = 4, line = 3, 'Futures Price')
abline(h = data_subset$Strike[1], col = 'darkgrey', lty = 2)
legend('topright' ,legend = c('Market', 'Black\'s ', 'Future', 'Strike'),
       fill = c('blue', 'lightblue3', 'black','darkgrey'), bty = 'n', lty = c(1,1,1,2), title = 'Price')

#Plotting % difference between blacks price and observed price

perc_diff = log(data_subset$CALL_PRICE/blacks_price)
plot(x = data_subset$Date, y = perc_diff*100, type = 'l',
     xlab = 'Date', ylab = '% Difference', col = 'blue', ylim = c(-45,45))
abline(h = 0)


#Plotting implied volatility
Black_optim<- function(sd, Ct, Ft, K, t, r, call = TRUE)
{
   B_Ct = Black(Ft =Ft, K = K, 
                t = t, r = r/100, 
                sd =sd, call = call)
   return(abs(Ct - B_Ct))
}
vol_index = which(AIH8$DAY_TILL_EXPIRY == 83)
#View(data[,c(1:8)])
length(vol_index)

vol_data = AIH8[vol_index,]
order(vol_data$Strike)
vols = c()

for(i in order(vol_data$Strike)){
  res <- optim(par = 0.16, fn = Black_optim, Ct = vol_data$CALL_PRICE[i], r = log((1+vol_data$JIBA3M[i]/4)^4),
             Ft = vol_data$AI1_LAST[i], t = vol_data$DAY_TILL_EXPIRY[i],
             K = vol_data$Strike[i], call = T,
             method = 'Brent', lower = 0, upper = 1)
  vols = c(vols, res$par)
}
#cbind(vol_data$CALL_PRICE, vol_data$PUT_PRICE)
ordered_strikes = vol_data$Strike[order(vol_data$Strike)]
plot(y = vols, x = ordered_strikes, 
     xlab = 'Strike', ylab = 'Implied Volatility', 
     #main = 'Implied Volatility at 83 days tills Expiration',
     ylim = c(0,0.25), type = 'l', lwd = 3)
  grid()
  #grid(nx = 8, ny = 7, col = "lightgray",
   #  lwd = par("lwd"), equilogs = TRUE)
abline(v =  vol_data$AI1_LAST[i], col = 'blue', lty = 2)
legend('bottomleft', legend = c( 'Implied Volatility', 'Future Price'),
       fill =  c('black', 'blue'), lty = c(1,2), bty = 'n')




##Ploting put call parity c+Ke^-rt = p + Fe^-rt
vol_data1 = vol_data[order(vol_data$Strike),]

call_side = vol_data1$CALL_PRICE + vol_data1$Strike*exp(-vol_data1$DAY_TILL_EXPIRY*
                                                          (log((1+vol_data1$JIBA3M/400)^4) /365))
put_side = vol_data1$PUT_PRICE + vol_data1$AI1_LAST*exp(-vol_data1$DAY_TILL_EXPIRY*
                                                          (log((1+vol_data1$JIBA3M/400)^4) /365))

cbind(call_side, put_side)

plot(y = call_side, x = vol_data1$Strike,type = 'l',
     lwd =2, col = 'blue', ylim = c(52000, 56000),
     ylab = 'Portfolio value', xlab = 'Strike', 
     xaxs ='i', yaxs = 'i')
lines(y = put_side, x = vol_data1$Strike,type = 'l', 
      lwd =2, col = 'grey')
grid(nx = 5, ny = 8, col = "lightgray",
     lwd = par("lwd"), equilogs = TRUE)
abline(v = vol_data1$AI1_LAST[1], lty= 2)
legend('topleft', legend = c(expression(paste('c + K',e^-rT, sep = '')),
                             expression(paste('p + F',e^-rT, sep = ''))),
      col = c('blue', 'grey'), bty = 'n', lty = c(1,1))

pc_perc_diff = log(call_side/put_side)
plot(x = vol_data1$Strike, y = pc_perc_diff*100, col = 'royalblue',
     xlab = 'Strike', ylab = '% Difference', type = 'l', lwd = 2,
     ylim = c(-1, 1), xaxs ='i', yaxs = 'i')
grid(nx = 5, ny = 8, col = "lightgray",
     lwd = par("lwd"), equilogs = TRUE)
abline(v = vol_data1$AI1_LAST[1], lty = 2)
abline(h = 0)
legend('topleft', legend = c('% Error', 'Future value'),
       col = c('royalblue', 'black'), lty = c(1,2), bty = 'n')











