#####################################################
#Hutchinson Data Simulation
#####################################################


rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/Technical Paper/Hutchingson")

library(dplyr)


#################################
#Model
#################################

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

HUTCHINSON_BS_SIMULATION <- function(S0, mu, r, vol)
{ #require dplyr
  #S0 = 50; r =0.1; vol = 0.2
  
  df = data.frame(Day = 0, TimeTillExpiry = 0, StockPrice =0, Strike = 0, CallPrice = 0, PutPrice = 0)
  St = numeric(505)
  St[1] = S0
  
  expiry = seq(from = 21, to = 504, by =21)
  DAYS_TILL_EXPIRY_1 = c(21,42,63,126)
  DAYS_TILL_EXPIRY_2 = c(21,42,105,168)
  DAYS_TILL_EXPIRY_3= c(21,42,84,147)
  
  POTENTIAL_STRIKES = seq(25,200, by = 5)
  TIMES = list(DAYS_TILL_EXPIRY_1, DAYS_TILL_EXPIRY_1, DAYS_TILL_EXPIRY_1)
  strikes = c(45,50,55); strike_upper = 45; strike_lower = 55
  times = TIMES
  
  for(day in 1:504){
    
    St[day+1] = St[day]*exp(rnorm(1, mean = mu/252, sd = vol/sqrt(252))) 
    
    #Removing strikes that have expired
    strikes_to_remove = c()
    for(i in 1:length(times)){
      if (any(times[[i]] == 1) & (length(times[[i]])==1)){
        strikes_to_remove = c(strikes_to_remove, i)
      } 
    }
    
    if(length(strikes_to_remove)>0){
      strikes = strikes[-strikes_to_remove]
      times = times[-c(strikes_to_remove)]
    }
    
    
    #Decreaing time by 1 and removing expiry dates that come
    for(i in 1:length(times)){
      times[[i]] =times[[i]] - 1
      if (any(times[[i]] == 0)){
        times[[i]] = times[[i]][-which(times[[i]]==0)]
      }
    }
    
    
    
    if(day%%21 == 0){ #Create new strikes each month
      if(day%%63 == 0){
        days_till_expiry = DAYS_TILL_EXPIRY_1
      }else if((day%%63)%%42 == 0){
        days_till_expiry = DAYS_TILL_EXPIRY_3
      }else{
        days_till_expiry = DAYS_TILL_EXPIRY_2
      }
      
      
      #Setting up new strikes for current exiprtation
      strike_ind = which(abs(POTENTIAL_STRIKES - St[day+1]) < 6)
      new_strikes = POTENTIAL_STRIKES[strike_ind]
      new_times = replicate(length(new_strikes), days_till_expiry, simplify=FALSE)
      #Identify new upper and lower strike prices
      strike_upper = max(new_strikes)
      strike_lower = min(new_strikes)
      
      #Append strikes currently not expired and not being rewritten
      
      if(sum(!strikes%in%new_strikes)>0){
        strike_keeping_index = which(!strikes%in%new_strikes)
        new_strikes = c(new_strikes, strikes[strike_keeping_index])
        for (j in strike_keeping_index){
          new_times[[length(new_times)+1]] =  times[[j]]
        }
      }
      
      #Update old strikes to new strikes
      strikes = new_strikes
      times = new_times
      
    }
    
    #Checking if vector of if need to add new strike due to move out of price bracket
    #Finding what dates to set the new strike for if price exceeds brackets
   
    #if (day%%21 > 16){
    #  if(day%%63 > 42){
    #    days_till_expiry = DAYS_TILL_EXPIRY_1 + (21 - day%%21)
    #  }else if((day%%63)%%63 > 21){
    #    days_till_expiry = DAYS_TILL_EXPIRY_3 + (21 - day%%21)
    #  }else{
    #    days_till_expiry = DAYS_TILL_EXPIRY_2 + (21 - day%%21)
    #  }
    #} else{
    #  days_till_expiry = times[[1]]
    #}
    days_till_expiry = times[[1]] #delete
    
    if (day%%21 < 16)
    {
    if(St[day+1] > strike_upper){
      if((strike_upper + 5)%in%strikes){
        #Update available strike dates
        index = which(strikes == (strike_upper + 5))
        times[[index]] = unique(c(times[[index]],days_till_expiry))
      }else {
        #Create new strike and give times
        strikes = c(strikes, (strike_upper+5))
        times[[length(times)+1]] = days_till_expiry
      }
      strike_upper = strike_upper +5
      
    } else if(St[day+1] < strike_lower){
      if((strike_lower - 5)%in%strikes){
        #Update available strike dates
        index = which(strikes == (strike_lower - 5))
        times[[index]] = unique(c(times[[index]],days_till_expiry))
        
      }else {
        #Create new strike and give times
        strikes = c(strikes, (strike_lower-5))
        times[[length(times)+1]] = days_till_expiry
      }
      strike_lower = strike_lower -5
    }
    }
    
    
    #calculating option prices and adding to store data frame
    for (k in 1:length(strikes)){
      call_price = BlackScholes(S = St[day+1], K = strikes[k], t = times[[k]],
                                r = r, sd = vol, call = TRUE)
      put_price = BlackScholes(S = St[day+1], K = strikes[k], t = times[[k]],
                               r = r, sd = vol, call = FALSE)
      
      to_append = cbind(Day = day, TimeTillExpiry = times[[k]], StockPrice = St[day+1], Strike = strikes[k], 
                        CallPrice = call_price, PutPrice = put_price)
      
      df = rbind(df, to_append)
      #Create vector to append to data frame
    }
    
  }
  
  
  df = df[-1,] #removing the initial row of zeros created
  OutOfScope <- which(((504-df$Day) - df$TimeTillExpiry)<0)#Removing observations that expire after final time.
  df = df[-OutOfScope,]
  df1 <- df%>%group_by(Strike)%>%mutate(TerminalDate = TimeTillExpiry + Day)%>%
    arrange(Strike, TerminalDate, desc(TimeTillExpiry))%>%ungroup()
  no_price_paths = nrow(count(df1, Strike, TerminalDate))
  #Sort the data so that repective strikes are together in order or date
  return(list(df = as.data.frame(df1) , no_price_paths = no_price_paths))
}

HUTCHINSON_BS_SIMULATION_TEST <- function(S0, mu, r, vol)
{ #require dplyr
  #S0 = 50; r =0.1; vol = 0.2
  
  df = data.frame(Day = 0, TimeTillExpiry = 0, StockPrice =0, Strike = 0, CallPrice = 0, PutPrice = 0)
  St = numeric(505)
  St[1] = S0
  
  expiry = seq(from = 21, to = 504, by =21)
  DAYS_TILL_EXPIRY_1 = c(21,42,63,126)
  DAYS_TILL_EXPIRY_2 = c(21,42,105,168)
  DAYS_TILL_EXPIRY_3= c(21,42,84,147)
  
  POTENTIAL_STRIKES = seq(25,200, by = 5)
  TIMES = list(DAYS_TILL_EXPIRY_1, DAYS_TILL_EXPIRY_1, DAYS_TILL_EXPIRY_1,DAYS_TILL_EXPIRY_1,DAYS_TILL_EXPIRY_1)
  strikes = c(40,45,50,55, 60); strike_upper = 40; strike_lower = 60
  times = TIMES
  
  for(day in 1:504){
    
    St[day+1] = St[day]*exp(rnorm(1, mean = mu/252, sd = vol/sqrt(252))) 
    
    #Removing strikes that have expired
    strikes_to_remove = c()
    for(i in 1:length(times)){
      if (any(times[[i]] == 1) & (length(times[[i]])==1)){
        strikes_to_remove = c(strikes_to_remove, i)
      } 
    }
    
    if(length(strikes_to_remove)>0){
      strikes = strikes[-strikes_to_remove]
      times = times[-c(strikes_to_remove)]
    }
    
    
    #Decreaing time by 1 and removing expiry dates that come
    for(i in 1:length(times)){
      times[[i]] =times[[i]] - 1
      if (any(times[[i]] == 0)){
        times[[i]] = times[[i]][-which(times[[i]]==0)]
      }
    }
    
    
    
    if(day%%21 == 0){ #Create new strikes each month
      if(day%%63 == 0){
        days_till_expiry = DAYS_TILL_EXPIRY_1
      }else if((day%%63)%%42 == 0){
        days_till_expiry = DAYS_TILL_EXPIRY_3
      }else{
        days_till_expiry = DAYS_TILL_EXPIRY_2
      }
      
      
      #Setting up new strikes for current exiprtation
      strike_ind = which(abs(POTENTIAL_STRIKES - St[day+1]) < 6)
      new_strikes = POTENTIAL_STRIKES[strike_ind]
      new_times = replicate(length(new_strikes), days_till_expiry, simplify=FALSE)
      #Identify new upper and lower strike prices
      strike_upper = max(new_strikes)
      strike_lower = min(new_strikes)
      
      #Append strikes currently not expired and not being rewritten
      
      if(sum(!strikes%in%new_strikes)>0){
        strike_keeping_index = which(!strikes%in%new_strikes)
        new_strikes = c(new_strikes, strikes[strike_keeping_index])
        for (j in strike_keeping_index){
          new_times[[length(new_times)+1]] =  times[[j]]
        }
      }
      
      #Update old strikes to new strikes
      strikes = new_strikes
      times = new_times
      
    }
    
    #Checking if vector of if need to add new strike due to move out of price bracket
    #Finding what dates to set the new strike for if price exceeds brackets
    
    #if (day%%21 > 16){
    #  if(day%%63 > 42){
    #    days_till_expiry = DAYS_TILL_EXPIRY_1 + (21 - day%%21)
    #  }else if((day%%63)%%63 > 21){
    #    days_till_expiry = DAYS_TILL_EXPIRY_3 + (21 - day%%21)
    #  }else{
    #    days_till_expiry = DAYS_TILL_EXPIRY_2 + (21 - day%%21)
    #  }
    #} else{
    #  days_till_expiry = times[[1]]
    #}
    days_till_expiry = times[[1]] #delete
    
    if (day%%21 < 16)
    {
      if(St[day+1] > strike_upper){
        if((strike_upper + 5)%in%strikes){
          #Update available strike dates
          index = which(strikes == (strike_upper + 5))
          times[[index]] = unique(c(times[[index]],days_till_expiry))
        }else {
          #Create new strike and give times
          strikes = c(strikes, (strike_upper+5))
          times[[length(times)+1]] = days_till_expiry
        }
        strike_upper = strike_upper +5
        
      } else if(St[day+1] < strike_lower){
        if((strike_lower - 5)%in%strikes){
          #Update available strike dates
          index = which(strikes == (strike_lower - 5))
          times[[index]] = unique(c(times[[index]],days_till_expiry))
          
        }else {
          #Create new strike and give times
          strikes = c(strikes, (strike_lower-5))
          times[[length(times)+1]] = days_till_expiry
        }
        strike_lower = strike_lower -5
      }
    }
    
    
    #calculating option prices and adding to store data frame
    for (k in 1:length(strikes)){
      call_price = BlackScholes(S = St[day+1], K = strikes[k], t = times[[k]],
                                r = r, sd = vol, call = TRUE)
      put_price = BlackScholes(S = St[day+1], K = strikes[k], t = times[[k]],
                               r = r, sd = vol, call = FALSE)
      
      to_append = cbind(Day = day, TimeTillExpiry = times[[k]], StockPrice = St[day+1], Strike = strikes[k], 
                        CallPrice = call_price, PutPrice = put_price)
      
      df = rbind(df, to_append)
      #Create vector to append to data frame
    }
    
  }
  
  
  df = df[-1,] #removing the initial row of zeros created
  OutOfScope <- which(((504-df$Day) - df$TimeTillExpiry)<0)#Removing observations that expire after final time.
  df = df[-OutOfScope,]
  df1 <- df%>%group_by(Strike)%>%mutate(TerminalDate = TimeTillExpiry + Day)%>%
    arrange(Strike, TerminalDate, desc(TimeTillExpiry))%>%ungroup()
  no_price_paths = nrow(count(df1, Strike, TerminalDate))
  #Sort the data so that repective strikes are together in order or date
  return(list(df = as.data.frame(df1) , no_price_paths = no_price_paths))
}
#################################
#Data Simualtion
#################################

TRAIN_DATA <- vector(mode = "list", length = 10)
TEST_DATA <- vector(mode = "list", length = 500)

set.seed(1995)

pp = numeric(10)
for(i in 1:10){
  hutch = HUTCHINSON_BS_SIMULATION(S0 = 50, r = 0.04681, mu = 0.1 , vol = 0.2)
  TRAIN_DATA[[i]] = hutch$df
  pp[i] = hutch$no_price_paths
}


pp1 = numeric(500)
for(i in 1:500){
  hutch = HUTCHINSON_BS_SIMULATION_TEST(S0 = 50, r = 0.04681, mu = 0.1 , vol = 0.2)
  TEST_DATA[[i]] = hutch$df
  pp1[i] = hutch$no_price_paths
  print(i)
}

min(pp);mean(pp);max(pp)

hist(pp1)

save(TRAIN_DATA, file = paste('Hutch_train_001_',Sys.Date(), ".RData", sep = ''))
save(TEST_DATA, file = paste('Hutch_test_001_',Sys.Date(), ".RData", sep = ''))


#################################
#Exploratory Statistics
#################################

prac_mean = numeric(500)
for ( i in 1:500){
  
  prac_data =  TEST_DATA[[i]]
  prac_data$StockPrice = as.character(prac_data$StockPrice)
  price_path <- prac_data%>%group_by(StockPrice)%>%
    filter(row_number() ==1)%>%ungroup%>%
    arrange(Day)
  #View(price_path)
  
  #plot(price_path$StockPrice, type = 'l')
  
  prac_returns = log(as.numeric(price_path$StockPrice[-1])/as.numeric(price_path$StockPrice[-503]))
  #hist(prac_returns)
  
  #length(price_path$StockPrice)
  #View(price_path)
  
  
  prac_mean[i] = mean(prac_returns)*252
}

hist(prac_mean)
abline(v = mean(prac_mean))


