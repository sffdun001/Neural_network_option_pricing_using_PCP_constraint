#Duncan Saffy Dissertation
#Market Data Explotation

rm(list = ls(all = TRUE))
setwd("~/UCT Academic/Masters/Dissertation/GitHub/Market Data")
library(dplyr)
library(lubridate)
library(tidyverse)

AIM9 = read.csv('AIM9_processed_data_2019-10-07.csv', header = T)
AIM8 = read.csv('AIM8_processed_data_2019-07-22.csv', header = T)
AIH8 = read.csv('AIH8_processed_data_2019-07-22.csv', header = T)
AIH9 = read.csv('AIH9_processed_data_2019-07-22.csv', header = T)
AIU8 = read.csv('AIU8_processed_data_2019-07-22.csv', header = T)
AIU9 = read.csv('AIU9_processed_data_2019-10-07.csv', header = T)
AIZ8 = read.csv('AIZ8_processed_data_2019-07-22.csv', header = T)

Data = list(AIH8, AIM8, AIU8, AIZ8, AIH9, AIM9, AIU9)
data_names = c('AIH8', 'AIM8', 'AIU8', 'AIZ8', 'AIH9', 'AIM9', 'AIU9')
names(Data) <- data_names


GEN = read.csv('BLB_GEN_DATA_Oct.csv', header = T)
nrow(AIM9)
nrow(AIM8)
#View(AIM8)

ncol(AIH8)
AIH8_Ex = AIH8[,-c(1,19)]
AIH8_Ex[1,] = names(AIH8_Ex)
names(AIH8_Ex) <- NULL
write.csv(rbind(as.matrix(AIH8_Ex[1:5,1:6]),
                as.matrix(AIH8_Ex[1:5,7:12]),
                as.matrix(AIH8_Ex[1:5,13:18]),
                as.matrix(AIH8_Ex[1:5,19:24])),
          file = 'data_example.csv')

head(AIH8_Ex)
######################################################################################################
#Plotting continuation contract over period from which data was collected
######################################################################################################


expiration_dates <- c('2017-12-21','2018-03-15','2018-06-21','2018-09-20',
                      '2018-12-20','2019-03-21','2019-06-20','2019-09-19')

gen_data = GEN%>%mutate(date = as_date(date))%>%
                        filter(date >= as.Date("2017-12-21") & date <= as.Date("2019-09-19"))

min(gen_data$PX_LAST.x)
max(gen_data$PX_LAST.x)

######################################################################################################
#Check if data needs to be cleaned
######################################################################################################







####################

for( i in 1:length(Data)){
  print(c(min(Data[[i]]$Strike, na.rm = T), 
          max(Data[[i]]$Strike, na.rm = T)))
}


min(Data[[i]]$Strike, na.rm = T)


#Plotting the Continuation contract over period collected 
plot(y = gen_data$PX_LAST.x, x = gen_data$date,
     type = 'l', xaxt = 'n', xlab = 'Date',
     ylab = 'Price'
     )
axis.Date(1, at = seq(from = as.Date('2018-03-02'), to = as.Date('2019-09-02'), 
                      length.out=7), format= "%Y-%m")

#grid(nx = 7, ny = 12, col = "lightgray",
#     lwd = par("lwd"), equilogs = TRUE)

for(d in expiration_dates[-1]){
  abline(v =  as.Date(d), col = 'darkgrey', lty = 2)
  
}
legend('bottomleft', legend = c('AI1 Index', 'Expirations'), 
        lty = c(1,2), bty = 'n', col = c('black', 'grey'), lwd = 2)


######################################################################################################
#Plotting number of contracts available as a function of days to expiry
######################################################################################################

#This removes observations where both put and call prices are NA
AIH8.Clean <- AIH8%>%filter(!is.na(CALL_PRICE) | !is.na(PUT_PRICE))

#select first 8 columns of data and by each strike find the first day of trading and the last as per dataset


AIH8.Clean.Subset = AIH8.Clean%>%select(1:8)%>%
               group_by(Strike)%>%
               mutate(start_time = max(DAY_TILL_EXPIRY), end_time = min(DAY_TILL_EXPIRY))%>%
               filter(row_number() == 1)%>%ungroup()

strike_counts = as.matrix(table(AIH8.Clean.Subset$start_time))

no_strikes = cumsum(rev(strike_counts))
starttimes = unique(AIH8.Clean.Subset$start_time)
creation_days = starttimes[rev(order(starttimes))]

contracts = list()
for (i in 1:length(Data)){
  df = Data[[i]]
  df.Clean <- df%>%filter(!is.na(CALL_PRICE) | !is.na(PUT_PRICE))
  df.Clean.Subset = df.Clean%>%select(1:8)%>%
    group_by(Strike)%>%
    mutate(start_time = max(DAY_TILL_EXPIRY), end_time = min(DAY_TILL_EXPIRY))%>%
    filter(row_number() == 1)%>%ungroup()
  #print(nrow(df.Clean.Subset))
  df.strike_counts = as.matrix(table(df.Clean.Subset$start_time))
  
  df.no_strikes = cumsum(rev(df.strike_counts))
  df.starttimes = unique(df.Clean.Subset$start_time)
  df.creation_days = df.starttimes[rev(order(df.starttimes))]
  print(c(min(df.Clean$Strike, na.rm = T), 
          max(df.Clean$Strike, na.rm = T)))
  contracts[[i]] = cbind(df.no_strikes, df.creation_days)
}


length(no_strikes)

plot(0,xlim  = c(95, 0), ylim = c(0,200), xlab = 'Days till expiration',
     ylab = 'Number of contracts', type = 'n')
for (i in 1:length(contracts)){
  lines(x = c(contracts[[i]][,2],0), y = c(contracts[[i]][,1], max(contracts[[i]][,1]))+ rnorm(1,0,2),
        type = 'S', col = colors()[50+i*5], lwd = 2)
}

legend('bottomright', legend = data_names, fill = c(colors()[seq(5,35, by = 5)+50]),
       bty = 'n')

AIH8.Clean.Subset%>%group_by(start_time)

colors()[seq(10,70, by = 10)]



######################################################################################################
#Comparing affect of removing NA values
######################################################################################################
Option_suffix = c('H8', "M8", "U8", "Z8", 'H9', "M9", "U9")


data_main <- vector(mode = "list", length = 7)

for(i in 1:7){ #not including last 2 for test set
  if (i ==6|i ==7){
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
  
  quarter_data.clean$InterestRate = log((1+quarter_data.clean$InterestRate/400)^4)
  
  
  data_main[[i]]= quarter_data.clean
}

Rows_removing = vector(mode = "list", length =  7)

for (i in 1:7){
  data = data_main[[i]]
  Rows_removing[[i]] = matrix(NA, nrow = 2000, ncol = 4)
  l = 1
  for (j in 1:nrow(data)){
    if(data$CallPrice[j] == 0){
      #print(paste(i, ":",j, sep = ""))
      if(sum(data$CallPrice[which((data$Strike == data$Strike[j]) & (data$TimeTillExpiry >=data$TimeTillExpiry[j]))]) == 0 ){
        info1 = c(j, data$StockPrice[j], data$Strike[j], data$TimeTillExpiry[j])
        Rows_removing[[i]][l,] =  info1
        print(paste(i, ":",j, sep = ""))
        l = l+1
      }
    }
    
  }
}


###############################################################################
#head(quarter_data.clean)
#length(na.omit(Rows_removing[[7]][,1]))

#We remove this data point as neven had call price
data_main[[7]] = data_main[[7]][-na.omit(Rows_removing[[7]][,1]),]
data_main[[3]] = data_main[[3]][-na.omit(Rows_removing[[3]][,1]),]





#######################
d2 = vector(mode = 'list', length = 7)
for (quarter in c(1:7)){
  #quarter = 1
  data <-  data_main[[quarter]]%>%filter(!is.na(CallPrice) | !is.na(PutPrice))
  Strikes = sort(unique(data$Strike))
  d2[[quarter]] = matrix(NA, nrow = length(Strikes), ncol = 3)
  for (Strike in 1:length(Strikes)){
    #Strike = 1
    d1 = data[which(data$Strike == Strikes[Strike]),]
    a1 = d1[which(d1$TimeTillExpiry==max(d1$TimeTillExpiry)),]
    d2[[quarter]][Strike,1] = a1$TimeTillExpiry
    d2[[quarter]][Strike,2] = a1$Strike
    d2[[quarter]][Strike,3] = a1$StockPrice/a1$Strike
  }
}

head(d2[[1]])
Moneyness_table = matrix(NA, nrow = 7, ncol = 4)
for (quarter in c(1:7)){
  Moneyness_table[quarter,] =  c( nrow(d2[[quarter]]), 
                                  length(which(d2[[quarter]][,3]>1.03)),
                                  length(which((d2[[quarter]][,3]<1.03)&(d2[[quarter]][,3]>0.97))), 
                                  length(which((d2[[quarter]][,3]<0.97))) 
  )
}
#head(data_main1[[1]])
colSums(Moneyness_table)
write.csv(as.data.frame(Moneyness_table), file = "MoneynessTable.csv")









