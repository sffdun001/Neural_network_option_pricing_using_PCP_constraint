#Duncan Saffy
#Data compacting 

###############################################
#Desciption
###############################################
#Having collected option price data this code compacts data of 2 sources into 90 day sections
#Note that you will need to change the file names and expiration dates for this code to run appropriately



setwd("F:/")# Set working Directory 

#Read in with names of various option contracts names, which can be retrieved from bloomberg
M9 <- read.csv("Options_M9.csv")
U9 <- read.csv("Options_U9.csv")

gen_data = read.csv("BLB_GEN_DATA.csv", header = T) 


#########################################
#Combining option price and general data
########################################


option_prices_names <- c('Date', "CALL_PRICE", "PUT_PRICE")
expiration_dates = c("2019-03-21", "2019-06-20", "2019-09-19") #3rd thurday of mar/jun/sep/dec
opt_list = list(M9, U9)#Just removing as already have other data
for (i in 1:2)
{ 
  #i =1
  combined_data  = c()
  ind1 = which(as.Date(gen_data$Date)>as.Date(expiration_dates[i])& (as.Date(gen_data$Date)<=as.Date(expiration_dates[i+1])))
  Index_data = gen_data[ind1,]
  Index_data$Date = as.Date(Index_data$Date)
  #i=1;j=1
  for (j in 1:nrow(opt_list[[i]])){
    #reading in option file
    #j =1
    option_file_name =  paste(substr(opt_list[[i]][j,2], start = 1, stop = 4), substr(opt_list[[i]][j,2], start = 6, stop = 18),".csv", sep = "")
    option_data = read.csv(option_file_name, header = T)
    option_prices <- data.frame("Date" = as.Date(option_data$date), "CALL_PRICE" = option_data$PX_LAST.x, 
                                "PUT_PRICE" = option_data$PX_LAST.y)

    
    ind2 = which(option_prices$Date>as.Date(expiration_dates[i])& (option_prices$Date<=as.Date(expiration_dates[i+1])))
    option_prices = option_prices[ind2,]
    option_prices$Strike = rep(substr(opt_list[[i]][j,2], start = 7, stop = 11), nrow(option_prices))
    option_prices$Expiry = expiration_dates[i+1]
    time_till_expiry = as.numeric(difftime(as.Date(expiration_dates[i+1]), option_prices$Date), units = "days")
    option_prices$DAY_TILL_EXPIRY = time_till_expiry
    combined_genData_optionPrices  = merge(option_prices, Index_data[,-1], by = "Date", all.x = TRUE, all.y = TRUE)
    combined_data = rbind(combined_data, combined_genData_optionPrices)
  }
  file_name = paste(substr(opt_list[[i]][j,2], start = 1, stop = 4),"_processed_data_",format(Sys.time(), "%Y-%m-%d"), ".csv",sep ="")
  write.csv(x = combined_data, file = file_name)
  print(i)
}


