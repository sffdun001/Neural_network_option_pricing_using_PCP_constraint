#Bloomberg option price data extraction
#Duncan Saffy dissertation

#####################################################
### Collecting Option Specific data
#####################################################

#install.packages("Rblpapi")
library(Rblpapi)
con <- blpConnect() 

setwd("E:/")
#General Data----------------------------------------------------------------------------------------------------------
# Here we collect general data for the first 2 years of Top40, SAVI and AI1 INDEX
# Choose appropriate ticker and field names to pull correct data
ticker_names <- c("AI1 INDEX","AI1 INDEX","AI1 INDEX","AI1 INDEX","AI1 INDEX","AI1 INDEX",
                  "AI1 INDEX","AI1 INDEX", "SAVIT40 INDEX", "JIBA3M INDEX",
                  "TOP40 INDEX", "TOP40 INDEX", "TOP40 INDEX",
                  "TOP40 INDEX", "TOP40 INDEX", "TOP40 INDEX","TOP40 INDEX", "TOP40 INDEX", "TOP40 INDEX")
field_names <- c("PX_LAST", "PX_OPEN", "PX_HIGH", "PX_LOW", "PX_MID", "PX_VOLUME", 
                 "PX_BID","PX_ASK","PX_LAST", "PX_LAST",
                 "BEST_DIV_YLD", "GROSS_AGGTE_DVD_YLD", "EQY_DVD_YLD_12M",
                 "PX_LAST", "PX_OPEN", "PX_HIGH", "PX_LOW", "PX_MID", "PX_VOLUME") 



data_store = data.frame()

data_store = bdh(ticker_names[1],field_names[1], Sys.Date()-365*3)

#Retrieve general data from bloomberg on AI1Index
for (i in 2:length(ticker_names)){
  columns <- bdh(ticker_names[i],field_names[i], Sys.Date()-365*3)
  data_store = merge(data_store, columns, by = "date", all.x =TRUE)
}
gen_data_names = c("no.","Date","AI1_LAST","AI1_OPEN","AI1_HIGH",
                   "AI1_LOW","AI1_MID","AI1_VOLUME",
                   "AI1_BID","AI1_ASK", "SAVIT40", 
                   "TOP40", "JIBA3M","TOP40_BEST_DIV_YIELD",
                   "TOP40_GROSS_AGGTE_DVD_YLD", "EQY_DVD_YLD_12M",
                   "TOP40_LAST","TOP40_OPEN","TOP40_HIGH",
                   "TOP40_LOW","TOP40_MID","TOP40_VOLUME"
)
names(data_store) = gen_data_names

write.csv(x = data_store, file = "BLB_GEN_DATA.csv") #Read in file of General Data



#####################################################
### Collecting Option Specific data
#####################################################
#Collecting, put/call prices, strikes and Expiries

#Read in with names of various option contracts names, which can be retrieved from bloomberg
M9 <- read.csv("Options_M9.csv")
U9 <- read.csv("Options_U9.csv")



opt_list = list(M9, U9)
#Retrieve data specific to each option contract
for (i in 1:2)
{
  fut_ind = substr(opt_list[[i]][1,2],start = 1, stop = 4)#future ind, i.e. AIM9
  fut_px =  bdh(paste(fut_ind, "INDEX", sep = " "),"PX_LAST", Sys.Date()-365.25*3)
  fut_vol =  bdh(paste(fut_ind, "INDEX", sep = " "),"PX_VOLUME", Sys.Date()-365.25*3)
  fut_high =  bdh(paste(fut_ind, "INDEX", sep = " "),"PX_HIGH", Sys.Date()-365.25*3)
  fut_low =  bdh(paste(fut_ind, "INDEX", sep = " "),"PX_LOW", Sys.Date()-365.25*3)
  fut_mid =  bdh(paste(fut_ind, "INDEX", sep = " "),"PX_MID", Sys.Date()-365.25*3)
  fut_open =  bdh(paste(fut_ind, "INDEX", sep = " "),"PX_OPEN", Sys.Date()-365.25*3)
  #putting all items in a list
  fut_data = list(fut_px, fut_vol, fut_high, fut_low, fut_mid, fut_open)
  
  for (j in 1:nrow(opt_list[[i]])){
    strk = substr(opt_list[[i]][j,2],start = 6, stop = 17) #strike_INDEX i.e 55000 INDEX
    #Reading in data from bloomberg
    opt_data_put_px = bdh(paste(fut_ind, "P", strk,sep = ""),"PX_LAST", Sys.Date()-365.25*3)
    opt_data_put_vol = bdh(paste(fut_ind, "P", strk,sep = ""),"PX_VOLUME", Sys.Date()-365.25*3)
    opt_data_call_px = bdh(paste(fut_ind, "C", strk,sep = ""),"PX_LAST", Sys.Date()-365.25*3)
    opt_data_call_vol = bdh(paste(fut_ind, "C", strk,sep = ""),"PX_VOLUME", Sys.Date()-365.25*3)
    opt_data_exp = bdp(paste(fut_ind, "C", strk,sep = ""),"LAST_TRADEABLE_DT")
    opt_data_strike = as.numeric(substr(opt_list[[i]][j,2],start = 7, stop = 11))
    
    #Merging all data into one file
    opt_data = data.frame(NA, nrow = max(nrow(opt_data_put_px), nrow(opt_data_call_px)),ncol = 8)
   
    #merging option price data
    opt_data = opt_data_call_px
    opt_data = merge(x = opt_data, y = opt_data_call_vol, by = "date", all.x =TRUE, all.y = TRUE)
    opt_data = merge(x = opt_data, y = opt_data_put_px, by = "date", all.x =TRUE, all.y = TRUE)
    opt_data = merge(x = opt_data, y = opt_data_put_vol, by = "date", all.x =TRUE, all.y = TRUE)
    
    #merging futures data
    for (lst in fut_data){
      opt_data = merge(x = opt_data, y = lst, by = "date", all.x =TRUE, all.y = TRUE)
    }
    #appending strike and expriry
    Strike= matrix(rep(opt_data_strike, nrow(opt_data)), nrow(opt_data), 1)
    Expiry= matrix(rep(as.character(opt_data_exp$LAST_TRADEABLE_DT), nrow(opt_data)), nrow(opt_data), 1)
    opt_data = cbind(opt_data, Strike, Expiry)
    
    write.csv(x = opt_data, file = paste(fut_ind, strk,".csv",sep =""))
    #print(j)
  }
}

