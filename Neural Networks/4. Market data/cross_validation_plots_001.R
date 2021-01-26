#Duncan Saffy: Dissertation Plots

#CV Error and training errors
setwd("~/UCT Academic/Masters/Dissertation/GitHub/Market Data")


#Sets 1-8 are hidden = 4;Sets 9-16 are hidden = 8;Sets 17-24 are hidden = 4,4;Sets 25-32 are hidden = 8,8
#the first 4 of each are learning rate = 0.1 and second 4 are learning rate = 0.01
#Of each of these subsets the lambdas  = 0.001,0.0001,0.00001,0

cv_files = c("mlp_market_cv_001_2020-06-21.RData",
             "mlpSC_market_cv_001_2020-08-14.RData",
             "mlpHC3_market_cv_001_2020-06-26.RData",
             "mlpHC4_market_cv_001_2020-07-02.RData")

load(cv_files[1])
load(cv_files[2])
load(cv_files[3])
load(cv_files[4])

cv_data = mlpSC_market_cv


cols = c("black", "darkblue", "blue", "cyan")

x = c(0.000001,0.00001,0.0001,0.001)
for (h in c(1:4)){
  for (lr in c(1:2)){    
    y = c(cv_data[[(h-1)*8+(lr-1)*4+4]][7,1],
          cv_data[[(h-1)*8+(lr-1)*4+3]][7,1],
          cv_data[[(h-1)*8+(lr-1)*4+2]][7,1],
          cv_data[[(h-1)*8+(lr-1)*4+1]][7,1])
    
    if((h==1)&(lr ==1)){
      plot(x,y, log= "x", ylim = c(0,0.0003), type = "b", xaxt = "n",
           yaxt = "n", col = cols[h], xlab="", ylab="MSE")
      axis(1, at = c(0.000001,0.00001,0.0001,0.001), 
           labels = c("0", "0.00001", "0.0001", "0.001"))
      axis(2, at = c(0.0001, 0.0002, 0.0003),
           labels = c("0.0001", "0.0002", "0.0003"))
      mtext(expression(lambda), side=1, line=2, las=1, cex =1.5)
    }else{
      
      lines(x,y, type = "b", col = cols[h], lty = lr)
      
    }
  }
}
leg =   c("h = 4| lr = 0.1",
          "h = 4| lr = 0.01",
          "h = 8| lr = 0.1",
          "h = 8| lr = 0.01",
          "h = 4,4| lr = 0.1",
          "h = 4,4| lr = 0.01",
          "h = 8,8| lr = 0.1",
          "h = 8,8| lr = 0.01")
legend(x = 0.0001, y = 0.00011, legend = leg,
       col =  c("black", "black", "darkblue", "darkblue",
                "blue",   "blue", "cyan", "cyan"),
       lty = c(1,2,1,2,1,2,1,2),
       bty = "n", cex = 0.8)












