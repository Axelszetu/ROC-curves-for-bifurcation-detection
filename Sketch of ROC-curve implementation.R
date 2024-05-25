#Sketch of computing ROC-curves
#Using variables and data as defined in AMOCestimation.Rmd
library(MASS) 
library(reshape2) 
library(reshape) 

#Starting with variance
no_paths <- dim(var.matrix)[2]
grid <- seq(from = 0.5, to = 1, by = 0.1)
thresholds <- qnorm(grid)
no_thresh <- length(thresholds)
positive_rates <- numeric(length = no_thresh)

for (i in (1:no_thresh)){
  thresh <- gam0 + thresholds[i] * sqrt(vargam0)
  no_triggered <- 0
  for(j in (1:no_paths)){
    path <- var.matrix[,j]
    triggered <- max(path > thresh)
    no_triggered <- no_triggered + triggered
  }
  positive_rates[i] <- no_triggered/no_paths
}

#This is written to take a matrix with a path in each column
#Can we meld xx.data to such a format?
#Is there an issue wit a matrix having 1000 columns?
#I mean, this one has more than a million rows, so recasting it would give us 10,000 by 1000.
#The number of cells doesn't increase

#We also need to do some amount of computation on the data:
#The current data is the raw observations.
#We need to compute running estimates of var and rho, and that data is to be fed to the positive rates computer
#I guess Susannes code for computing Var.matrix also takes data of the format we have

EWS.data.H0 <- dcast(data = xx.data.H0, formula = "repetition")

var.matrix.H0 = matrix(-1, ncol = rep, nrow = n - Tw/dt)
rho.matrix.H0 = matrix(-1, ncol = rep, nrow = n - Tw/dt)
for(i in 1:(n - Tw/dt)){
  for(j in 1:rep){
    data.i = X[i:(i+Tw/dt),j] ## Get data in running window
    data.i = data.i[data.i > -2.3] ## Remove data if tipped
    time.i = (1:length(data.i))*dt ## Get time interval for data in running window
    temp   = lm(data.i ~ time.i)   ## Get linear trend
    data.i = data.i - temp$fitted.values ## Subtract linear trend
    var.matrix.H0[i, j] = var(data.i) ## Variance estimate in running window
    rho.matrix.H0[i, j] = acf(data.i, lag.max = 1, plot = FALSE)$acf[2] ## Autocorrelation estimate in running window
  }
}