#Sketch of computing ROC-curves
#Using variables and data as defined in AMOCestimation.Rmd

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

