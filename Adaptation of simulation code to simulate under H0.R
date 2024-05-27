#The following script is a modification of the code in Supplementary Material of Ditlevsen et al 2023
#It is used to generate realizations of the sst model under stationaly dynamics(no change in lambda over time)

### SIMULATIONS ###
set.seed(123) ## For reproducibility
nsim = 1000   ## Number of simulations

#Define parameters
t0 <- 2020.99 #Postpoining ramping until simulation is over, thus getting unramped dynamics in the whole period

T0      = t0-1870      ## Time length of stationary dynamics (lambda constant)
sig     = sqrt(sigma2) ## Infinitisimal standard deviation: sigma
time    = AMOC.data$time 
n       = length(time)
lam.seq = lambda0*(1-pmax(0,(time-t0))/tau) ## Time varying lambda
alpha.seq = 2*sqrt(-a*lam.seq) ## Time varying alpha
mu.seq  = m + sqrt(-lam.seq/a) ## Time varying mu
var0    = sigma2/(2*alpha0)    ## Time varying variance
Delta_t = 1/12                 ## Time step in years (one measurement per month)
nloop   = 20                   ## Number of simulation steps between observation steps 
tp      = 2020 - t0            ## Time length of non-stationary dynamics (lambda linearly increasing)
nt      = ceiling(tp/Delta_t) ## Total number of points from t_0 to today

##Simulate process with lambda -> 0, starting at time t_0
xx = matrix(NA, n, nsim)

for (isim in 1:nsim){
  #Initial condition is random from stationary distribution
  x0 = m + sqrt(-lambda0/a) + rnorm(1, mean = 0, sd = sig/sqrt(2*alpha0)) 
  
  xxx = X.traj(sigma = sig, lambda0 = lambda0, tau = tau, m = m, a = a,
               T0 = T0, X0 = x0, dt = Delta_t/nloop, 
               Ymax = nloop*n-1)
  nx = length(xxx$X) ##If tipping happens before, trajectory is shorter than nt
  xxxx = xxx$X[seq(1, nx, nloop)] ##Only keep points at observed times
  nxx = length(xxxx)
  xx[1:(nxx), isim] = xxxx
}

xxtime <- time[1:n]
xx.data.H0 = data.frame(X = c(xx), time = rep(xxtime, nsim), 
                     repetition = rep(1:nsim, each = n),
                     model = rep(m + sqrt(-lam.seq[1:n]/a), nsim))

## End of simulation of nsim paths

## Save for later use
save(xx.data.H0, nsim, Delta_t, file = "SimulatedTracesH0.Rdata")




