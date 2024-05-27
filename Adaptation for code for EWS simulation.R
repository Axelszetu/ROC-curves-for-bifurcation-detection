### Simulation with ramping
### Simulations adaoted from EWS section

########################################################
## Figures with realizations under changing lambda    ##
## Running estimation of variance and autocorrelation ##
########################################################

## Parameter values used in the simulation
tau  = 200  ## 
l0   = -3   ## Baseline control parameter lambda
sigm = 0.25  ## 
s2   = sigm^2  ## 

## Settings for simulation and estimation
rep = 10   ## Number of repetitions of simulated traces
Tw  = 50   ## Window size for variance and autocorrelation estimation
t0  = 100  ## Start of change in control parameter lambda
dt  = 0.1  ## Time step between observations
nloop = 10 ## Substeps between observation steps to decrease discretization error 
alpha0 = 2*sqrt(-l0)      ## Baseline alpha in OU approximation
rho0 = exp(-alpha0 * dt)    ## Baseline autocorrelation in OU approximation
gam0 = s2/(2*alpha0)        ## Baseline variance in OU approximation
vargam0 = 2*gam0^2/(alpha0 * Tw)  ## Baseline variance of variance estimator
varrho0 = 2*alpha0*dt^2/Tw  ## Baseline variance of autocorrelation estimator
q95 = qnorm(0.95)  ## 95% quantile of normal distribution

X0   = sqrt(-l0)  ## Starting value for simulation
Tend = t0 + tau   ## Length of simulations
time = seq(0, Tend, dt) - t0  ## Time for simulations
n    = length(time)  ## Number of observation points in simulations
lt   = c(rep(l0, t0/dt), (l0*(1 - (0:(tau/dt))/(tau/dt)))) ## Evolution of lambda
mu   = sqrt(-lt)     ## Evolution of mu
alpha = 2 * sqrt(-lt)    ## Evolution of alpha
rho  = exp( - alpha * dt)  ## Evolution of autocorrelation
gam  = s2/(2*alpha)        ## Evolution of variance
vargam = 2*gam^2/(alpha*Tw) ## Evolution of variance of var-estimator
varrho = 2*alpha*dt^2/Tw   ## Evolution of variance of corr-estimator

## Simulations
X = matrix(-10, ncol = rep, nrow = n)

for (i in 1:rep){
  #Initial condition is random from stationary distribution
  x0 = X0 + rnorm(1, mean = 0, sd = sigm/sqrt(2*alpha0)) 
  
  xx = X.traj(sigma = sigm, lambda0 = l0, tau = tau, m = 0, a = 1,
              T0 = t0, X0 = x0, dt = dt/nloop, 
              Ymax = nloop*n-1)
  nx = length(xx$X) ##If tipping happens before, trajectory is shorter than nt
  xxx = xx$X[seq(1, nx, nloop)] ##Only keep points at observed times
  nxx = min(n, length(xxx))
  X[1:nxx, i] = xxx[1:nxx]
}

plotdata = data.frame(X = c(X), time = rep(time, rep), rep = rep(1:rep, each = n))
plotdata$X[plotdata$X < -9.9] = NA

ggplot(data = plotdata, aes(x = time, y = pmax(X,-2.2), group = as.factor(rep))) +
  geom_line(color = "gray") +
  geom_line(data = data.frame(X = mu, time = time), color = "blue", size = 0.8) +
  geom_line(data = data.frame(X = -mu, time = time), color = "blue", size = 0.8, linetype = "dashed") +
  ylim(-2.2, 2.2) +
  ylab("Realizations") + xlab("Time (years)")

### Simulations without ramping
### Simulations adapted from EWS section

X.traj <- function(sigma = 0.1, lambda0 = -2, tau = 1000, m = 0, a = 1,  
                   T0 = 0, X0 = sqrt(2), dt = 0.1, Ymax = 1000000){
  ##T0: Time before ramping starts
  ##Ymax: Max number of simulated points, if tipping does not happen
  Y = 0  ## Counting integration steps up to tipping
  xbarrier = m - 2 ## One smaller than crossing point at start
  Xtraj = X0
  X = X0
  ## Simulation during stationary period, constant lambda
  while(X > xbarrier & Y < T0/dt){
    X = S.onestep(sigma = sigma, lambda = lambda0,  
                  m = m, a = a, X0 = X, dt = dt)
    Xtraj = c(Xtraj, X)
    Y = Y+1
  }
  Y = Y*dt
  return(list(FPT = Y, X = Xtraj))
}

########################################################
## Figures with realizations under changing lambda    ##
## Running estimation of variance and autocorrelation ##
########################################################

## Parameter values used in the simulation
tau  = 0  ## 
l0   = -3   ## Baseline control parameter lambda
sigm = 0.25  ## 
s2   = sigm^2  ## 

## Settings for simulation and estimation
rep = 10   ## Number of repetitions of simulated traces
Tw  = 50   ## Window size for variance and autocorrelation estimation
t0  = 300  ## Start of change in control parameter lambda
dt  = 0.1  ## Time step between observations
nloop = 10 ## Substeps between observation steps to decrease discretization error 
alpha0 = 2*sqrt(-l0)      ## Baseline alpha in OU approximation
rho0 = exp(-alpha0 * dt)    ## Baseline autocorrelation in OU approximation
gam0 = s2/(2*alpha0)        ## Baseline variance in OU approximation
vargam0 = 2*gam0^2/(alpha0 * Tw)  ## Baseline variance of variance estimator
varrho0 = 2*alpha0*dt^2/Tw  ## Baseline variance of autocorrelation estimator
q95 = qnorm(0.95)  ## 95% quantile of normal distribution

X0   = sqrt(-l0)  ## Starting value for simulation
Tend = t0 + tau   ## Length of simulations
time = seq(0, Tend, dt) - t0  ## Time for simulations
n    = length(time)  ## Number of observation points in simulations
lt   = c(rep(l0, t0/dt), (l0*(1 - (0:(tau/dt))/(tau/dt)))) ## Evolution of lambda
mu   = sqrt(-lt)     ## Evolution of mu
alpha = 2 * sqrt(-lt)    ## Evolution of alpha
rho  = exp( - alpha * dt)  ## Evolution of autocorrelation
gam  = s2/(2*alpha)        ## Evolution of variance
vargam = 2*gam^2/(alpha*Tw) ## Evolution of variance of var-estimator
varrho = 2*alpha*dt^2/Tw   ## Evolution of variance of corr-estimator

## Simulations
X = matrix(-10, ncol = rep, nrow = n)

for (i in 1:rep){
  #Initial condition is random from stationary distribution
  x0 = X0 + rnorm(1, mean = 0, sd = sigm/sqrt(2*alpha0)) 
  
  xx = X.traj(sigma = sigm, lambda0 = l0, tau = tau, m = 0, a = 1,
              T0 = t0, X0 = x0, dt = dt/nloop, 
              Ymax = nloop*n-1)
  nx = length(xx$X) ##If tipping happens before, trajectory is shorter than nt
  xxx = xx$X[seq(1, nx, nloop)] ##Only keep points at observed times
  nxx = min(n, length(xxx))
  X[1:nxx, i] = xxx[1:nxx]
}

plotdata = data.frame(X = c(X), time = rep(time, rep), rep = rep(1:rep, each = n))
plotdata$X[plotdata$X < -9.9] = NA

ggplot(data = plotdata, aes(x = time, y = pmax(X,-2.2), group = as.factor(rep))) +
  geom_line(color = "gray") +
  geom_line(data = data.frame(X = mu, time = time), color = "blue", size = 0.8) +
  geom_line(data = data.frame(X = -mu, time = time), color = "blue", size = 0.8, linetype = "dashed") +
  ylim(-2.2, 2.2) +
  ylab("Realizations") + xlab("Time (years)")

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

plotdata.var.rho.H0 = data.frame(var = c(var.matrix.H0), rho = c(rho.matrix.H0),
                              time = - t0 + Tw/2 + rep(dt*seq(1,n - Tw/dt,1), rep), 
                              rep = rep(1:rep, each = n - Tw/dt))

## Variance
ggplot(data = plotdata.var.rho.H0, aes(x = time, y = var, group = as.factor(rep))) + 
  geom_line(color = "gray", size = 0.3) +
  geom_line(data = data.frame(var = gam, time = time), color = "blue", size = 1) +
  geom_vline(xintercept = tau, size = 1.5) + 
  geom_hline(yintercept = gam0, size = 1) + 
  geom_hline(yintercept = gam0 + q95 * sqrt(vargam0), size = 0.7, linetype = "dashed") + 
  geom_hline(yintercept = gam0 - q95 * sqrt(vargam0), size = 0.7, linetype = "dashed") + 
  geom_line(data = data.frame(var = gam + q95 * sqrt(vargam), time = time), size = 0.7, linetype = "dashed") +
  geom_line(data = data.frame(var = gam - q95 * sqrt(vargam), time = time), size = 0.7, linetype = "dashed") +
  ylim(0, 0.05) + xlim(-Tw/2,tau) +
  ylab("Variance") + xlab("Time") +
  geom_segment(aes(x = 0,  y = 0.025, xend = Tw, yend = 0.025), lineend = "butt", size = 1) +
  geom_segment(aes(x = 0,  y = 0.024, xend = 0,  yend = 0.026), lineend = "butt", size = 1) +
  geom_segment(aes(x = Tw, y = 0.024, xend = Tw, yend = 0.026), lineend = "butt", size = 1) +
  annotate("text", x=Tw/2, y=0.028, label= "window size", size = 5)  

## Autoorrelation
ggplot(data = plotdata.var.rho.H0, aes(x = time, y = rho, group = as.factor(rep))) +
  geom_line(color = "gray", size = 0.3) +
  geom_line(data = data.frame(rho = rho, time = time), color = "blue", size = 1) +
  geom_vline(xintercept = tau, size = 1.5) + 
  geom_hline(yintercept = rho0, size = 1) + 
  geom_hline(yintercept = rho0 + q95 * sqrt(varrho0), size = 0.7, linetype = "dashed") + 
  geom_hline(yintercept = rho0 - q95 * sqrt(varrho0), size = 0.7, linetype = "dashed") + 
  geom_line(data = data.frame(rho = rho - q95 * sqrt(varrho), time = time), 
            size = 0.7, linetype = "dashed") +
  geom_line(data = data.frame(rho = rho + q95 * sqrt(varrho), time = time), 
            size = 0.7, linetype = "dashed") +
  ylim(0.6, 1) + xlim(-Tw/2,tau) +
  ylab("Autocorrelation") + xlab("Time") +
  geom_segment(aes(x = 0,  y = 0.9,  xend = Tw, yend = 0.9),  lineend = "butt", size = 1) +
  geom_segment(aes(x = 0,  y = 0.89, xend = 0,  yend = 0.91), lineend = "butt", size = 1) +
  geom_segment(aes(x = Tw, y = 0.89, xend = Tw, yend = 0.91), lineend = "butt", size = 1) +
  annotate("text", x=Tw/2, y=0.925, label= "window size", size = 5)

#Where we stand at the date of bifurcation
no_paths <- dim(var.matrix.H0)[2]
grid <- seq(from = 0.95, to = 1, length.out = 10)
thresholds <- qnorm(grid)
no_thresh <- length(thresholds)
false_positive_rates <- numeric(length = no_thresh)
thresh_scaled <- gam0 + thresholds*sqrt(vargam0)

for (i in (1:no_thresh)){
  thresh <- gam0 + thresholds[i] * sqrt(vargam0)
  no_triggered <- 0
  for(j in (1:no_paths)){
    path <- var.matrix.H0[,j]
    triggered <- max(path > thresh)
    no_triggered <- no_triggered + triggered
  }
  false_positive_rates[i] <- no_triggered/no_paths
}

false_positive_rates

roc.plot.data <- data.frame(false_positive_rates, true_positive_rates)
ggplot(data = roc.plot.data, aes(x = false_positive_rates, y = true_positive_rates)) +
  geom_step() +
  xlim(0,1) + ylim(0,1)

#Now, where we stand at the time where ramping starts, t=

grid <- seq(from = 0.85, to = 1, length.out = 100)
thresholds <- qnorm(grid)
no_thresh <- length(thresholds)
thresh_scaled <- gam0 + thresholds*sqrt(vargam0)

true_positive_rates <- numeric(length = no_thresh)
for (i in (1:no_thresh)){
  thresh <- gam0 + thresholds[i] * sqrt(vargam0)
  no_triggered <- 0
  for(j in (1:no_paths)){
    path <- var.matrix[1:1500,j]
    triggered <- max(path > thresh)
    no_triggered <- no_triggered + triggered
  }
  true_positive_rates[i] <- no_triggered/no_paths
}


false_positive_rates <- numeric(length = no_thresh)
for (i in (1:no_thresh)){
  thresh <- gam0 + thresholds[i] * sqrt(vargam0)
  no_triggered <- 0
  for(j in (1:no_paths)){
    path <- var.matrix.H0[1:1500,j]
    triggered <- max(path > thresh)
    no_triggered <- no_triggered + triggered
  }
  false_positive_rates[i] <- no_triggered/no_paths
}

roc.plot.data <- data.frame(false_positive_rates, true_positive_rates)
ggplot(data = roc.plot.data, aes(x = false_positive_rates, y = true_positive_rates)) +
  geom_step() +
  xlim(0,1) + ylim(0,1)
