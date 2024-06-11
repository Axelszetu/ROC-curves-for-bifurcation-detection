#This script is a sketch of how we implement computation of T_w, the window of observation for running estimation of EWS.
A <- 1
alpha_computer <- function(lambda, A){
  alpha <- 2*sqrt(A*abs(lambda))
  alpha
}

alpha_computer(-3, A)
#Okay, good shit.

T_w_var_computer <- function(lambda_0, lambda_t, q, A = 1){
  alpha_0 <- alpha_computer(lambda_0, A)
  alpha_t <- alpha_computer(lambda_t, A)
  T_w <- 2*q^2*((alpha_t/sqrt(alpha_0) + alpha_0/sqrt(alpha_t))/(alpha_0-alpha_t))^2
  T_w
}
T_w_var_computer(lambda_0 = -1.5, lambda_t = -1.5/2, q = 1, A = 0.75)
#Okay, I get som reasonable numbers. Good jobbo.

T_w_rho_computer <- function(lambda_0, lambda_t, q, A, rho){
  alpha_0 <- alpha_computer(lambda_0, A)
  alpha_t <- alpha_computer(lambda_t, A)
  T_w <- 2*q^2*((sqrt(alpha_0) + sqrt(alpha_t))/(alpha_0 - alpha_t))^2*rho^(-2)
  T_w
}
T_w_rho_computer(lambda_0 = -1.5, lambda_t = -1.5/2, q = 1, A = 0.75, rho = rho0)
