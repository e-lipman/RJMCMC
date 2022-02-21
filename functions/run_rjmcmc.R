# run MCMC

run_mcmc <- function(y, k, w_init, z_init, mu_init, sig2_init,
                     sweeps = 1000, burn = 1000){
  
  
  # constants from data for priors
  M <- (max(dat$y) + min(dat$y))/2    # midpoint
  R2 <- (max(dat$y) - min(dat$y))^2   # range
  
  # outputs
  z <- matrix(nrow=sweeps, ncol=length(y))
  w <- vector(mode = "list", length = sweeps)
  mu <- vector(mode = "list", length = sweeps)
  sig2 <- vector(mode = "list", length = sweeps)
  
  # current values
  z_curr <- z_init
  w_curr <- w_init
  mu_curr <- mu_init
  sig2_curr <- sig2_init
  b_curr <- b_init
  
  for (i in 1:(burn+sweeps)){
    # step a
    w_curr <- update_weights(z_curr, k, del=configs$delta)
    
    # seto b
    mu_curr <- update_mu(y, k, z_curr, mu_curr, sig2_curr,   
                         xi=M, kappa=configs$kappa_mult/R2)
    sig2_curr <- update_sig2(y, k, z_curr, mu_curr, b_curr, a=configs$a)
    
    # step c
    z_curr <- update_z(y, w_curr, mu_curr, sig2_curr)
    
    # step d
    b_curr <- update_b(sig2_curr, 
                         a=configs$a, 
                         g=configs$g,
                         h=configs$h_mult/R2)
    
    # store outputs
    if (i>burn){
      idx <- i-burn
      z[idx,] <- z_curr
      w[[idx]] <- w_curr
      mu[[idx]] <- mu_curr
      sig2[[idx]] <- sig2_curr
    }
  }
  
  list(z=z, w=w, mu=mu, sig2=sig2)
}
