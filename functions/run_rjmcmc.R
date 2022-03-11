# run MCMC

run_mcmc <- function(y, k, w_init, z_init, mu_init, sig2_init,
                     sweeps = 1000, burn = 1000){
  
  # constants from data for priors
  M <- (max(dat$y) + min(dat$y))/2    # midpoint
  R2 <- (max(dat$y) - min(dat$y))^2   # range
  
  # outputs
  k <- rep(NA, sweeps)
  z <- matrix(nrow=sweeps, ncol=length(y))
  w <- vector(mode = "list", length = sweeps)
  mu <- vector(mode = "list", length = sweeps)
  sig2 <- vector(mode = "list", length = sweeps)
  
  accept_split <- c()
  accept_combine <- c()
  
  # current values
  k_curr <- k_init
  z_curr <- z_init
  w_curr <- w_init
  mu_curr <- mu_init
  sig2_curr <- sig2_init
  b_curr <- b_init
  
  for (i in 1:(burn+sweeps)){
    # step a
    w_curr <- update_weights(z_curr, k_curr, del=configs$delta)
    
    # seto b
    mu_curr <- update_mu(y, k_curr, z_curr, mu_curr, sig2_curr,   
                         xi=M, kappa=configs$kappa_mult/R2)
    sig2_curr <- update_sig2(y, k_curr, z_curr, mu_curr, 
                             b_curr, a=configs$a)
    
    # step c
    z_curr <- update_z(y, w_curr, mu_curr, sig2_curr)
    
    # step d
    b_curr <- update_b(sig2_curr, 
                         a=configs$a, 
                         g=configs$g,
                         h=configs$h_mult/R2)
    
    # step e
    ## choose whether to do split or combine step
    choose_split <- (sample(0:1, 1)==1) & k_curr!=configs$k_max
    if (choose_split){
      split_comb <- split_step(y, k_curr, z_curr, 
                          w_curr, mu_curr, sig2_curr, 
                          a=configs$a, b=b_curr, 
                          xi=M, kappa=configs$kappa_mult/R2,
                          delta=configs$delta)
      accept_move_e <- split_comb$accept
      accept_split <- c(accept_split, accept_move_e)
    } else {
      split_comb <- combine_step(y, k_curr, z_curr, 
                          w_curr, mu_curr, sig2_curr, 
                          a=configs$a, b=b_curr, 
                          xi=M, kappa=configs$kappa_mult/R2,
                          delta=configs$delta)
      accept_move_e <- split_comb$accept
      accept_combine <- c(accept_combine,accept_move_e)
    }
    if (accept_move_e){
      k_curr <- split_comb$k
      z_curr <- split_comb$z
      w_curr <- split_comb$w
      mu_curr <- split_comb$mu
      sig2_curr <- split_comb$sig2
    }
    
    # step f
    
    # store outputs
    if (i>burn){
      idx <- i-burn
      k[idx] <- k_curr
      z[idx,] <- z_curr
      w[[idx]] <- w_curr
      mu[[idx]] <- mu_curr
      sig2[[idx]] <- sig2_curr
    }
  }
  
  list(k=k, z=z, w=w, mu=mu, sig2=sig2,
       accept_split=accept_split, accept_combine=accept_combine)
}
