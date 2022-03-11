# move type f

laccept_prob_birth <- function(w_new, k, k0, n, delta){
  # k is k for smaller configuration (pre-birth or post-death)
  
  prob_params <-
    -lbeta(k*delta, delta) +
    (delta-1)*log(w_new) +
    (n+k*(delta-1))*log(1-w_new) +
    log(k+1)
  
  hastings_ratio <- 
    -log(k0+1) -
    dbeta(w_new,1,k, log=T) +
    k*log(w_new)
  
  prob_params + hastings_ratio
}

death_step <- function(k, z, w, mu, sig2, delta){
  is_empty <- which(!(1:k %in% z))
  if (length(is_empty)==0){
    return(list(accept=FALSE))
  }
  
  j <- as.numeric(sample(as.character(is_empty), 1))
  laccept_prob <- 
    laccept_prob_birth(w[j], k, 
                       k0=length(is_empty), n=length(z), delta=delta)
  accept = (log(runif(1))<(-laccept_prob))
  
  out <- list(accept=accept)
  # if accepted, create new values
  if (accept){
    out$k <- k - 1
    out$w <- w[-j]/(1-w[j])
    out$mu <- mu[-j]
    out$sig2 <- sig2[-j]
    out$z <- ifelse(z>=j, z-1, z)
  }
  return(out)
}

birth_step <- function(k, z, w, mu, sig2,
                       a, b, xi, kappa, delta){
  # generate new parameters
  w_new <- rbeta(1, 1, k) 
  mu_new <- rnorm(1, xi, 1/sqrt(kappa))
  sig2_new <- 1/rgamma(1,a,b)
  
  # accept or reject
  num_empty <- k - length(unique(z))
  laccept_prob <- 
    laccept_prob_birth(w_new, k, 
                       k0=num_empty, n=length(z), delta=delta)
  accept = (log(runif(1))<laccept_prob)
  
  out <- list(accept=accept)
  # if accepted, create new values
  if (accept){
    out$k <- k + 1
    out$mu <- sort(c(mu, mu_new))
    
    j_new <- which(mu_out==mu_new)
    out$sig2 <- append(sig2, sig2_new, after=j_new-1)
    out$w <- append(w*(1-w_new), w_new, after=j_new-1)
    out$z <- ifelse(z>=j_new, z+1, z)
  }
  return(out)
}