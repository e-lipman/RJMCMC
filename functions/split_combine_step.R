# helpers 
laccept_prob_split <- function(yj,
                              # params from small k model
                              k_sm, 
                              w_sm, mu_sm, sig2_sm,
                              # params from large k model
                              z_lg, 
                              w_lg, mu_lg, sig2_lg,
                              # helpful quantties
                              u, lp_alloc,
                              #hyperparams,
                              a, b, xi, kappa, delta){
  # note: arguments with suffix sm apply to the setting 
  #       with smaller k (i.e. before splitting) and arguments 
  #      with suffix lg apply to larger k (after splitting)
  
  # helpful quantities
  nj <- sum_by_z(2, z_lg)
  
  # Note: all quantities below are on the LOG SCALE
  #       and most are RATIOS of lg/sm
  
  # part1: P(x')/P(x)
  LR <- sum(dnorm(yj, mu_lg[z_lg], sig2_lg[z_lg], log=T)) - 
    sum(dnorm(yj, mu_sm, sig2_sm, log=T))
  
  p_wz_given_k <- 
    (delta - 1 + nj[1])*log(w_lg[1]) +
    (delta - 1 + nj[2])*log(w_lg[2]) -
    (delta - 1 + sum(nj))*log(w_sm) -
    lbeta(delta, delta*k_sm)
  
  p_mu_given_k <-
    log(k_sm+1) + # factor from permutattions
    dnorm(mu_lg[1], xi, sqrt(1/kappa), log=T) +
    dnorm(mu_lg[2], xi, sqrt(1/kappa), log=T) -
    dnorm(mu_sm, xi, sqrt(1/kappa), log=T)
  
  p_sig2_given_k <-
    a*log(b) - lgamma(a) -
    (a+1)*log(prod(sig2_lg)/sig2_sm) -
    b*(1/sum(1/sig2_lg)-sig2_sm)
  
  part1_sum <- LR + p_wz_given_k + 
    p_mu_given_k + p_sig2_given_k
  
  # part 2: Hastings ratio
  ## note: our implementation assumes bk=dk for all valid moves
  
  move_prob <- 
    - lp_alloc -
    sum(dbeta(u, c(2,2,1), c(2,2,1), log=T))
  
  jacobian <-
    log(w_sm*abs(mu_lg[1]-mu_lg[2])*prod(sig2_lg)) -
    log(u[2]*(1-u[2]^2)*u[3]*(1-u[3])*sig2_sm)
  
  part2_sum <- move_prob + jacobian
  
  # putting it together
  part1_sum + part2_sum
}

# move type e
combine_step <- function(y, k, z, w, mu, sig2,
                         a, b, xi, kappa, delta){
  # choose component to split
  j1 <- sample(1:(k-1), 1)
  j2 <- j1+1
  
  # find values for new component
  w_new <- sum(w[j1:j2])
  mu_new <- weighted.mean(mu[j1:j2], w[j1:j2])
  sig2_new <- weighted.mean(mu[j1:j2]^2 + sig2[j1:j2], 
                            w[j1:j2]) - mu_new^2
  
  # acompute transition variables u
  u1 <- w[j1]/w_new
  u2 <- sqrt(w[j1]/w[j2])*(mu_new-mu[j1])/sqrt(sig2_new)
  u3 <- (sig2[j1]/sig2_new)*(w[j1]/w_new)*(1/(1-u2^2))
  if (any(c(u1,u2,u3)>1))browser()
  
  # compute arguments to acceptance prob function 
  yj <- y[z %in% c(j1,j2)]
  z_lg <- z[z %in% c(j1,j2)] - j1
  probs <- get_allocation_probs(yj,
                                w[j1:j2],mu[j1:j2],sig2[j1:j2])
  lp_alloc <- sum(dbinom(z_lg, 1, probs[,2],log=T))
  
  # acceptance
  laccept_prob <- laccept_prob_split(
    yj,  k_sm=k-1, 
    w_sm=w_new, mu_sm=mu_new, sig2_sm=sig2_new, 
    z_lg=z_lg, 
    w_lg=w[j1:j2], mu_lg=mu[j1:j2], sig2_lg=sig2[j1:j2], 
    u=c(u1,u2,u3), lp_alloc=lp_alloc, 
    a=a, b=b, xi=xi, kappa=kappa, delta=delta)
  accept = (log(runif(1))<laccept_prob)
  out <- list(accept=accept)
  
  # if accepted, create new values if needed
  if (accept){
    out$k <- k -1
    out$z <- ifelse(z>j1, z-1, z)
    
    out$w <- append(w[-c(j1:j2)], w_new, after=j1-1)
    out$mu <- append(mu[-c(j1:j2)], mu_new, after=j1-1)
    out$sig2 <- append(sig2[-c(j1:j2)], sig2_new, after=j1-1)
  }
  
  return(out) 
}

split_step <- function(y, k, z, w, mu, sig2,
                       a, b, xi, kappa, delta){
  # choose component to merge and draw latent variables
  j <- sample(1:k, 1)
  u <- rbeta(3, c(2,2,1), c(2,2,1))
  
  # compute params for new component
  w1 <- w[j]*u[1]
  w2 <- w[j] - w1
  mu1 <- mu[j] - u[2]*sqrt(sig2[j])*sqrt(w2/w1)
  mu2 <- mu[j] + u[2]*sqrt(sig2[j])*sqrt(w1/w2)
  sig21 <- u[3]*(1-u[2]^2)*sig2[j]*(w[j]/w1)
  sig22 <- (1-u[3])*(1-u[2]^2)*sig2[j]*(w[j]/w2)
  
  # check that mu satisfies the ordering conditions
  conditions <- logical(0)
  if (j>1){
    conditions <- c(conditions, mu1 > mu[j-1])
  }
  if (j<k){ 
    conditions <- c(conditions, mu2 < mu[j+1])
  }
  if (!all(conditions)){
    return(list(accept=F))  # EARLY RETURN
  }
  
  # compute allocation probabilities 
  probs <- get_allocation_probs(y[z==j],c(w1,w2),
                                c(mu1,mu2),c(sig21,sig22))
  alloc <- rbinom(nrow(probs), 1, probs[,2])
  lp_alloc <- sum(dbinom(alloc, 1, probs[,2],log=T))
  
  # acceptance
  laccept_prob <- laccept_prob_split(
    y[z==j],  k_sm=k, 
    w_sm=w[j], mu_sm=mu[j], sig2_sm=sig2[j], 
    z_lg=alloc+1, 
    w_lg=c(w1,w2), mu_lg=c(mu1,mu2), sig2_lg=c(sig21, sig22), 
    u=u, lp_alloc=lp_alloc, 
    a=a, b=b, xi=xi, kappa=kappa, delta=delta)
  accept = (log(runif(1))<laccept_prob)
  out <- list(accept=accept)
  
  # if accepted, create new values
  if (accept){
    out$k <- k + 1
    out$z <- ifelse(z>j, z+1, z)
    out$z[z==j] = out$z[z==j] + alloc
  
    out$w <- append(w[-c(j)], c(w1,w2), after=j-1)
    out$mu <- append(mu[-c(j)], c(mu1,mu2), after=j-1)
    out$sig2 <- append(sig2[-c(j)], c(sig21,sig22), after=j-1)
  }
  
  return(out)
}

