library(DirichletReg)
library(abind)

# helpers
sum_by_z <- function(k, z, vec=NULL){
  if (is.null(vec)){ map_dbl(1:k, ~sum(z==.x)) }
  else { map_dbl(1:k, ~sum(vec[z==.x])) }
}

get_allocation_probs <- function(y,w,mu,sig2){
  # conditional distribution of z for each y
  ## likelihood times weight
  probs <- map(1:length(w), 
               ~( w[.x] * dnorm(y, mu[.x], sqrt(sig2[.x])))) %>% 
    abind(along=0) %>% t()
  probs/rowSums(probs)
}

# move type a
update_weights  <- function(z, k, 
                            delta){
  # alpha_i = delta + n_i
  alpha <- rep(delta, k) + sum_by_z(k,z)
  as.numeric(rdirichlet(1, alpha))
}

# move type b
update_mu <- function(y, k, z, mu, sig2, 
                      xi, kappa){
  
  # conditional dist
  n <- sum_by_z(k,z)
  y_sum <- sum_by_z(k, z, y) 
  cond_var <- 1/(n/sig2 + kappa)
  cond_mean <- (y_sum/sig2 + kappa*xi)*cond_var
  
  # propose and accept
  mu_proposal <- rnorm(k, cond_mean, cond_var)
  
  ## accept if and only if order is right
  if (all(order(mu_proposal)==1:k)){
    mu_proposal
  } else {
    mu
  }
}

update_sig2 <- function(y, k, z, mu, b, a){
  # conditional dist
  a_cond <- a + 0.5*sum_by_z(k,z)
  b_cond <- b + 0.5*sum_by_z(k, z, (y-mu[z])^2)
  
  1/rgamma(k, a_cond, b_cond)  
}

# move type c
update_z <- function(y, w, mu, sig2){ 
  probs <- get_allocation_probs(y,w,mu,sig2)
  map_int(1:length(y), ~sample(1:length(w), 1, prob = probs[.x,]))
}

# move type d
update_b <- function(sig2, a, g, h){
  
  a_cond <- g + a*length(sig2)
  b_cond <- h + sum(1/sig2)
  
  rgamma(1, a_cond, b_cond)
}


