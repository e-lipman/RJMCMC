library(tidyverse)
library(yaml)

configs <- read_yaml("config.yaml")

for (f in list.files("functions")){
  print(f)
  source(file.path("functions",f))
}

# generate toy example
N <- 200
mu <- c(0,5,10)
sig2 <- c(.5,1,.5)

dat <- tibble(z = sample(1:3, N, rep=T)) %>%
  mutate(mu = mu[z],
         sig = sig[z],
         y = rnorm(n(), mu, sig)) %>%
  arrange(y)

hist(dat$y, breaks=50)

# run MCMC with fixed k=3

## hyperprior constants from data
xi <- (max(dat$y) + min(dat$y))/2
R2 <- (max(dat$y) - min(dat$y))^2

## initial values
w_init <- c(.3,.4,.3)
z_init <- rep(c(1:3), each=ceiling(nrow(dat)/3))[1:nrow(dat)]
mu_init <- 1:3
sig2_init <- c(1,1,1)
k=3

set.seed(123)
b_init <- rgamma(1, configs$g, configs$h_mult/R2)

# test run
out <- run_mcmc(dat$y, k, w_init, z_init, mu_init, sig2_init)

# interpret results
make_cluster_plot(out)

post_means_fixed_k(out)

post_pred <- post_pred_dist(out)
hist(dat$y, breaks=50, probability=T)
hist(post_pred, breaks=50, probability=T)

