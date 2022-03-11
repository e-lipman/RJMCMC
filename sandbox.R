library(tidyverse)
library(yaml)

configs <- read_yaml("config.yaml")

for (f in list.files("functions")){
  print(f)
  source(file.path("functions",f))
}

# generate toy example
N <- 200
w <- rep(1/3,3)
mu <- c(0,5,10)
sig2 <- c(.5,1,.5)

dat <- tibble(z = sample(1:3, N, prob=w, rep=T)) %>%
  mutate(mu = mu[z],
         sig2 = sig2[z],
         y = rnorm(n(), mu, sqrt(sig2))) %>%
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
k_init=3

set.seed(123)
b_init <- rgamma(1, configs$g, configs$h_mult/R2)

# test run
set.seed(12345)
out <- run_mcmc(dat$y, k_init, w_init, z_init, mu_init, sig2_init,
                sweeps=1000, burn=1000)

# interpret results

## posterior for k
prop.table(table(out$k))
plot(out$k)

# acceptance fpr split and combine
accept_stats <- tibble(label = c("split","combine"),
                 res = list(out$accept_split, 
                            out$accept_combine)) %>%
  mutate(x=map(res, sum),
         n=map(res, length),
         p=map(res, mean)) %>%
  select(-res) %>%
  unnest(c(x,n,p))
accept_stats

## posterior means
post_means_fixed_k(out, k=4)
 
## cluster plot
make_cluster_plot(out)
 
## Predictive densities (in progress)
plot_posterior_densities(dat$y, out, k=3:5)

## How often are components empty?
z_count <- apply(out$z, 1, function(x){length(unique(x))})
table(out$k, z_count)
