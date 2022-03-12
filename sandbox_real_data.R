library(tidyverse)
library(yaml)

configs <- read_yaml("config.yaml")

for (f in list.files("functions")){
  print(f)
  source(file.path("functions",f))
}

y <- sort(MASS::galaxies)

k_init <- 10
w_init <- 1/k_init
z_init <- ceiling(k_init*(1:length(y))/length(y))
mu_init <- quantile(y, seq(0,1,length=k_init))
sig2_init <- rep(var(y), k_init)

# test run
set.seed(12345)
out <- run_mcmc(y, k_init, w_init, z_init, mu_init, sig2_init,
                sweeps=1000, burn=1000)

# interpret results

## posterior for k
prop.table(table(out$k))
plot(out$k, type="l")

# acceptance fpr split and combine
accept_stats <- tibble(label = c("split","combine",
                                 "birth","death"),
                       res = list(out$accept_split, 
                                  out$accept_combine,
                                  out$accept_birth,
                                  out$accept_death)) %>%
  mutate(x=map(res, sum),
         n=map(res, length),
         p=map(res, mean)) %>%
  select(-res) %>%
  unnest(c(x,n,p))
accept_stats

## posterior means
post_means_by_k <-
  tibble(k=unique(out$k)) %>%
  mutate(post_means = map(k, post_means_fixed_k, out=out)) %>%
  unnest(post_means)

post_means_by_k %>% filter(k==3)

post_means_by_k %>%
  filter(mu==max(mu) | mu==min(mu))

## cluster plot
make_cluster_plot(out)

1## Predictive densities (in progress)
plot_posterior_densities(y, out, k=3:5)

## How often are components empty?
z_count <- apply(out$z, 1, function(x){length(unique(x))})
table(out$k, out$k-z_count)
