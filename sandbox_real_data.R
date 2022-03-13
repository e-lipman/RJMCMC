library(tidyverse)
library(yaml)
library(tictoc)

configs <- read_yaml("config.yaml")

for (f in list.files("functions")){
  print(f)
  source(file.path("functions",f))
}

dat_name = "enzyme"

if (dat_name == "enzyme"){
  library(BNPdensity)
  data(enzyme)
  y <- sort(enzyme)
} else if (dat_name == "acidity"){
  library(multimode)
  data(acidity)
  y <- sort(acidity)
} else if (dat_name == "galaxy"){
  y <- sort(MASS::galaxies)/1000
}

hist(y, breaks=50)

k_init <- 2
w_init <- rep(1/k_init, k_init)
z_init <- ceiling(k_init*(1:length(y))/length(y))
mu_init <- quantile(y, seq(0,1,length=k_init))
sig2_init <- rep(max(y)-min(y), k_init)

# test run
tic()
set.seed(12345)
out <- run_mcmc(y, k_init, w_init, z_init, mu_init, sig2_init,
                sweeps=1000, burn=1000)
toc()

# interpret results

## posterior for k
prop.table(table(out$k))
plot(out$k, type="l")

## Predictive densities (in progress)
plot_posterior_densities(y, out, k=1:9)

# acceptance fpr split and combine
get_accept_stats(out, agg=T)

## posterior means
post_means_by_k <-
  tibble(k=unique(out$k)) %>%
  mutate(post_means = map(k, post_means_fixed_k, out=out)) %>%
  unnest(post_means)

post_means_by_k %>% filter(k==2)

## cluster plot
make_cluster_plot(out)

## How often are components empty?
z_count <- apply(out$z, 1, function(x){length(unique(x))})
table(out$k, out$k-z_count)
mean(out$k-z_count)

