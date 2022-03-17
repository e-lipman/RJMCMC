library(tidyverse)
library(yaml)
library(tictoc)

configs <- read_yaml("config.yaml")

for (f in list.files("functions")){
  print(f)
  source(file.path("functions",f))
}

dat_name = configs$dat_name

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
} else{
  stop("Bad dat_name")
}

hist(y, breaks=50)

# run sampler
tic()
set.seed(71835)
out <- run_rjmcmc(y, sweeps=configs$sweeps, burn=configs$burnin,
                  progress = 1000)
toc()

# results

## posterior for k
prop.table(table(out$k))      # posterior for k
plot(out$k, type="l")         # traceploft for k

k_post <- tibble(k=names(table(out$k)),
                 prob=table(out$k)) %>%
  mutate(prob = prob/sum(prob))
if (configs$save){
  saveRDS(k_post,
          file.path("results", paste0(dat_name, "_kpost.RDS")))
}

## cluster plot
cluster_plot <- make_cluster_plot(out)
cluster_plot

if (configs$save){
  saveRDS(cluster_plot,
          file.path("results", paste0(dat_name, "_cluster.RDS")))
  ggsave(width = 5, 
         height = 3, 
         file.path("figures", paste0(dat_name, "_cluster.jpeg")))
}

## posterior predictive densities
top_n_k <- arrange(k_post, desc(prob)) %>%
  pull(k) %>% .[1:configs$num_k_plot] %>% as.numeric()
pp_plot <- plot_posterior_densities(y, out, k=top_n_k)
pp_plot

if (configs$save){
  saveRDS(pp_plot,
          file.path("results", paste0(dat_name, "_pp.RDS")))
  ggsave(width=5, height=2, 
         file.path("figures", paste0(dat_name, "_pp.jpeg")))
}

## posterior means
post_means_by_k <-
  tibble(k=unique(out$k)) %>%
  mutate(post_means = map(k, post_means_fixed_k, out=out)) %>%
  unnest(post_means)

k_mode <- arrange(k_post, desc(prob)) %>%
  pull(k) %>% .[1]
post_means_by_k %>% filter(k==k_mode)

## acceptance rates and empty components
accept_stats <- get_accept_stats(out, agg=T) 
accept_stats                                  # acceptance rates

z_count <- apply(out$z, 1, function(x){length(unique(x))})
mean(out$k-z_count)                           # empty components

if (configs$save){
  mix_stats <- tibble(accept_e = accept_stats$p[1], 
                      accept_f = accept_stats$p[2], 
                      num_empty = mean(out$k-z_count))
  print(mix_stats)
  saveRDS(mix_stats,
          file.path("results", paste0(dat_name, "_stats.RDS")))
}
