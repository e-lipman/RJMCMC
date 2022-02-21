# functions to interprest results

make_cluster_plot <- function(out){
  z <- out$z
  
  ij <- cross_df(list(i=1:ncol(z), j=1:ncol(z))) %>% 
    filter(i<j) %>%
    mutate(n_ij = map2(i,j, ~sum(z[,.x]==z[,.y]))) %>% 
    unnest(n_ij) %>%
    mutate(p = n_ij/nrow(z)) %>%
    select(-n_ij)
  ji <- ij %>% set_names(c("j","i","p"))
  ii <- tibble(i=1:ncol(z)) %>%
    mutate(j=i, p=1)
  
  bind_rows(ij, ji, ii) %>%
    ggplot(aes(x=i, y=j, col=p)) +
    geom_tile() +
    theme_bw()
}

# posterior predictive density
post_pred_dist_one <- function(i, w, mu, sig2){
  z <- sample(length(w[[i]]), 1, prob = w[[i]])  
  rnorm(1, mu[[i]][z], sig2[[i]][z])
}

post_pred_dist <- function(out){
  map_dbl(1:length(out$w), post_pred_dist_one,
          w=out$w, mu=out$mu, sig2=out$sig2)  
}

make_cluster_plot(out$z)

# posterior for fixed k
post_means_fixed_k <- function(out, k=3){
  w <- abind(out$w, along=0) %>% colMeans()
  mu <- abind(out$mu, along=0) %>% colMeans()
  sig2 <- abind(out$sig2, along=0) %>% colMeans()
  
  tibble(j = 1:k,
         w=w, mu=mu, sog2=sig2)
}



