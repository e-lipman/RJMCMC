# functions to interpret results

# plot clusters
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
    ggplot(aes(x=i, y=j, fill=p)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme_bw()
}

# posterior predictive density
normal_mix_dens <- function(w, mu, sig2, x_grid){
  map(1:length(w), 
      ~( w[[.x]]*dnorm(x_grid, mu[[.x]], sqrt(sig2[[.x]])) )) %>%
    abind(along=0) %>% colSums()
}

post_pred_density <- function(out, x_grid){
  dens_per_iter <- map(1:length(out$w),
                       ~normal_mix_dens(out$w[[.x]], 
                                        out$mu[[.x]], out$sig2[[.x]],
                                        x_grid))
  dens_per_iter %>% abind(along=0) %>% 
    colMeans()
}


# posterior for fixed k
post_means_fixed_k <- function(out, k=3){
  w <- abind(out$w, along=0) %>% colMeans()
  mu <- abind(out$mu, along=0) %>% colMeans()
  sig2 <- abind(out$sig2, along=0) %>% colMeans()
  
  tibble(j = 1:k,
         w=w, mu=mu, sig2=sig2)
}



