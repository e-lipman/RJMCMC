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
    xlab("") + ylab("") +
    theme_bw()
}

# posterior predictive density
normal_mix_dens <- function(w, mu, sig2, x_grid){
  map(1:length(w), 
      ~( w[[.x]]*dnorm(x_grid, mu[[.x]], sqrt(sig2[[.x]])) )) %>%
    abind(along=0) %>% colSums()
}

post_pred_density <- function(out, x_grid, k=NA){
  if (!is.na(k)){
    k_idx <- which(out$k==k)
  } else {
    k_idx <- 1:length(out$k)
  }
  dens_per_iter <- map(k_idx,
                       ~normal_mix_dens(out$w[[.x]], 
                                        out$mu[[.x]], out$sig2[[.x]],
                                        x_grid))
  dens_per_iter %>% abind(along=0) %>% 
    colMeans()
}

plot_posterior_densities <- function(y, out, k, combined=T,
                                     x_range = "data"){
  k <- intersect(k, out$k)
  
  stopifnot(x_range %in% c("data","mu"))
  if (x_range=="data"){
    x_grid <- seq(min(y), max(y), length=100)
  } else {
    range_mu <- range(unlist(out$mu))
    x_grid <- seq(range_mu[1],range_mu[2],length=100)
  }
  
  dens_by_k <- tibble(k=k) %>%
    mutate(dens=map(k,
                    ~tibble(x=x_grid,
                            y=post_pred_density(out, x_grid, k=.x)))) %>%
    unnest(dens) %>%
    mutate(k=as.character(k))
  if (combined){
    dens_combined <- tibble(k="all",
           x=x_grid,
           y=post_pred_density(out, x_grid))
  }
  
  out <- dens_by_k %>%
    ggplot(aes(x=x, y=y, color=k)) +
    geom_histogram(data=tibble(y=y),
                   aes(x=y, y = ..density..), 
                   alpha=.4,
                   inherit.aes=F) +
    geom_line(linetype=ifelse(combined, "dashed","solid")) +
    xlab("") + ylab("") +
    theme_bw()
  
  if (combined){
    out <- out + 
      geom_line(data = dens_combined, 
                color="black", linetype="solid", size=.25)
  }
  return(out)
}


# posterior for fixed k
post_means_fixed_k <- function(out, k=3){
  k_idx <- which(out$k==k)
  w <- abind(out$w[k_idx], along=0) %>% colMeans()
  mu <- abind(out$mu[k_idx], along=0) %>% colMeans()
  sig2 <- abind(out$sig2[k_idx], along=0) %>% colMeans()
  
  tibble(j = 1:k,
         w=w, mu=mu, sig2=sig2)
}

# get acceptance stats
get_accept_stats <- function(out, agg = F){
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
  
  if (agg){
    accept_stats %>%
      mutate(group=rep(c("e","f"), each=2)) %>%
      group_by(group) %>%
      summarise(x=sum(x), n=sum(n), p=x/n)
  } else {
    accept_stats
  }
}



