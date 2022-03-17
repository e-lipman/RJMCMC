library(tidyverse)
library(yaml)
library(pander)

configs <- read_yaml("config.yaml")

# summarise results across data sets

dat_names <- c("Enzyme","Acidity","Galaxy")

k_posterior_plot <- function(){
   kpost <- tibble(label=dat_names) %>% 
    mutate(kpost = map(label,
                       ~readRDS(file.path("results",
                                          paste0(tolower(.x), 
                                                 "_kpost.RDS"))))) %>%
    unnest(kpost) %>%
    mutate(k=as.numeric(k),
           prob=as.numeric(prob),
           label = fct_relevel(as.factor(label), 
                               dat_names)) %>%  
     filter(prob>=.001) %>%
    
    # color top n bars
    arrange(label,desc(prob)) %>%
    group_by(label) %>%
    mutate(top=( (1:n())<=configs$num_k_plot)) %>%
    ungroup()
     
   
   kpost %>% 
     ggplot(aes(x=k, y=prob, fill=top)) +
     geom_col(position=position_dodge(1),
              show.legend = FALSE) +
     facet_grid(.~label) +
     scale_fill_manual(values=c("#56B1F7","#132B43")) +
     ylab("") +
     theme_bw()
}

mixing_table <- function(){
  tibble(label=dat_names) %>% 
    mutate(mix = map(label, 
                     ~readRDS(file.path("results", 
                                        paste0(tolower(.x),  
                                               "_stats.RDS"))))) %>%
    unnest(mix) %>%
    mutate_if(is.numeric, ~round(100*.))
}

# results summaries

k_posterior_plot()
ggsave(width=5, height=2, 
       file.path("figures", "combined_kpost.jpeg"))

mixing_table() %>%
  pandoc.table()
