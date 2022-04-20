library(readr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)

source(here::here('R/Fun/data_pipelines.R'))
source(here::here('R/Fun/plot/plot_styles.R'))

rt_clouds <- function(cloud_data){
  cloud_data %>%
    ggplot(aes(x=t, y=rt_sample, color=model)) +
    geom_point(size=0.5, shape=1) +
    geom_hline(yintercept=180, color="#333333") +
    model_colors + 
    facet_wrap(vars(model), ncol=2) +
    xlab("Trial") +
    ylab("Response Time (msec)") +
    style +
    guides(color=FALSE)
}

rt_cloud_marginals <- function(cloud_data){
  data_cloud <- cloud_data %>% subset(model=="Data") %>% subset(rt_sample > 180)
  cloud_data %>%
    subset(model != "Data") %>%
    rbind(data_cloud) %>%
    ggplot(aes(x=rt_sample, color=model, fill=model)) +
    geom_density(alpha=0.6) +
    model_fills +
    model_colors +
    #facet_grid(rows=vars(model)) + 
    xlab("Response Time (msec") +
    ylab("Density") +
    style + 
    theme(legend.position="bottom") + 
    guides(fill=guide_legend(title="Model"),
           color=FALSE) +
    coord_cartesian(expand = FALSE, xlim=c(120,500), ylim=c(0,0.012))
}

rt_cloud_fig <- function(rt_cloud_data){
  plot_cloud <- rt_clouds(rt_cloud_data)
  plot_marginals <- rt_cloud_marginals(rt_cloud_data)
  
  plot_grid(plot_cloud, plot_marginals, rel_heights=c(2,1), nrow=2, ncol=1)
}

#rt_cloud_data <- read_csv(here::here("/Python/Output/832a60f/rt_clouds.csv")) %>% relabel_and_order_models()
#print(rt_cloud_fig(rt_cloud_data))
