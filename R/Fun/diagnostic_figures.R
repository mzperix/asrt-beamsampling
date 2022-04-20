library(ggplot2)
library(magrittr)
library(dplyr)
library(viridis)

colors <- c("darkblue"= "#0000A1",
            "mediumblue" = "#1F6ED4",
            "lightblue" = "#39BAE8")

base_style <- theme_minimal()

mean_response_times_style <- base_style
parameter_diagnostic_plot_style <- base_style

figure_chain_diagnostics <- function(data){
  p <- data %>%
    ggplot(aes(y=ini, color=eff_sample_size)) +
    geom_point(aes(x=mean)) +
    geom_errorbarh(aes(xmin=q05, xmax=q95)) +
    facet_grid(epoch_start ~ participant_short, scales="free") +
    parameter_diagnostic_plot_style +
    theme(axis.text.x = element_text(angle = 90)) + 
    guides(color="none") +
    xlab(data$parameter) +
    ylab("Chain name") #+
    scale_color_viridis_d()
  return(p)
}

figure_residuals <- function(data){
  plot_residuals <- data %>%
    subset(filters==TRUE) %>%
    subset(correct_response==1) %>%
    relabel_and_order_models() %>%
    group_by(participant_test, model, e_test) %>%
    mutate(z_score = (r-mean(r)) / sd(r)) %>%
    ggplot(aes(color=model)) +
    #stat_bin(aes(y=..density..), bins=20, color="white", size=0.3) +
    geom_abline(slope=1, intercept=0, color="#747474", linetype="dashed") +
    stat_qq(aes(sample=z_score), geom="line", size=0.8) +
    #stat_density(adjust=0.5) +
    facet_rep_wrap("participant_test", ncol=5) +
    style +
    model_colors +
    coord_fixed(ratio=1, xlim=c(-2,2), ylim=c(-2,2)) + 
    ggtitle(paste('Residuals (r)','Train:',unique(data$e_train),'Test:',unique(data$e_test))) +
    theme(legend.position = "bottom"
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank()
          ) +
    guides(color=guide_legend(title="Model"))
  return(plot_residuals)
}