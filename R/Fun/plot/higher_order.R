source(here::here('R/Fun/data_pipelines.R'))

fig_higher_order <- function(data){
  ho_stats <- compute_higher_order_stats(data)
  il_mean <- ho_stats[["il_mean"]]
  correlations <- ho_stats[["correlations"]] %>% ungroup() %>% relabel_and_order_models()
  
  fig_ho_ihmm <- il_mean %>%
    mutate(rt_significant=(sign(rt_diff_mean-rt_diff_se*2)==sign(rt_diff_mean+rt_diff_se*2)))%>%
    ungroup() %>%
    relabel_and_order_models() %>%
    ggplot(aes(x=rt_diff_mean, y=rt_predicted_diff_mean)) +
    geom_hline(yintercept=0.0, size=0.7, color=pantone_2014[['grey']])+
    geom_vline(xintercept=0.0, size=0.7, color=pantone_2014[['grey']]) +
    geom_abline(intercept = 0.0, slope = 1.0, color=pantone_2014[['grey']]) +
    style +
    scale_x_continuous(minor_breaks = c(-10,0,10), breaks = c(-20,-10,0,10,20)) +
    scale_y_continuous(minor_breaks = c(-10,0,10,20), breaks = c(-10,0,10,20)) +
    coord_fixed(xlim=c(-25,25), ylim=c(-25,25)) +
    #coord_fixed() +
    xlab('Higher-Order Statistical Learning \n (Response Time diff. in msec)') +
    ylab('Higher-Order Statistical Learning \n (Model prediction diff. in msec)') +
    theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = pantone_2014[['light_grey']]),
          panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                          colour = pantone_2014[["light_grey"]]),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=rt_predicted_diff_mean-rt_predicted_diff_se*2,
                      ymax=rt_predicted_diff_mean+rt_predicted_diff_se*2), color=pantone_2014[['light_grey']]) +
    geom_errorbarh(aes(xmin=rt_diff_mean-rt_diff_se*2, xmax=rt_diff_mean+rt_diff_se*2), color=pantone_2014[['light_grey']]) +
    geom_point(size=1.3, aes(color=rt_significant)) +
    scale_color_manual(values = c(pantone_2014[['grey']],pantone_2014[['orange']]))+
    guides(colour = FALSE) +
    geom_text(data=correlations,
              mapping = aes(x = -Inf, y = Inf,
                            label=paste('r=', format(correlations$correlation,digits=2))),
              hjust   = -0.1,
              vjust   = 1) +
    # geom_text(data=correlations,
    #           mapping = aes(x = -Inf, y = Inf,
    #                         label=paste('p=', format(correlations$p_value,digits=5))),
    #           hjust   = -0.1,
    #           vjust   = 2.2) +
    facet_grid(session_test~model)

  return(fig_ho_ihmm)
}
