fig_ihmm_markov_triplet <- function(data){
  p <- data %>%
    add_session_info("elarasztas") %>%
    compute_ihmm_markov_triplet_diff() %>%
    ungroup() %>%
    mutate(session_test=ordered(
      session_test,
      levels=c(1,2,3,4,5,6,7,8),
      labels=paste('Day',c(1,2,3,4,5,6,7,8))
    )) %>%
    ggplot(aes(x=Triplet, y=iHMM-Markov, text=paste("participant:",participant_test))) +
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_point(color=pantone_2014[['dark_blue']], size=1) +
    facet_wrap(vars(session_test), ncol=4) +
    coord_fixed() +
    xlab(bquote('Trigram performance'~(R^2))) +
    ylab(bquote('CT perf.- Markov perf.'~(R^2))) +
    style
  return(p)
}

fig_ihmm_markov_io <- function(data){
  p <- data %>%
    add_session_info("elarasztas") %>%
    compute_ihmm_markov_io_diff() %>%
    ungroup() %>%
    mutate(session_test=ordered(
      session_test,
      levels=c(1,2,3,4,5,6,7,8),
      labels=paste('Day',c(1,2,3,4,5,6,7,8))
    )) %>%
    ggplot(aes(x=GroundTruth, y=iHMM-Markov, text=paste("participant:",participant_test))) +
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_point(color=pantone_2014[['dark_blue']], size=1) +
    facet_wrap(vars(session_test), ncol=4) +
    coord_fixed() +
    xlab(bquote('Ideal Observer performance'~(R^2))) +
    ylab(bquote('CT perf.- Markov perf.'~(R^2))) +
    style
  return(p)
}

fig_ihmm_markov_triplet_hist <- function(data){
  d <- compute_ihmm_markov_triplet_diff(data)

  pvalue1 <- t.test(subset(d, session_test %in% c(6,7,8))$diff,
                    alternative="greater")[["p.value"]]

  p1 <- d %>%
    subset(session_test %in% c(6,7,8)) %>%
    ggplot(aes(x=diff)) +
    annotate('text', x=Inf, y=Inf, hjust=1, vjust=1,
             label=paste("p=",format(pvalue1,digits=4),sep=""),
             size=6) +
    geom_histogram(bins=12, fill=pantone_2014[["violet"]]) +
    xlab('iHMM-Markov-Triplet') +
    ylab('Count')+
    ggtitle('Day 6-8')+
    style

  pvalue2 <- t.test(subset(d, session_test %in% c(8))$diff,
                    alternative="greater")[["p.value"]]
  p2 <- d %>%
    subset(session_test %in% c(8)) %>%
    ggplot(aes(x=diff)) +
    geom_histogram(bins=12, fill=pantone_2014[["violet"]]) +
    annotate('text', x=Inf, y=Inf, hjust=1, vjust=1,
             label=paste("p=",format(pvalue2,digits=4),sep=""),
             size=6) +
    xlab('iHMM-Markov-Triplet') +
    ylab('Count')+
    ggtitle('Day 8')+
    style
  return(plot_grid(p1,p2))
}
