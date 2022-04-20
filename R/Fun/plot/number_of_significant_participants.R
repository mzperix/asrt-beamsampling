source(here::here('R/Fun/data_pipelines.R'))
source(here::here('R/Fun/plot/plot_styles.R'))

io_vs_markov_number_of_significant_participants <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(., dataset = "elarasztas") %>%
    relabel_and_order_models(.) %>%
    filter_and_average_runs(., 100, apply_filters=TRUE) %>%
    group_by(participant_test, e_train, session_test, model) %>%
    summarize(performance = cor(rt, rt_predicted)) %>%
    subset(model %in% c("Markov","Ideal Observer")) %>%
    ungroup() %>%
    mutate(session_test = as.integer(session_test)) %>%
    subset(session_test <= 8) %>%
    spread(model, performance) %>%
    mutate(better = ifelse(`Ideal Observer`>Markov,"Ideal Observer", "Markov")) %>%
    group_by(session_test, better) %>%
    summarize(n_significant = n()) %>%
    drop_na()
  
  df %>%
    drop_na() %>%
    ggplot(aes(x=better, y=n_significant, color=better, fill=better)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=n_significant), size=3, position=position_dodge(width=0.9), vjust=1.25, color="white") +
    facet_grid(.~session_test) +
    xlab("Day") +
    ylab(expression("Number of participants with larger R"^2)) +
    guides(fill=guide_legend(title="Model"), color=FALSE) +
    style + 
    model_fills +
    model_colors +
    #ylim(0,25) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank())
}


ihmm_markov_trigram_number_of_significant_participants <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(., dataset = "elarasztas") %>%
    relabel_and_order_models() %>%
    filter_and_average_runs(., 100, apply_filters=TRUE) %>%
    group_by(participant_test, e_train, session_test, model) %>%
    summarize(performance = cor(rt, rt_predicted),
              significance = cor.test(rt, rt_predicted, method="pearson")[['p.value']])
  df %>%
    mutate(significant = significance < 0.01) %>%
    group_by(session_test, model, significant) %>%
    summarize(n_significant = n()) %>%
    ggplot(aes(x=significant, y=n_significant, fill=model, group=model, color=model)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=n_significant), position=position_dodge(width=0.9), vjust=1.25, color="white") +
    facet_grid(model ~ session_test) +
    style +
    model_fills +
    model_colors
}