library(readr)
library(magrittr)
library(ggplot2)

source("R/Fun/plot/learning_curves.R")
source("R/Fun/result_figures.R")
source("R/Fun/data_pipelines.R")
source("R/Fun/plot/plot_styles.R")
source("R/Fun/kl_divergence.R")
source("R/Fun/plot/number_of_significant_participants.R")
source("R/Fun/prediction_spaces.R")
source("R/Fun/diagnostic_figures.R")
source("R/Fun/plot/fingerprint.R")
source("R/Fun/plot/prior_predictive_checks.R")

####################
###  INPUT DATA  ###
####################

input_dir <- Sys.getenv("DATA_DIR")
filename_rt_prediction_data <- "866066d_LEARNING_CURVES_elarasztas_pre_submission"
filename_probability_prediction_data <- "866066d_PREDICTED_PROBABILITIES_elarasztas_pre_submission"
filename_across_sequences_prediction_data <- "866066d_ACROSS_SEQUENCES_elarasztas_pre_submission"
filename_probability_prediction_cross_entropy_data <- "31f7b92_PREDICTED_PROBABILITIES_elarasztas_ploscb_rev_pred_probs"
filename_residual_data <- "839179e_RESIDUALS_elarasztas"
filename_rt_cloud_data <- "rt_clouds"
filename_fingerprint_data <- "866066d_FINGERPRINT3_elarasztas_pre_submission"

load_data <- function(filename){
    return(read_csv(paste0(input_dir, "/", filename, ".csv")))
}

rt_prediction_data <- load_data(filename_rt_prediction_data) %>%
    add_session_info(., dataset = "elarasztas")

ihmm_markov_triplet_data <- rt_prediction_data %>%
      mutate(session_train = e_train) %>%
      mutate(session = session_test)

probability_prediction_data <- load_data(filename_probability_prediction_data)
across_sequences_data <- load_data(filename_across_sequences_prediction_data)
probability_prediction_cross_entropy_data <- load_data(filename_probability_prediction_cross_entropy_data)
residual_data <- load_data(filename_residual_data)
rt_cloud_data <- load_data(filename_rt_cloud_data)
fingerprint_data <- load_data(filename_fingerprint_data)
artificial_asrt_params <- read_csv("artificial_asrt_params")
######################
### CREATE FIGURES ###
######################

final_fig2c <- rt_prediction_data %>%
    relabel_and_order_models() %>%
    subset(participant_test == 119) %>%
    subset(e_test == "176_185") %>%
    subset(model == "CT") %>%
    fig_individual_trials() +
    style

final_fig4 <- rt_prediction_data %>%
    filter(model %in% c("iHMM", "GroundTruth")) %>%
    #relabel_and_order_models() %>%
    #filter(model %in% c("CT", "Ideal Observer")) %>%
    fig_model_comparison()

final_fig5 <- probability_prediction_data %>%
    subset(model %in% c("iHMM", "GroundTruth")) %>%
    fig_error_e(., participants = c(102, 110))

final_fig5e <- across_sequences_data %>%
    fig_across_sequences()

final_fig6a <- fig3a(ihmm_markov_triplet_data)
final_fig6e <- fig_learning_curves_individual(rt_prediction_data)
final_fig6g <- fig6gh(probability_prediction_cross_entropy_data, "26_35", c(102,110,119))
final_fig6h <- fig6gh(probability_prediction_cross_entropy_data, "176_185", c(102,110,119))

final_fig7a <- io_vs_markov_number_of_significant_participants(rt_prediction_data)
final_fig7bde <- fig4b(ihmm_markov_triplet_data, rt_prediction_data)
final_fig7c <- fig_predictability(rt_prediction_data)

### Supplementary ###
final_figS4 <- figure_residuals(residual_data)
final_figS5 <- fig2b(rt_prediction_data, 119, rt_cloud_data)
final_figS6 <- fingerprint_data %>%
      subset(e_train == "186_195") %>%
      mutate(e_train = paste0("e_train_", e_train)) %>%
      fig_fingerprint_predicted_rt(min_rt_count=5)

final_figS7 <- probability_prediction_data %>%
    filter(model %in% c("iHMM", "GroundTruth")) %>%
    fig_error_roc_all(ncol=5)

final_figS8 <- across_sequences_data %>% fig_across_sequences_individual()
final_figS10 <- ihmm_markov_triplet_data %>% fig_ihmm_markov_io()
final_figS11 <- fig_prior_predictive_checks(artificial_asrt_params)

####################
### SAVE FIGURES ###
####################

output_path <- "Figures/figlist"
fig_sizes <- c(
    fig2c = c(width = 3, height = 5),
    fig4 = c(width = 7, height = 3),
    fig5 = c(width = 5.5, height = 10),
    fig5e = c(width = 8.5, height = 5),
    fig6a = c(width = 7.5, height = 3.5),
    fig6e = c(width = 10.5, height = 8.5),
    fig6g = c(width = 3.5, height = 3.5),
    fig6h = c(width = 3.5, height = 3.5),
    fig7a = c(width = 5.5, height = 3.0),
    fig7bde = c(width = 5.5, height = 8.5),
    fig7c = c(width = 5.5, height = 4.0),
    figS4 = c(width = 8.5, height = 8.5),
    figS5 = c(width = 11.0, height = 6.0),
    figS6 = c(width = 7.0, height = 14.0),
    figS7 = c(width = 10.0, height = 10.0),
    figS8 = c(width = 8.5, height = 8.5),
    figS9 = c(width = 10.5, height = 8.5),
    figS10 = c(width = 9.5, height = 6.5),
    figS11 = c(width = 12.5, height = 12.5)
)

export_fig <- function(input_name, figure_name, figure) {
    ggsave(
        paste0(figure_name, "_", input_name, ".pdf"),
        figure,
        device = "pdf",
        path = output_path,
        width = getElement(fig_sizes, paste(figure_name, "width", sep = ".")),
        height = getElement(fig_sizes, paste(figure_name, "height", sep = "."))
    )
}

export_fig(filename_rt_prediction_data, "fig2c", final_fig2c)
export_fig(filename_rt_prediction_data, "fig4", final_fig4)
export_fig(filename_probability_prediction_data, "fig5", final_fig5)
export_fig(filename_across_sequences_prediction_data, "fig5e", final_fig5e)
export_fig(filename_rt_prediction_data, "fig6a", final_fig6a)
export_fig(filename_rt_prediction_data, "fig6e", final_fig6e)
export_fig(filename_probability_prediction_cross_entropy_data, "fig6g", final_fig6g)
export_fig(filename_probability_prediction_cross_entropy_data, "fig6h", final_fig6h)
export_fig(filename_rt_prediction_data, "fig7a", final_fig7a)
export_fig(filename_rt_prediction_data, "fig7bde", final_fig7bde)
export_fig(filename_rt_prediction_data, "fig7c", final_fig7c)

export_fig(filename_residual_data, "figS4", final_figS4)
export_fig(filename_rt_prediction_data, "figS5", final_figS5)
export_fig(filename_fingerprint_data, "figS6", final_figS6)
export_fig(filename_probability_prediction_data, "figS7", final_figS7)
export_fig(filename_across_sequences_data, "figS8", final_figS8)
export_fig(filename_rt_prediction_data, "figS9", final_fig6e) # Same as fig6e
export_fig(filename_rt_prediction_data, "figS10", final_figS10)
export_fig("artificial_asrt_params", "figS11", final_figS11)
