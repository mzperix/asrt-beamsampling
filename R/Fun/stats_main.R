library(readr)
source(here::here('R/Fun/stats.R'))

## LEARNING CURVES
rt_prediction_data <- read_csv(here::here('Python/Output/832a60f/866066d_LEARNING_CURVES_elarasztas_pre_submission.csv'))
cat(learning_curve_stats(rt_prediction_data, 119))
cat(higher_order_stats(rt_prediction_data))

## ERROR PREDICTION
error_prediction_data <- read_csv(here::here('Python/Output/832a60f/866066d_PREDICTED_PROBABILITIES_elarasztas_pre_submission.csv'))
cat(error_prediction_tests(error_prediction_data))

## ACROSS PARTICIPANTS
across_participants_data <- read_csv(here::here('Python/Output/832a60f/866066d_ACROSS_PARTICIPANTS_elarasztas_pre_submission.csv'))
cat(across_participant_tests(across_participants_data))

## ACROSS SEQUENCES
across_sequences_data <- read_csv(here::here('Python/Output/832a60f/866066d_ACROSS_SEQUENCES_elarasztas_pre_submission.csv'))
cat(across_sequence_comparisons(across_sequences_data))
