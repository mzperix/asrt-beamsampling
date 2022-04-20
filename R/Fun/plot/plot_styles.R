library(ggplot2)

colors <- c("darkblue"= "#0000A1",
            "mediumblue" = "#1F6ED4",
            "lightblue" = "#39BAE8")

base_style <- theme_minimal()

mean_response_times_style <- base_style
parameter_diagnostic_plot_style <- base_style

pantone_2014 = c('cayenne' = '#E66665',
                 'light_cayenne' = '#ED9393',
                 'dark_cayenne' = '#A14746',
                 'dark_blue' = '#076BB6',
                 'darker_blue' = '#044A7F',
                 'sand' = '#CDB48C',
                 'light_sand' = '#DEC59D',
                 'dark_sand' = '#8C734B',
                 'green' = '#9ECEB4',
                 'light_blue' = '#5197CB',
                 'yellow' = '#FFD602',
                 'purple' = '#AE70AF',
                 'violet' = '#9295CA',
                 'grey' = '#A9B2B1',
                 'orange' = '#F47D43',
                 'light_grey' = '#D9E2E1',
                 'dark_grey' = '#343434')

entropy_colors <- c('blue' = '#4CA3DD',
                    'dark_blue' = '#2d6184',
                    'green' = '#66CDAA',
                    'dark_green' = '#478F76')

model_colors_list <- c("CT"=pantone_2014[['cayenne']],
                       "CT Day 9"=pantone_2014[['light_cayenne']],
                       "CT Day 10"=pantone_2014[['dark_cayenne']],
                       "Markov"=pantone_2014[['dark_blue']],
                       "Markov Day 9"=pantone_2014[['light_blue']],
                       "Markov Day 10"=pantone_2014[['darker_blue']],
                       "Trigram"=pantone_2014[['sand']],
                       "Trigram Day 8"=pantone_2014[['dark_sand']],
                       "Trigram Day 9"=pantone_2014[['light_sand']],
                       "Ideal Observer"=entropy_colors[['green']],
                       "Expert model"=entropy_colors[['green']],
                       "Data"=pantone_2014[['dark_grey']])

model_colors <- scale_color_manual(c("CT", "CT Day 9", "CT Day 10",
                                     "Markov", "Markov Day 9", "Markov Day 10",
                                     "Trigram", "Trigram Day 9", "Trigram Day 8",
                                     "Ideal Observer", "Expert model","Data"), 
                                   values=model_colors_list,
                                   drop=TRUE,
                                   limits=force) 
model_fills <-  scale_fill_manual(c("CT", "CT Day 9", "CT Day 10",
                                    "Markov", "Markov Day 9", "Markov Day 10",
                                    "Trigram", "Trigram Day 9", "Trigram Day 8",
                                    "Ideal Observer", "Expert model","Data"), 
                                  values=model_colors_list,
                                  drop=TRUE,
                                  limits=force)

model_across_colors <- scale_color_manual(
  values=c(pantone_2014[['cayenne']],
           pantone_2014[['orange']],
           pantone_2014[['light_blue']],
           pantone_2014[['dark_blue']],
           pantone_2014[['light_sand']],
           pantone_2014[['dark_sand']]))

style <- theme_classic() +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = pantone_2014[["light_grey"]]), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        text = element_text(family='Helvetica', color = "#343434"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),
        strip.text = element_text(size=10),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.6, color = "#747474")
  )