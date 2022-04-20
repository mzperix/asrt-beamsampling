library(diagram)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(glue)
library(jsonlite)
library(grid)

colors <- c("#37324f",
  "#6697c3",
  "#ffd00d",
  "#fb8b00")

model_figure_colors <- c("arrow" = "#8E8E8E",
                         "dot"   = "#8E8E8E",
                         "stim1" = "#DF422A",
                         "stim2" = "#FDD262",
                         "stim3" = "#D3DDDC",
                         "stim4" = "#34A5DA")

test_model_figure <- function(){
  N <- 5
  
  Pi <- matrix(runif(N^2), ncol=N)
  Phi <- matrix(runif(4*N), ncol=4)
  Pi <- Pi / rowSums(Pi)
  Phi <- Phi / rowSums(Phi)
  draw_model(Pi,Phi)  
}

reorder_states <- function(Pi,Phi){
  n <- dim(Pi)[1]
  Pi2 <- Pi[1:n,1:n]
  diag(Pi2) <- -1
  s <- 1
  l <- c(s)
  O <- matrix(rep(0,n*n),ncol=n)
  O[1,1] <- 1
  for (i in 2:n){
    Pi2[,s] <- -1
    next_s <- which.max(Pi2[s,])
    O[next_s,i] <- 1
    l <- c(l,next_s)
    s <- next_s
  }
  return(list(t(O)%*%Pi%*%O,Phi[l,]))
}

bound_transitions <- function(Pi,include_state_bound,include_edge_bound){
  n <- dim(Pi)[1]
  if (n==1) return(matrix(c(1,0,0,0.01), ncol=2))
  Pi <- Pi[1:n,1:n]
  Pi <- Pi[colSums(Pi)>include_state_bound, 
           colSums(Pi)>include_state_bound]
  Pi <- Pi*(Pi>include_edge_bound)
  return(Pi / rowSums(Pi))
}

state_pos <- function(n,r){
  return(data.frame(x=r*cos((1:n-1)*2*pi/n), y=r*sin((1:n-1)*2*pi/n)))
}

draw_emissions <- function(Phi, pp){
  dist <- 0.03
  n <- dim(Phi)[1]
  pos <- state_pos(n,2.7)
  emissions <- cbind(pos, data.frame(Phi)) %>%
    gather(key="Y",value="p", -c(x,y)) %>%
    mutate(p = sqrt(p)/2)
  pp + 
    geom_rect(data=subset(emissions, Y=="X1"), aes(xmin=x+dist, xmax=x+p+dist, ymin=y+dist, ymax=y+p+dist), fill=model_figure_colors[["stim1"]]) +
    geom_rect(data=subset(emissions, Y=="X2"), aes(xmin=x-p-dist, xmax=x-dist, ymin=y+dist, ymax=y+p+dist), fill=model_figure_colors[["stim2"]]) +
    geom_rect(data=subset(emissions, Y=="X3"), aes(xmin=x-p-dist, xmax=x-dist, ymin=y-p-dist, ymax=y-dist), fill=model_figure_colors[["stim3"]]) +
    geom_rect(data=subset(emissions, Y=="X4"), aes(xmin=x+dist, xmax=x+p+dist, ymin=y-p-dist, ymax=y-dist), fill=model_figure_colors[["stim4"]]) %>%
    return()
}

draw_transitions <- function(Pi, p){
  n <- dim(Pi)[1]
  r <- 1.8
  pos <- state_pos(n,r)
  pos <- cbind(pos, data.frame(self_transition=diag(Pi)))
  p1 <- cbind(state_pos(n,r),data.frame(state=1:n))
  p2 <- merge(p1, expand.grid(1:n,1:n), by.x="state",by.y="Var1")
  p3 <- merge(p2, p1, by.x="Var2", by.y="state") %>%
    arrange(state,Var2) %>%
    cbind(data.frame(p=as.vector(t(Pi)))) %>%
    subset(state != Var2)
  w <- 0.23
  p3 <- p3 %>%
    mutate(x.y_arrow = w*x.x+(1-w)*x.y,
           y.y_arrow = w*y.x+(1-w)*y.y,
           s = p) %>%
    subset(s>0)
  p <- p + 
    geom_point(data=pos, aes(x=x, y=y, size=sqrt(pos$self_transition)), color=model_figure_colors[["dot"]]) +
    geom_curve(data=p3,
      aes(x=x.x, xend=x.y_arrow, y=y.x, yend=y.y_arrow, size=s/12), 
      curvature=0.0,
      lineend = "butt",
      arrow = arrow(length = unit(p3$s/17,"npc"), angle=23, type="closed", ends = "last"),
      color=model_figure_colors[["arrow"]]) +
    geom_curve(data=p3,
      aes(x=x.x, xend=x.y, y=y.x, yend=y.y, size=s/22),
      curvature=0.0,
      lineend = "butt",
      color=model_figure_colors[["arrow"]]) +
    guides(size=FALSE) +
    scale_size_continuous(range=c(unit(0,"npc"),unit(5,"npc")), limits=c(0,1))
  return(p)
}

draw_model <- function(Pi, Phi, include_state_bound, include_edge_bound){
  p <- ggplot() + 
    theme_classic() + 
    coord_fixed() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  if (dim(Pi)[1] == 1) return(ggplot())
  Pi <- bound_transitions(Pi, include_state_bound, include_edge_bound)
  #print(Pi)
  reordered <- reorder_states(Pi,Phi)
  Pi <- reordered[[1]]
  Phi <- reordered[[2]]
  #print(Pi)
  p <- draw_transitions(Pi, p)
  p <- draw_emissions(Phi, p)
  return(p)
}

figure_models <- function(data, b1,b2){
  participants = unique(data$participant_train)
  plots <- list()
  i = 0
  for (p in participants){
    p_data <- data %>%
      subset(participant_train == p) %>%
      group_by(participant_train, ini, model) %>%
      subset(iteration==max(iteration))
    Pi <- p_data$pi[1][[1]]
    Phi <- p_data$phi[1][[1]]
    i <- i+1
    plots[[i]] <- draw_model(Pi,Phi,b1,b2)
  }
  return(plot_grid(plotlist=plots))
}

sequence_plot <- function(sequence){
  colors <- c(
    model_figure_colors[["stim1"]],
    model_figure_colors[["stim2"]],
    model_figure_colors[["stim3"]],
    model_figure_colors[["stim4"]]
  )
  data.frame(x=factor(sequence+1, levels=c(1,2,3,4)), pos=c(1,2,3,4)) %>%
    ggplot(aes(x=pos, color=x)) +
    geom_point(y=0, size=10) +
    ylim(-1,1) +
    xlim(0,5) +
    ggtitle("Sequence") +
    guides(color=FALSE) +
    theme_void() +
    scale_color_manual(values=colors)
}

fingers_plot <- function(data){
  colors <- c(
    model_figure_colors[["stim1"]],
    model_figure_colors[["stim2"]],
    model_figure_colors[["stim3"]],
    model_figure_colors[["stim4"]]
  )
  data.frame(x=as.factor(c(1,2,3,4)), pos=c(1,2,4,5)) %>%
    ggplot(aes(x=pos, color=x)) +
    geom_point(y=0, size=10) +
    ylim(-1,1) +
    xlim(-0.5,6.5) +
    ggtitle("Fingers") +
    guides(color=FALSE) +
    theme_void() +
    scale_color_manual(values=colors)
}

fig_example_model <- function(){
  Phi <- matrix(c(0.1,0.0,0.1,0.8, 
                  0.3,0.5,0.2,0.0,
                  0.25,0.25,0.25,0.25,
                  0.2,0.8,0.0,0.0,
                  0.8,0.1,0.1,0.0,
                  0.0,0.2,0.8,0.0), 
                ncol=4, byrow=TRUE)
  Pi <-  matrix(c(0.1,0.7,0.2,0.0,0.0,0.0,
                  0.0,0.1,0.6,0.3,0.0,0.0,
                  0.2,0.2,0.2,0.4,0.0,0.0,
                  0.0,0.0,0.0,0.1,0.4,0.2,
                  0.0,0.0,0.0,0.0,0.5,0.5,
                  0.6,0.1,0.2,0.0,0.0,0.1),
                ncol=6, byrow=TRUE)
  df <- data.frame(Pi=I(list(Pi)), Phi=I(list(Phi)),
                   participant_train=c(1),
                   ini=c('manual'),
                   model=c('iHMM'),
                   iteration=c(1))
  return(figure_models(df,0.0,0.0))
}

#df <- fromJSON("/Projects/asrt-beamsampling/Python/Output/4e9c1cf/d8cf5a2_MODEL_SAMPLES_elarasztas_model_samples_markov.json")
#plots <- df %>%
#  subset(e_train=='186_195') %>%
#  subset(participant_train=='elarasztas_119') %>%
#  figure_models(., 0.2)
#plots

