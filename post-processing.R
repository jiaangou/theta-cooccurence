#Post-processing 
library(dplyr)
library(tidyr)
library(ggplot2)


#load simulations
sim <- readRDS("simulation_runs.rds")


#data strucutre: netsted list
#first level: optima distribution (uniform,normal,skew; length = 3) 
sim%>%length
optima <- c('uniform','normal', 'skew')

#second level: parallel simulations (4 simulations)
sim[[1]]%>%length
iterations <- 1:4


#third level: `mcomsimr` simulation outputs (dynamics.df, landscape, env.df, env_traits, disp_mat, interaction)
sim[[1]][[1]]%>%names

#Combine each level into a single dataframe ------
dynamics.df.all <- list()

for (i in 1:length(optima)){
  
  dynamics.df.iterations <- list()
  for (j in iterations){
    #extract each iteration, label, and put in a list
    dynamics.df.iterations[[j]] <- sim[[i]][[j]]$dynamics.df%>%
      mutate(optima = optima[i])%>%
      mutate(iteration = j)
    
  }
  #combine the labeled iteration dfs into a single df by binding rows
  dynamics.df.all[[i]] <- dynamics.df.iterations%>% 
    bind_rows()

  #name the list elements by optima distribution
  names(dynamics.df.all) <- optima[1:i]
  
}

#COMBINE
dynamics.df.all <- dynamics.df.all%>%
  bind_rows(.id = 'optima_distribution')

