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
simulation <- 1:4


#third level: `mcomsimr` simulation outputs (dynamics.df, landscape, env.df, env_traits, disp_mat, interaction)
sim[[1]][[1]]%>%names

#Combine each level into a single dataframe ------
dynamics.df.all <- list()

for (i in 1:length(optima)){
  
  dynamics.df.simulations <- list()
  for (j in simulation){
    #extract each iteration, label, and put in a list
    dynamics.df.simulations[[j]] <- sim[[i]][[j]]$dynamics.df%>%
      mutate(optima = optima[i])%>%
      mutate(simulation = j)
    
  }
  #combine the labeled iteration dfs into a single df by binding rows
  dynamics.df.all[[i]] <- dynamics.df.simulations%>% 
    bind_rows()

  #name the list elements by optima distribution
  names(dynamics.df.all) <- optima[1:i]
  
}

#COMBINE into a nested tibble grouped by optima distribution and simulation 
dynamics.df.all <- dynamics.df.all%>%
  bind_rows(.id = 'optima_distribution')%>%
  group_by(optima_distribution, simulation)%>% 
  nest()

#Calculate theta for each simulated metacommunity------
#Ia. convert dynamics.df to wide format (sample x species)
spe_compo_wide <- function(x){
  out <- x%>%
    mutate(species = paste0('sp',species))%>%
    pivot_wider(id_cols = c(patch, env), names_from = 'species', values_from = N, values_fill = 0) #fills NAs with 0
  return(out)
}

dynamics.df.all$data[[1]]%>%
  spe_compo_wide()


#Ib. extract spp (fundamental) parameters (optima, niche width, max_r)
spe_parameters <- function(x){
  out <- x%>%
    select(species, optima, env_niche_breadth, max_r)%>%
    distinct()
  return(out)
}

dynamics.df.all$data[[1]]%>%
  spe_parameters()

#II. Calculate theta
devtools::install_github ('zdealveindy/theta')
library(theta)





