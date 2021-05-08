#Post-processing 
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)


#load simulations

sim <- here("data","simulation_runs.rds")%>%
  readRDS()


#data strucutre: netsted list
#first level: optima distribution (uniform,normal,skew; length = 3) 
sim%>%length
optima_dist <- c('uniform','normal', 'skew')

#second level: parallel simulations (4 simulations)
sim[[1]]%>%length
simulation <- 1:4


#third level: `mcomsimr` simulation outputs (dynamics.df, landscape, env.df, env_traits, disp_mat, interaction)
sim[[1]][[1]]%>%names

#Combine each level into a single dataframe ------
dynamics.df.all <- list()

for (i in 1:length(optima_dist)){
  
  dynamics.df.simulations <- list()
  for (j in simulation){
    #extract each iteration, label, and put in a list
    dynamics.df.simulations[[j]] <- sim[[i]][[j]]$dynamics.df%>%
      mutate(optima_dist = optima_dist[i])%>%
      mutate(simulation = j)
    
  }
  #combine the labeled iteration dfs into a single df by binding rows
  dynamics.df.all[[i]] <- dynamics.df.simulations%>% 
    bind_rows()

  #name the list elements by optima distribution
  names(dynamics.df.all) <- optima_dist[1:i]
  
}

#COMBINE into a nested tibble grouped by optima distribution and simulation 
dynamics.df.all <- dynamics.df.all%>%
  bind_rows(.id = 'optima_distribution')%>%
  group_by(optima_distribution, simulation)%>% 
  nest()


#save simulation data
saveRDS(object = dynamics.df.all, file = 'simulation_data.rds')

#Calculate theta for each simulated metacommunity------
#Ia. convert dynamics.df to wide format (sample x species)
spe_compo_wide <- function(x){
  out <- x%>%
    mutate(species = paste0('sp',species))%>%
    pivot_wider(id_cols = c(patch, env), names_from = 'species', values_from = N, values_fill = 0) #fills NAs with 0
  return(out)
}

spec.comp <- dynamics.df.all$data%>%
  lapply(function(x)x%>%
           spe_compo_wide()%>%
           select(-patch,-env))


#Ib. extract spp (fundamental) parameters (optima, niche width, max_r)
spe_parameters <- function(x){
  out <- x%>%
    select(species, optima, env_niche_breadth, max_r)%>%
    distinct()%>%
    mutate(species = paste0('sp',species))
  return(out)
}

fundamental <- dynamics.df.all$data%>%
  lapply(FUN = function(x)spe_parameters(x))%>%
  tibble(spp_parameters = .)%>%
  bind_cols(dynamics.df.all[,1:2],. )


#Spp paramaters set for simulations
saveRDS(fundamental, file = 'Spp_parameters.rds')


#II. Calculate theta
#devtools::install_github ('zdealveindy/genspe')
library(genspe)

theta <- dynamics.df.all$data[[1]]%>%
  spe_compo_wide()%>%
  select(-patch,-env)%>%
  calculate.theta()%>%
  rename(`species` = 'sci.name')

#BD metrics table
metrics <- read.delim ('theta metrics for evaluation.txt', row.names = 1, stringsAsFactors = F)


#Parallelize
library(parallel)


cl <- makeCluster (6)
theta.out <- c()
clusterExport (cl, c ('theta.out', 'spec.comp','metrics'))
theta.output <- parLapply(cl, spec.comp, function(x){
  
  
  #loop through metrics and apply each to theta function
  for (j in rownames(metrics)){
    
    theta.out[[j]] <- genspe::calculate.theta(x,
                                              method = metrics[j,"method"],
                                              beta.div.method = metrics[j,"beta.div.method"],
                                              beta.div.sqrt.D = metrics[j, "beta.div.sqrt.D"],
                                              force.subsample = metrics[j, "force.subsample"],
                                              reps =100)}
 
  #bind rows with new column as metric name
  out <- dplyr::bind_rows(theta.out, .id = 'metric')
  
  return(out)
  })

stopCluster (cl)



#calculate correlation between theta and niche width (fundamental)


#combine fundmanetal and theta tables
theta_fundamental <- mapply(x = theta.output, y = fundamental, FUN = function(x,y)x%>%
         rename('species' = `sci.name`)%>%
         left_join(y, by = 'species'), SIMPLIFY = FALSE)


#combined theta-fundamental table into master nested tibble -------
#create a theta tibble with the same optima disitrbution and simulation variables as master tibble
theta_tibble <- tibble(optima_distribution = dynamics.df.all$optima_distribution,
       simulation = dynamics.df.all$simulation, 
       theta_metrics = theta_fundamental)
#join the 2 nested tibbles
master_tibble <- dynamics.df.all%>%
  left_join(theta_tibble)
  
#Calculate correlation by optima distribution, simulation, metric
cor_table <- master_tibble%>%
  unnest(theta_metrics)%>%
  group_by(optima_distribution, simulation, metric)%>%
  summarise(correlation = cor(theta, env_niche_breadth, method = 'spearman'))
  
#Visualize data and summary statistics ----------
library(RColorBrewer)
cols <- brewer.pal(3, name = 'Set1')

cor_table%>%
  group_by(optima_distribution, metric)%>%
  summarise(mean_correlation = mean(correlation), sd_correlation = sd(correlation))%>%
  ggplot(aes(x = metric, y = mean_correlation, col = optima_distribution)) + 
  geom_bar(aes(fill = optima_distribution), position = position_dodge(), stat = 'identity')+
  scale_fill_manual(values = cols)+
  theme_bw()



################################################
#Saturated VS unsaturated communities ----------------
################################################

#Un/saturate communities by # species and by # individual
unsaturated_comm <- dynamics.df.all$data%>%
  lapply(FUN = function(x)saturation_sampling(dynamics.df = x, k = 0.5, individuals = 100))


un_saturated_tibble <- tibble(optima_distribution = dynamics.df.all$optima_distribution,
                             simulation = dynamics.df.all$simulation,
                             unsaturated = lapply(unsaturated_comm, function(x)x$unsaturated)%>%
                               lapply(function(x)x%>%
                                        rename_all(~paste0('sp', .x))),
                             saturated = lapply(unsaturated_comm, function(x)x$saturated)%>%
                               lapply(function(x)x%>%
                                        rename_all(~paste0('sp', .x))))


#long tibble
un_saturated_theta <- un_saturated_tibble%>%
  pivot_longer(cols = 3:4, names_to = 'spe_comp')%>%
  mutate(theta.out = parallel_theta(spe_comp = spe_comp, metrics = metrics)%>%
           nest())


#Parallelize
parallel_theta <- function(spe_comp, metrics){
  require(parallel)
  
  cl <- makeCluster (6)
  theta.out <- c()
  spec.comp <- spe_comp
  clusterExport (cl, c ('theta.out', 'spec.comp','metrics'))
  theta.output <- parLapply(cl, spec.comp, function(x){
    
    
    #loop through metrics and apply each to theta function
    for (j in rownames(metrics)){
      
      theta.out[[j]] <- genspe::calculate.theta(x,
                                                method = metrics[j,"method"],
                                                beta.div.method = metrics[j,"beta.div.method"],
                                                beta.div.sqrt.D = metrics[j, "beta.div.sqrt.D"],
                                                force.subsample = metrics[j, "force.subsample"],
                                                reps =100)}
    
    #bind rows with new column as metric name
    out <- dplyr::bind_rows(theta.out, .id = 'metric')
    
    return(out)
  })
  
  stopCluster (cl)
  return(theta.output)
}


#Calculate theta in parallel ------
unsaturated_theta <- un_saturated_tibble$unaturated%>%
  parallel_theta(., metrics = metrics)

saturated_theta <- un_saturated_tibble$saturated%>%
  parallel_theta(., metrics = metrics)


#Combine into a single tibble -------
community_tibble <- un_saturated_tibble%>%
  pivot_longer(cols = c(3,4), names_to = 'saturation', values_to = 'community_matrix')

theta_tibble <- un_saturated_tibble%>%
  select(-unsaturated, -saturated)%>%
  mutate(unsaturated_theta = unsaturated_theta)%>%
  mutate(saturated_theta = saturated_theta)%>%
  pivot_longer(cols = 3:4, names_to = 'saturation', values_to = 'theta_table')%>%
  mutate(saturation = plyr::mapvalues(saturation, from = 'unsaturated_theta', to = 'unsaturated'))%>%
  mutate(saturation = plyr::mapvalues(saturation, from = 'saturated_theta', to = 'saturated'))
  
  
#Complete tibble

community_theta_tibble <- community_tibble%>%
  left_join(theta_tibble)

saveRDS(community_theta_tibble, file = 'theta_tibble.rds')
