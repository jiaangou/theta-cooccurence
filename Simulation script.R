#Scenario setup
library(dplyr)
library(tidyr)
library(ggplot2)

##############################
#Useful functions for data wrangling
##############################
#Community matrix: converts dynamics.df to community matrix
community_mat_fun <- function(dynamics.df){
  
  #converts simmulate_MC output into sampleXspecies matrix with time points split into list elements
  
  require(dplyr)
  require(tidyr)
  
  
  comm_mat <- dynamics.df%>%
    filter(time > 0)%>%
    select(time, patch, env, N, species)%>%
    pivot_wider(id_cols = c(patch, env, time), names_from = species, values_from = N, values_fill = 0)
  
  
  return(comm_mat)
  
}


###############################################################
#1. Uneven sampling: 1) Skew (beta distribution) 2) Random
# Note: this is set PRIOR to metacommunity simulation
###############################################################
#Generate sampling probabilities along env gradient
sampling_eveness <- function(env_df, shape1 = 0.5, shape2 = 5){

  #bin and probability dataframe
  env_prob <- data.frame(env1 = env_df$env1%>%sort)%>%
    mutate(skew_prob = rbeta(n = n(), shape1 = shape1, shape2 = shape2)%>%sort) #skew probability
  
  return(env_prob)
}

#Example:
#sample(df$env1, size = 1000, replace = TRUE, prob = df$skew_prob)%>%
#  hist(breaks = 100)

#Visualization function
sampling_visual <- function(env_df, sample_out, sample_size = 1000){
  
  gg_sample_evennes <- function(env_df, skew_samp, uni_samp){
    require(ggplot2)
    env_df%>%
      ggplot(aes(x = env1, y = ..scaled..))+geom_density()+
      geom_density(data = skew_samp, aes(x = env1, y = ..scaled..), fill = 'red', alpha = 0.3)+
      geom_density(data = uni_samp, aes(x = env1, y = ..scaled..), fill = 'blue',alpha = 0.3)+
      theme_classic()
  }
  
  skew_samp <- sample(x = sample_out$env1, size = sample_size, replace = TRUE)%>%
    data.frame(env1 = .)
  
  uni_samp <- sample(x = sample_out$env1, size = sample_size, replace = TRUE, prob = samp_even$skew_prob)%>%
    data.frame(env1 = .)
  
  p <- gg_sample_evennes(env_df = env_df, skew_samp = skew_samp, uni_samp = uni_samp)
  
  samples <- data.frame(random = uni_samp$env1, skew = skew_samp$env1)
  
  out <- list(p, samples)
  return(out)
}


#Example: 
#samp_even <- sampling_eveness(env_df = env_df, shape1 = 0.8)
#sampling_visual(env_df = env_df, sample_out = samp_even, sample_size = 100)


###############################################################
#2. Species optima: i) Skew, ii) Random
#Note: 1) restrict optima within a range, 2) is set PRIOR to metacommunity simulation
###############################################################

optima_distribution <- function(no.spe, skew = c('left','right')){
  
  #Skew: left or right
  if(skew == 'left'){
  skew_optim <- rbeta(no.spe, shape1 = 1.5, shape2 = 5)
      
  }else if(skew == 'right'){
  skew_optim  <- rbeta(no.spe, shape1 = 5, shape2 = 1.5)
      
  }else{
    stop("Specify skew direction")
  }
    

  
  #uniform
  uni_optim <- runif(no.spe)
  
  #normal
  norm_optim <- rnorm(no.spe)%>%
    scales::rescale(to = c(0,1))
  
  out <- list(skew = skew_optim, uniform = uni_optim, normal = norm_optim)
  
  return(out)
  }


###############################################################
#3. Species pool: i) Non saturated, ii) Saturated
#Unsaturated/Linear  - sample 50% of species pool (Manthey & Fridley 2008)
#Saturated/Curvelinear - sample 50 individuals
#Note: this is done AFTER metacommunity simulation 
###############################################################


saturation_sampling <- function(dynamics.df, k = 0.5, individuals = 100){
  
  #convert to matrix
  community_mat_fun <- function(dynamics.df){
    
    #converts simmulate_MC output into sampleXspecies matrix with time points split into list elements
    
    require(dplyr)
    require(tidyr)
    
    
    comm_mat <- dynamics.df%>%
      filter(time > 0)%>%
      select(time, patch, env, N, species)%>%
      pivot_wider(id_cols = c(patch, env, time), names_from = species, values_from = N, values_fill = 0)
    
    
    return(comm_mat)
    
  }
  
  sampleXspecies <- dynamics.df%>%
    community_mat_fun
  
  spe_comp <- sampleXspecies%>%
    select(-1:-3)
  
  env_gradient <- sampleXspecies%>%
    select(1:2)
  
  #i) Sampling by species
  unsaturated_comm <- function(spe_comp, k = 0.5){
    require(dplyr)
    require(vegan)
    #presnce absence
    spe_pa <- spe_comp>0
    
    #probabilities (relative abundances)
    sp_prob <- vegan::decostand(spe_comp,method = 'total')
    sp_prob_gradient <- list()
    for(i in seq(nrow(sp_prob))){
      sp_prob_gradient[[i]] <- sp_prob[i,][spe_pa[i,]]
    }
    
    #Species pool along gradient (species IDs)
    sp_pool_gradient <- sp_prob_gradient%>%
      lapply(names)
    
    #Species pool size along gradient (# of species)
    sp_pool_size_gradient <- sp_prob_gradient%>%
      lapply(length)
    
    #Subsample a fraction (k) of species pool: fraction is set as mean of poisson variable
    subsample_pool_n <- function(k = 0.5, pool_size){
      subset_size <- round(k * pool_size) #no stochasticity to ensure subset size always smaller than pool
      return(subset_size)
    }
    
    sp_fraction <- sp_pool_size_gradient%>%
      lapply(function(x)subsample_pool_n(pool_size = x, k = k))
    
    #metrics
    #metrics <- list(spe_pool = unlist(sp_pool_gradient),
    #                pool_size = unlist(sp_pool_size_gradient),
    #                sp_prob = unlist(sp_fraction))
    
    
    #unstaurated community output
    out <- mapply(x = sp_pool_gradient, y = sp_fraction, z = sp_prob_gradient,
                  function(x,y,z)sample(x = x, size = y, replace = FALSE, prob = z))%>%
      lapply(table)%>%
      bind_rows()%>%
      mutate_all(~replace(., is.na(.), 0))%>%
      as.data.frame
    
    #out.list <- list(unsaturated_comm  = out, metrics = metrics)
    return(out)
    
  }
  
  
  #ii) Sampling by individuals
  saturated_comm <- function(spe_comp, individuals = 100){
    require(dplyr)
    #Presence/absence 
    spe_pa <- spe_comp>0
    
    #Species abundances
    sp_abund_gradient <- list()
    for(i in seq(nrow(spe_pa))){
      sp_abund_gradient[[i]] <- spe_comp[i,][spe_pa[i,]]
    }
    
    #Species pool
    sp_pool_gradient <- sp_abund_gradient%>%
      lapply(names)
    
    #convert to number of individuals
    sp_individual_gradient <- mapply(x = sp_pool_gradient, y = sp_abund_gradient, function(x,y)rep(x = x, times = y))
    
    #Subsample by individuals
    #convert to number of individuals
    
    community_individuals <- sp_individual_gradient%>%
      lapply(function(x)sample(x, size = individuals, replace = TRUE))%>%
      lapply(table)%>%
      bind_rows%>%
      mutate_all(~replace(., is.na(.), 0))%>%
      as.data.frame
    
    return(community_individuals)
    
  }
  
  #pars <- list(k = k, individuals = 100)
  #list_funs <- list(unsaturated_comm, saturated_comm)
  
  
  unsat_comm <- spe_comp%>%
    unsaturated_comm(spe_comp = ., k = k)
  

  sat_comm <- spe_comp%>%
    saturated_comm(individuals = individuals)
  
  out <- list(unsaturated = unsat_comm, saturated = sat_comm, env_gradient = env_gradient)
  return(out)
  
}



###############################################################
# Function testing 
###############################################################
#Species pool sampling------
#Test function using default arguments and defined arguments
sat_test <- simulated_metacomms$uniform[[4]]$dynamics.df%>%
  saturation_sampling()
sat_test2 <- simulated_metacomms$uniform[[4]]$dynamics.df%>%
  saturation_sampling(k = 0.3, individuals = 50)


#Richness~env-----
#unsaturated
sat_test$unsaturated%>%
  rowSums%>%
  data.frame(richness = ., env = sat_test$env_gradient$env)%>%
  ggplot(aes(x = env, y = richness))+geom_point()
#saturated
rowSums(sat_test$saturated>0)%>%
  data.frame(richness = ., env = sat_test$env_gradient$env)%>%
  ggplot(aes(x = env, y = richness))+geom_point()




###############################################################
#Scenario sampling
###############################################################

#Load simulated parallel metacommunities
#Optima distributions: 1. 'uniform', 2. 'normal', 3. 'skew'
#4 parallel simulations for each optima distribution
simulated_metacomms <- readRDS('simulation_runs.rds')
names(simulated_metacomms) <- c('uniform','normal','skew')

#Double check optima distributions
simulated_metacomms[[1]][[1]]$env_traits.df$optima%>%hist
simulated_metacomms[[2]][[1]]$env_traits.df$optima%>%hist
simulated_metacomms[[3]][[1]]$env_traits.df$optima%>%hist


#Visualize individual/richness gradient
rich_ind_gradient <- function(dynamics.df){
  
  rich <- dynamics.df%>%
    group_by(patch, env)%>%
    summarise(richness = n())
  
  individual <- dynamics.df%>%
    group_by(patch, env)%>%
    summarise(individuals = sum(N))
 
  out <- rich%>%
    left_join(individual, by = c('patch','env'))
   
  return(out)
}

#Uniform
simulated_metacomms$uniform[[4]]$dynamics.df%>%
  rich_ind_gradient()%>%
  pivot_longer(cols = 3:4, names_to = 'index')%>%
  ggplot(aes(x = env, y = value))+geom_point(size = 3, alpha = 0.5)+
  facet_wrap(~index, scales = 'free')+
  labs(x = 'Environmental value', y = 'Count')+
  theme_bw()

#Skew
simulated_metacomms$skew[[4]]$dynamics.df%>%
  rich_ind_gradient()%>%
  pivot_longer(cols = 3:4, names_to = 'index')%>%
  ggplot(aes(x = env, y = value))+geom_point(size = 3, alpha = 0.5)+
  facet_wrap(~index, scales = 'free')+
  labs(x = 'Environmental value', y = 'Count')+
  theme_bw()



#Species pool sampling -------------------
library(parallel)
simulated_metacomms$uniform[[1]]$dynamics.df
sp_pool_sampling <- list()

simulated_metacomms[[1]]%>%length
simulated_metacomms%>%length
names(simulated_metacomms)[1]

for(i in 1:length(simulated_metacomms)){
  sp_pool_sampling[[i]] <- mclapply(simulated_metacomms[[i]],
                                    FUN = function(x)saturation_sampling(dynamics.df = x$dynamics.df),
                                    mc.cores = 4)
  names(sp_pool_sampling)[i] <- names(simulated_metacomms)[i]
}
#Check output
length(sp_pool_sampling[[1]])
sp_pool_sampling[[1]]

#Convert nested list to list-columns
sc_list_column <- sp_pool_sampling%>%
  setNames(nm = c('uniform', 'normal','skew'))%>% #name optimas 
  as_tibble()%>%
  pivot_longer(cols = 1:3, names_to = 'optima')%>%
  unnest_auto(value)%>%
  tibble::rowid_to_column(var = 'ID')%>%
  pivot_longer(cols = c('unsaturated','saturated'), names_to = 'saturation', values_to = 'sp_comp')


theta_list_column <- sc_list_column%>%
  mutate(theta_table = purrr::map(sp_comp, calculate.theta))


theta_list_column$theta_table[1]
theta_list_column$theta_table
#

#Post hoc analysis: calcualte theta using Fridley algorithm with different beta-diversity metrics
library(genspe)

#beta-diversity metrics
metrics <- read.delim ('theta metrics for evaluation_reduced_pres-abs.txt', row.names = 1, stringsAsFactors = F)
?calculate.theta



clusterExport (cl, c ('scenarios', 'metrics', 'reps'))
cor.ij <- parLapply (cl, sc, fun = function (sc.temp){
  cor.j <- vector ('numeric', length = nrow (metrics))
  names (cor.j) <- rownames (metrics)
  for (j in rownames (metrics))
  {
    temp.theta <- genspe::calculate.theta (sc.temp$a.mat, #Sp Comp matrix
                                           method = metrics[j,"method"], #BD method
                                           beta.div.method = metrics[j,"beta.div.method"], #metric class
                                           beta.div.sqrt.D = metrics[j, "beta.div.sqrt.D"], #sqrt bd 
                                           force.subsample = metrics[j, "force.subsample"], #
                                           reps = reps)
    cor.j[j] <- cor (sc.temp$simul.comm$range[temp.theta$sci.name], temp.theta$theta, method = 'spearman')
  }
  cor.j
})
?calculate.theta





for (j in rownames (metrics))
{
  temp.theta <- theta::calculate.theta (sc.temp$a.mat,
                                        method = metrics[j,"method"],
                                        beta.div.method = metrics[j,"beta.div.method"],
                                        beta.div.sqrt.D = metrics[j, "beta.div.sqrt.D"],
                                        force.subsample = metrics[j, "force.subsample"],
                                        pa.transform = metrics[j, "pa.transform"],
                                        q = metrics [j, 'q'],
                                        reps = reps)
  cor.j[j] <- cor (sc.temp$simul.comm$range[temp.theta$sci.name], temp.theta$theta, method = 'spearman')
}
temp.theta <- theta::calculate.theta (sc.temp$a.mat, #species composition matrix
                                      method = metrics[j,"method"], #BD method
                                      beta.div.method = metrics[j,"beta.div.method"], #
                                      beta.div.sqrt.D = metrics[j, "beta.div.sqrt.D"],
                                      force.subsample = metrics[j, "force.subsample"],
                                      pa.transform = metrics[j, "pa.transform"],
                                      q = metrics [j, 'q'],
                                      psample = 5, #Threshold species frequency for calculating theta
                                      reps = reps) #number of random subsamples used to 



theta_out <- tibble_sample$uniform$`1`$unsaturated%>%
  calculate.theta()


simulated_metacomms$uniform[[1]]$env_traits.df$species

data.frame(species = as.integer(theta_out$sci.name), theta = theta_out$theta)%>%
  left_join(simulated_metacomms$uniform[[1]]$env_traits.df, by = 'species')%>%
  ggplot(aes(x = theta, y = env_niche_breadth))+geom_point(size = 3)


?calculate.theta



