library(parallel)
library(dplyr)
library(mcomsimr)

###############################################################
# Metacommunity simulation (based on Thompson et al 2020 Ecology Letters)
###############################################################

#Metacommunity initiation
metacomm_sim <- function(no.patch = 300, no.spe = 10, sp.optima = 'uniform',
                         time_steps = 1, burn_in = 50,
                         init = 20, autocorr = 500){
  require(dplyr)
  require(mcomsimr)
  require(scales)
  #New function remove progressbar for parallelization 
  m_sim_new <- function (patches, species, dispersal = 0.01, plot = FALSE, torus = FALSE, 
                         kernel_exp = 0.1, env1Scale = 500, timesteps = 1200, burn_in = 800, 
                         initialization = 200, max_r = 5, min_env = 0, max_env = 1, 
                         env_niche_breadth = 0.5, optima_spacing = "random", intra = 1, 
                         min_inter = 0, max_inter = 1, comp_scaler = 0.05, extirp_prob = 0, 
                         landscape, disp_mat, env.df, env_optima, int_mat) 
  {
    if (missing(landscape)) {
      landscape <- landscape_generate(patches = patches, plot = plot)
    }
    else {
      landscape <- landscape_generate(patches = patches, xy = landscape, 
                                      plot = plot)
    }
    if (missing(disp_mat)) {
      disp_mat <- dispersal_matrix(landscape = landscape, torus = torus, 
                                   kernel_exp = kernel_exp, plot = plot)
    }
    else {
      disp_mat <- dispersal_matrix(landscape = landscape, disp_mat = disp_mat, 
                                   torus = torus, kernel_exp = kernel_exp, plot = plot)
    }
    if (missing(env.df)) {
      env.df <- env_generate(landscape = landscape, env1Scale = env1Scale, 
                             timesteps = timesteps + burn_in, plot = plot)
    }
    else {
      env.df <- env_generate(landscape = landscape, env.df = env.df, 
                             env1Scale = env1Scale, timesteps = timesteps + burn_in, 
                             plot = plot)
    }
    if (missing(env_optima)) {
      env_traits.df <- env_traits(species = species, max_r = max_r, 
                                  min_env = min_env, max_env = max_env, env_niche_breadth = env_niche_breadth, 
                                  optima_spacing = optima_spacing, plot = plot)
    }
    else {
      env_traits.df <- env_traits(species = species, max_r = max_r, 
                                  min_env = min_env, max_env = max_env, env_niche_breadth = env_niche_breadth, 
                                  optima_spacing = optima_spacing, optima = env_optima, 
                                  plot = plot)
    }
    if (missing(int_mat)) {
      int_mat <- species_int_mat(species = species, intra = intra, 
                                 min_inter = min_inter, max_inter = max_inter, comp_scaler = comp_scaler, 
                                 plot = FALSE)
    }
    else {
      int_mat <- species_int_mat(species = species, int_mat = int_mat, 
                                 intra = intra, min_inter = min_inter, max_inter = max_inter, 
                                 comp_scaler = comp_scaler, plot = FALSE)
    }
    dynamics.df <- data.frame()
    N <- matrix(rpois(n = species * patches, lambda = 0.5), nrow = patches, 
                ncol = species)
    #pb <- txtProgressBar(min = 0, max = initialization + burn_in + 
    #                       timesteps, style = 3)
    for (i in 1:(initialization + burn_in + timesteps)) {
      if (i <= initialization) {
        if (i %in% seq(10, 100, by = 10)) {
          N <- N + matrix(rpois(n = species * patches, 
                                lambda = 0.5), nrow = patches, ncol = species)
        }
        env <- env.df$env1[env.df$time == 1]
      }
      else {
        env <- env.df$env1[env.df$time == (i - initialization)]
      }
      r <- max_r * exp(-(t((env_traits.df$optima - matrix(rep(env, 
                                                              each = species), nrow = species, ncol = patches))/(2 * 
                                                                                                                   env_traits.df$env_niche_breadth)))^2)
      N_hat <- N * r/(1 + N %*% int_mat)
      N_hat[N_hat < 0] <- 0
      N_hat <- matrix(rpois(n = species * patches, lambda = N_hat), 
                      ncol = species, nrow = patches)
      E <- matrix(rbinom(n = patches * species, size = N_hat, 
                         prob = rep(dispersal, each = patches)), nrow = patches, 
                  ncol = species)
      dispSP <- colSums(E)
      I_hat_raw <- disp_mat %*% E
      I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
      I_hat[is.nan(I_hat)] <- 1
      I <- sapply(1:species, function(x) {
        if (dispSP[x] > 0) {
          table(factor(sample(x = patches, size = dispSP[x], 
                              replace = TRUE, prob = I_hat[, x]), levels = 1:patches))
        }
        else {
          rep(0, patches)
        }
      })
      N <- N_hat - E + I
      N[rbinom(n = species * patches, size = 1, prob = extirp_prob) > 
          0] <- 0
      dynamics.df <- rbind(dynamics.df, data.frame(N = c(N), 
                                                   patch = 1:patches, species = rep(1:species, each = patches), 
                                                   env = env, time = i - initialization - burn_in))
      #setTxtProgressBar(pb, i)
    }
    #close(pb)
    dynamics.df <- left_join(dynamics.df, env_traits.df)%>%
      filter(N>0)%>%
      filter(time>0)
    
    env.df$time_run <- env.df$time - burn_in
    env.df_init <- data.frame(env1 = env.df$env1[env.df$time == 
                                                   1], patch = 1:patches, time = NA, time_run = rep(seq(-(burn_in + 
                                                                                                            initialization), -burn_in, by = 1), each = patches))
    env.df <- rbind(env.df_init, env.df)
    if (plot == TRUE) {
      sample_patches <- sample(1:patches, size = min(c(patches, 
                                                       6)), replace = FALSE)
      g <- dynamics.df %>% filter(time %in% seq(min(dynamics.df$time), 
                                                max(dynamics.df$time), by = 10)) %>% filter(patch %in% 
                                                                                              sample_patches) %>% ggplot(aes(x = time, y = N, group = species, 
                                                                                                                             color = optima)) + geom_line() + facet_wrap(~patch) + 
        scale_color_viridis_c() + geom_path(data = filter(env.df, 
                                                          patch %in% sample_patches), aes(y = -5, x = time_run, 
                                                                                          color = env1, group = NULL), size = 3)
      print(g)
    }
    return(list(dynamics.df = dynamics.df, landscape = landscape, 
                env.df = env.df, env_traits.df = env_traits.df, disp_mat = disp_mat, 
                int_mat = int_mat))
  }
  
  
  #require(dplyr)
  #require(ggplot2)
  #require(mcomsimr)
  
  #Total simulation time
  total_time <- time_steps + burn_in
  
  #Initialize landscape and envrionmental gradient -----------------
  #Customized landscape generation function
  custom_land_generation <- function(patches = 100){
    x <- sample(1:patches)
    y <- sample((patches+1) : (patches*2))
    land <- data.frame(x = x, y = y)%>%
      tibble::rowid_to_column(var = 'patch')
    return(land)
  }
  
  #Landscape
  land_sim <- custom_land_generation(patches = no.patch)
  
  #Envrionmental gradient with spatial autocorrelation (env1Scale argument)
  env_sim <- env_generate(landscape = land_sim,timesteps = 1,
                          env1Scale = autocorr, plot = FALSE)
  
  #Make environment constant over time
  stable_env <- lapply(1:total_time, function(x)x=env_sim[,1:2])%>%
    bind_rows(.id = 'time')%>%
    select(env1, patch, time)%>%
    mutate(time = as.integer(time))
  
  #Define species attributes ------------------
  #Niche width --------------
  #niche_width <- runif(min = 0.1, max = 0.5, n = no.spe) #species niche width are drawn from uniform distirbuion 
  niche_width <- rexp(no.spe, rate = 10^(-3))%>%scales::rescale(to = c(0,0.5))
  
  #Growth rates --------------
  spe_r <- rnorm(n = no.spe, mean = 5) #normal distribution with mean 5
  
  #Species optima --------------
  if (sp.optima == 'uniform'){
    optima <- runif(min = 0.2, max = 0.8,n = no.spe) #constrained within 0.2 and 0.8
  }
  if (sp.optima == 'normal'){
    optima <- rnorm(n = no.spe, mean = 0.5)%>% #constrained within 0.2 and 0.8
      scales::rescale(to = c(0.2,0.8))
  }
  if (sp.optima == 'skew'){
    optima <- rbeta(n = no.spe, shape1 = 1, shape2 = 5)%>%
      scales::rescale(to = c(0.2,0.8))
  }
  
  #Competition structure---------
  int_mat_coexist <- species_int_mat(species = no.spe, intra = 1,
                                     min_inter = 0, max_inter = 0.99, plot = FALSE)
  
  
  #SIMULATE----------
  sim <- m_sim_new(timesteps = time_steps, burn_in = burn_in, initialization = init,
                     patches = no.patch,
                     species = no.spe, max_r = spe_r,
                     env_optima = optima, env_niche_breadth = niche_width,
                     dispersal = 0.01, int_mat = int_mat_coexist,
                     env.df = stable_env, plot = FALSE)
  
  return(sim)
}


#Visualizing fundamental niches
#set up dataframe
density_ind_pars_skew <- env_traits(species = no.spe, max_r = spe_r, min_env = 0, max_env = 1,
                                    env_niche_breadth = niche_width, optima = skew_optim, plot = FALSE)
#plot function
sp_resp_curve <- function(env.trait.df, gr.length = 100){
  no.spe <- length(env.trait.df$species)
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  
  out <- sapply(X = 1:no.spe, FUN = function(x) {
    exp(-((env.trait.df$optima[x] - seq(0, 1, length = gr.length))/(2*env.trait.df$env_niche_breadth[x]))^2)
  })*rep(env.trait.df$max_r, each = gr.length)
  
  plot.out <- out%>%
    as_tibble()%>%
    tibble::rowid_to_column(var = 'Gradient')%>%
    pivot_longer(cols = -Gradient, names_to = 'species')%>%
    mutate(species = as.integer(substr(species, 2, nchar(species))))%>%
    left_join(par_sim[[1]]$env_traits.df%>%select(species, env_niche_breadth), by = 'species')%>%
    ggplot(aes(x = Gradient, y = value, group = species))+
    geom_line(alpha = 0.8, aes(col = env_niche_breadth), size = 2)+
    labs(y = 'r')+
    geom_point(data = env.trait.df, aes(x = optima*gr.length, y = 0),
               alpha = 0.5, size = 7, pch = '|')+guides(col = FALSE)+
    theme_classic()
  
  list.out <- list(out, plot.out)
  
  return(list.out)
}

#Realized niche----------------
real_fund_niche_fun <- function(dynamics.df, env_traits.df){
  require(dplyr)
  
  df <- dynamics.df%>%
    group_by(species)%>%
    summarise(iqr = IQR(env))%>%
    left_join(env_traits.df, by = 'species')
  
  g <- df%>%
    ggplot(aes(x = iqr, y = env_niche_breadth))+
    geom_abline(slope = 1, size = 1, linetype = 'dashed')+
    geom_point(size = 5, alpha = 0.7)+
    labs(x = 'Realized niche width (IQR)', y = "Fundamental niche width")+
    theme_classic()
  
  out <- list(table = df, plot = g)
  return(out)
  
}


#####################################
#Simulate 4 independent runs: use 4 cores
#Each optima scenario contains 4 replicates (replicate runs are parallelized)
optima_scenario <- c('uniform','normal','skew')
#clusterExport(cl, list('metacomm_sim', 'patches', 'optima_scenario'))
#run each scenario 4 times
optima_simulations <- list()
for (i in 1:length(optima_scenario)){
  optim <- optima_scenario[i]
  
  #Parallel
  cl <- makeCluster(4)  
  clusterExport(cl, list('metacomm_sim', 'optim', 'optima_scenario'))
  #Long
  #optima_simulations[[i]] <- parLapply(cl, 1:4, function(x)metacomm_sim(sp.optima =  optim,
  #                                                           no.patch = 500, no.spe = 50, burn_in = 100))
  #short
  optima_simulations[[i]] <- parLapply(cl, 1:4, function(x)metacomm_sim(sp.optima =  optim,
                                                                        no.patch = 100, no.spe = 10,
                                                                        burn_in = 10))
  
  
  stopCluster(cl)  
  
}




