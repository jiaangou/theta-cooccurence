#Realized and fundamental niches
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


#Load data -------
simulation_tibble <- here('data','simulation_data.rds')%>%
  readRDS()
theta_tibble <- readRDS(here('data','theta_tibble.rds'))
#spp_parameters <- readRDS(here('data','Spp_parameters.rds'))
metrics <- here('theta metrics for evaluation.txt')%>%
  read.delim(., row.names = 1, stringsAsFactors = F)



#Extract fundamental and realized niches -------
spe_parameters <- function(x){
  out <- x%>%
  select(species, optima, env_niche_breadth, max_r)%>%
  distinct()%>%
  mutate(species = paste0('sp',species))
return(out)
}
#Fundamental niche: defined a priori in the model 
#Realized niche calculation: Option 1 = IQR/percentiles; Option 2 = SD; 

#NOTE: realized niches for saturated VS unsaturated are the same as it is calculated from raw community matrix data (ie. before samplingby species)

#Realized niche needs to be weighted by abundance -----

#Custom function to weight env by abundance to be used for group_by summaries
weighted_realized_niche <- function(env, N, SD = TRUE, quantile_p = c(0.05, 0.95)){
  #replicate each env by number of individuals 
  weighted_env <- rep(env, times = N)
  
  if(SD == TRUE){
    out <- weighted_env%>%
      sd
    
  }
  else{
    quant <- weighted_env%>%
      quantile(probs = quantile_p)
    
    out <- quant[2] - quant[1]

  }

  return(out)
}

#calcualte realized niche and join fundamental niche into same tibble
niche_tibble <- simulation_tibble%>%
  tibble(., realized_niche = lapply(simulation_tibble$data, FUN = function(x)x%>%
           mutate(species = paste0('sp',species))%>%
           group_by(species)%>%
           summarise(quantile_range = weighted_realized_niche(env = env, N = N, SD = FALSE),
                     sd = weighted_realized_niche(env = env, N = N, SD = TRUE))),
         fundamental_niche = lapply(simulation_tibble$data, FUN = function(x)x%>%
                                      spe_parameters))


#Combine tibbles -----
complete_tibble <- theta_tibble%>%
  left_join(niche_tibble, by = c('optima_distribution','simulation'))%>%
  mutate(theta_table = theta_table%>%
           lapply(function(x)x%>%
                    rename(`species` = sci.name)))

#Unnest theta and realized niche separately -------

#Theta
theta_df <- complete_tibble%>%
  unnest(theta_table)%>%
  select(optima_distribution, simulation, saturation, metric, species, theta)

#Niches
#I. Realized
realized_df <- complete_tibble%>%
  unnest(realized_niche)%>%
  select(optima_distribution, simulation, saturation,species, quantile_range, sd)

#II. Fundamental (include fundamental niche for completeness)
fundamental_df <- complete_tibble%>%
  unnest(fundamental_niche)%>%
  select(optima_distribution, simulation, saturation,species, optima, env_niche_breadth)


#Join by optima distribution, simulation, saturation, species
theta_realized_df <- theta_df%>%
  left_join(realized_df, by = c('optima_distribution', 'simulation', 'saturation', 'species'))%>%
  left_join(fundamental_df, by = c('optima_distribution', 'simulation', 'saturation', 'species'))

#Plot 1: Correlation between fundamental and realized niche --------
realized_fundamental_p <- theta_realized_df%>%
  mutate(simulation = as.factor(simulation))%>%
  ggplot(aes(y = env_niche_breadth, x = sd, col = simulation,
             group = interaction(simulation, optima_distribution, saturation)))+
  geom_point(size = 2)+
  scale_color_manual(values = RColorBrewer::brewer.pal(4, name = 'Set2'),
                     name = "Simulation replicate")+
  geom_smooth(method = 'lm', se = FALSE)+
  facet_grid(optima_distribution~saturation)+
  labs(x = 'Realized niche (SD)', y = "Fundamental niche")+
  ylim(0, 0.5)+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
ggsave(file = 'niche-correlation.png')




#Calculate correlation between theta and realized niche by optima distirbution X simulation X saturation X metric
correlation_tibble <- theta_realized_df%>%
  group_by(optima_distribution, simulation, saturation, metric)%>%
  summarise(theta_realized_cor = cor(theta, sd, method = 'spearman'),
            theta_fundamental_cor = cor(theta, env_niche_breadth, method = 'spearman'))

#Exploratory analysis ------
#I. Plot 2: Correlation between theta-fundamental correlation VS theta-realized -----
niche_theta_correlations_p <- correlation_tibble%>%
  group_by(optima_distribution, saturation, metric)%>%
  summarise(theta_real_cor_mu = mean(theta_realized_cor), theta_fund_cor_mu = mean(theta_fundamental_cor))%>%
  ggplot(aes(x = theta_real_cor_mu, y = theta_fund_cor_mu, col = optima_distribution))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = RColorBrewer::brewer.pal(3, name = 'Set2'),
                     name = 'Optima distribution')+
  facet_wrap(~saturation)+
  labs(x = "Theta-Realized r", y = "Theta-Fundamental r")+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 12))

ggsave(file = 'theta-niche-correlation.png', width = 8, height = 5, units = 'in')

#II. PCA with metric as "species" (using just theta_realized correlation)
library(vegan)
#convert to wide by pivotting metric
correlation_wide <- correlation_tibble%>%
  pivot_longer(cols = c('theta_realized_cor', 'theta_fundamental_cor'), names_to = 'niche', values_to = 'r')%>%
  mutate(niche = plyr::mapvalues(niche, from = 'theta_realized_cor', to = 'realized'))%>%
  mutate(niche = plyr::mapvalues(niche, from = 'theta_fundamental_cor', to = 'fundamental'))%>%
  pivot_wider(names_from = "metric", values_from = "r")
  
#Calculate PCA
PCA <- correlation_wide%>%
  ungroup()%>%
  select(-optima_distribution, -simulation, -saturation, -niche)%>%
  rda()

#Axes explained variance
pca_explained <- eigenvals(PCA)/PCA$tot.chi

#Get data using scaling = 1
pca_scores <- scores(PCA, scaling =1)%>% #Species are more clear for scaling = 1
  lapply(.,as.data.frame)

#Metric groups
metric_df <- pca_scores$species%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var = 'metric')%>%
  mutate(metric_group = stringr::str_split(metrics$method, "[.]")%>%
         lapply(function(x)x[1])%>%
         unlist())

#Sample scores
sample_df <- correlation_wide%>%
  select(optima_distribution, simulation, saturation, niche)%>%
  cbind(., pca_scores$sites)


  
#Plot 3: PCA of BD metrics and scenarios --------- 
theta_pca_p <- sample_df%>%  
  ggplot(data = .)+
  geom_hline(aes(yintercept = 0),linetype = 'dashed')+
  geom_vline(aes(xintercept = 0), linetype = 'dashed')+
  geom_point(aes(x = PC1, y = PC2, pch = saturation), size = 3, alpha = 0.7)+
  scale_shape_manual(values = c(1,16), name = 'Saturation')+
  geom_segment(data = metric_df, 
               aes(x=0, xend = PC1, y =0, yend = PC2, col = metric_group),
               arrow = arrow(length = unit(0.2,'cm')),inherit.aes = FALSE)+
  scale_color_manual(values = RColorBrewer::brewer.pal('Dark2', n = 7),name = 'Metric category')+
  geom_text(data = pca_scores$species, aes(x= PC1, y = PC2),
            label = rownames(pca_scores$species))+
  labs(x = paste0('PC1',' (', round(pca_explained[1]*100,2), '%)'),
       y = paste0('PC2',' (', round(pca_explained[2]*100,2), '%)'))+
  ylim(-1,1) +
  theme_bw() + 
  #theme(legend.position = c(.9, .75))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size = 13))

#quartz()
#theta_pca_p
ggsave(file = 'theta-PCA.png')

#Rank plot
sce.rank.plot <- function(data, sce ='sc01'){
  p <- data%>%
    filter(scenario == sce)%>%
    arrange(desc(value))%>%
    mutate(metric = factor(metric, levels = metric))%>%
    ggplot(aes(x= metric,y = value,col = metric.group))+
    geom_point(size = 3)+ggtitle(label = sce)+
    labs(y = expression(italic('r')), x = "Metric")+
    ylim(0,1)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))
  return(p)
}

correlation_tibble%>%
  arrange(desc(theta_realized_cor))%>%
  ggplot(aes(x = metric, y = theta_realized_cor))+
  geom_point()+
  facet_grid(optima_distribution~saturation)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))

#Get summary stats for each group at i)  metric and ii) across metric level
#i) metric level
correlation_summ <- correlation_tibble%>%
  pivot_longer(cols = c('theta_realized_cor','theta_fundamental_cor'),
               names_to = 'niche', values_to = 'r')%>%
  mutate(niche = plyr::mapvalues(x = niche, from = 'theta_fundamental_cor',
                                 to = 'Fundamental'))%>%
  mutate(niche = plyr::mapvalues(x = niche, from = 'theta_realized_cor',
                                 to = 'Realized'))%>%
  group_by(optima_distribution, saturation, metric, niche)%>%
  summarise(mean_r = mean(r), sd_r = sd(r))

#ii) across metric
scenario_summ <- correlation_tibble%>%
  pivot_longer(cols = c('theta_realized_cor','theta_fundamental_cor'),
               names_to = 'niche', values_to = 'r')%>%
  mutate(niche = plyr::mapvalues(x = niche, from = 'theta_fundamental_cor',
                                 to = 'Fundamental'))%>%
  mutate(niche = plyr::mapvalues(x = niche, from = 'theta_realized_cor',
                                 to = 'Realized'))%>%
  group_by(optima_distribution, saturation, niche)%>%
  summarise(mean_r = mean(r), sd_r = sd(r))

#Plot 4: Theta (realized) performance
performance_realized_p <- correlation_summ%>%
  filter(niche == 'Realized')%>%
  ggplot(aes(x = metric, y = mean_r))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r), width = 0.3)+
  facet_grid(optima_distribution~saturation)+
  labs(x = 'Metric', y = 'Theta-Realized r')+
  geom_hline(data = scenario_summ%>%filter(niche == 'Realized'),
             aes(yintercept = mean_r), linetype = 'dashed')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))

ggsave(plot = performance_realized_p, file = 'realized-theta-performance.png')


#Plot 5: Theta (fundamental) performance
performance_fun_p <- correlation_summ%>%
  filter(niche == 'Fundamental')%>%
  ggplot(aes(x = metric, y = mean_r))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r), width = 0.3)+
  facet_grid(optima_distribution~saturation)+
  labs(x = 'Metric', y = 'Theta-Fundamental r')+
  geom_hline(data = scenario_summ%>%filter(niche == 'Fundamental'),
             aes(yintercept = mean_r), linetype = 'dashed')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))

ggsave(plot = performance_fun_p, file = 'fundamental-theta-performance.png')

