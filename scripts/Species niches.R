#Realized and fundamental niches
library(here)
library(dplyr)
library(tidyr)




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


niche_tibble <- simulation_tibble%>%
  tibble(., realized_niche = lapply(simulation_tibble$data, FUN = function(x)x%>%
           mutate(species = paste0('sp',species))%>%
           group_by(species)%>%
           summarise(quantile_range = quantile(env, probs = 0.95) - quantile(env, probs = 0.05), sd = sd(env))),
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

#Calculate correlation between theta and realized niche by optima distirbution X simulation X saturation X metric
correlation_tibble <- theta_realized_df%>%
  group_by(optima_distribution, simulation, saturation, metric)%>%
  summarise(theta_realized_cor = cor(theta, sd, method = 'spearman'),
            theta_fundamental_cor = cor(theta, env_niche_breadth, method = 'spearman'))

#Exploratory analysis ------
#I. Correlation between theta-fundamental correlation VS theta-realized
correlation_tibble%>%
  group_by(optima_distribution, saturation, metric)%>%
  summarise(theta_real_cor_mu = mean(theta_realized_cor), theta_fund_cor_mu = mean(theta_fundamental_cor))%>%
  ggplot(aes(x = theta_real_cor_mu, y = theta_fund_cor_mu, col = optima_distribution))+
  geom_point()+
  facet_wrap(~saturation)

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


  
#plot 
sample_df%>%  
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
  theme(legend.position = c(.9, .75))




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


