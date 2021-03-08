# Evaluation of individual theta metrics using different scenarios of simulate data
setwd ('c:\\Users\\zeleny\\Documents\\specialists and generalists on the whole database\\scripts\\')
scenarios <- read.delim ('scenarios of simulated model for evaluation of theta metrics.txt', row.names = 1, stringsAsFactors = F)
metrics <- read.delim ('theta metrics for evaluation_reduced.txt', row.names = 1, stringsAsFactors = F)
#devtools::install_github ('zdealveindy/genspe')
library (genspe)
totS <- 300
Np <- 300
reps <- 100

sc <- list ()
sce <- rep (rownames (scenarios), each = 10)
for (i in seq (1, length (sce)))
{
  which.scenario <- sce[i]
  sc[[i]] <- sample.comm (simul.comm = simul.comm (totS = totS, niche.type = scenarios[which.scenario,]$niche.type, spec.optima = scenarios[which.scenario,]$spec.optima), Np = Np, sample.x = scenarios[which.scenario,]$sample.x, based.on = scenarios[which.scenario,]$based.on)
}

cl <- makeCluster (6)
clusterExport (cl, c ('scenarios', 'metrics', 'reps'))
cor.ij <- parLapply (cl, sc, fun = function (sc.temp){
  cor.j <- vector ('numeric', length = nrow (metrics))
  names (cor.j) <- rownames (metrics)
  for (j in rownames (metrics))
  {
    temp.theta <- genspe::calculate.theta (sc.temp$a.mat, method = metrics[j,"method"], beta.div.method = metrics[j,"beta.div.method"], beta.div.sqrt.D = metrics[j, "beta.div.sqrt.D"], force.subsample = metrics[j, "force.subsample"], reps = reps)
    cor.j[j] <- cor (sc.temp$simul.comm$range[temp.theta$sci.name], temp.theta$theta, method = 'spearman')
  }
  cor.j
})
stopCluster (cl)

cor.ij.df <- do.call (rbind.data.frame, cor.ij)
names (cor.ij.df) <- names (cor.ij[[1]])
save (cor.ij.df, file = 'cor.ij.df.2.r')
