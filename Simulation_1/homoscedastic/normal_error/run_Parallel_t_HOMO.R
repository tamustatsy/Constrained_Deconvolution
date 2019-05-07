library(batch)

parseCommandArgs()

seed_list <- 1000 + c(1:10)*100
n.ind_list <- c(1000, 5000, 10000, 15000)     #sample sizes for simulations

dir.create(paste(getwd(),"result", sep="/"))
for (i in 1:length(n.ind_list))  {
  dir.create(paste(getwd(),"result",n.ind_list[i], sep="/"))
}

for(seed in seed_list)
  for(n.ind in n.ind_list)
    rbatch("Constrained_Bayes_t_HOMO_doOne.R", seed = seed, n = n.ind)

rbatch.local.run()

