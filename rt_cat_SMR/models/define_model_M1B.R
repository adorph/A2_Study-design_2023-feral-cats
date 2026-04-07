
####                    BUILD MODEL FOR NIMBLE                         ####
####                                                                       
####  This code was compiled based on repositories developed by Ben        
####  Augustine on github (from: https://github.com/benaug) and code from
####  Herrmann (2022) (https://doi.org/10.5281/zenodo.7092971).
####                                                                       
####  Modified here by A Dorph.                                                 
####                                                                       

NimModel <- nimbleCode({
  
  # Detection function priors:
  lam0.beta0 ~ dnorm(-3.5, 1 / ( 1.5 ^ 2)) 
  sigma.beta0 ~ dnorm(log(496),  1 / ( 0.578 ^ 2))# From data & literature, more diffuse
  
  # Data augmentation prior
  psi ~ dbeta(2,2) #dunif(0,1) 
  
  # Categorical ID covariate priors
  thin.beta0 ~ 	dnorm(0, 1)  # Intercept - black, bicolor false, long hair false - keeps values centered near 0 but allows variation.
  for(i in 1:4){
    thin.beta.coat[i] ~ dnorm(0, 1) 
  }
  thin.beta.bi ~ dnorm(0, 1)  # Bicolor TRUE offset

  for(m in 1:n.cat){
    alpha[m,1:n.levels[m]] <- 1 # prior parameters
    gammaMat[m,1:n.levels[m]]~ddirch(alpha[m,1:n.levels[m]])
  }
  
  # For the spatial mask area 
  lambda.cell[1:n.cells] <- InSS[1:n.cells]
  # pi.cell is a simulated population of activity centers
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  # Generate likelihoods (except for s priors)
  for(i in 1:M) {
    # Latent inclusion indicator for each individual (1 = real, 0 = augmented).
    z[i] ~ dbern(psi)
    # True categorical covariates for each individual
    for(m in 1:n.cat){
      G.true[i,m]~dcat(gammaMat[m,1:n.levels[m]])
    }
    # Generate activity centres from within mask
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    # Get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(n.cells.x-(xlim[2] - s[i,1])/res[1])+1,trunc(n.cells.y-(ylim[2]-s[i,2])/res[2])+1]
    # Categorical likelihood for this cell, also disallowing s's in non-habitat:
    dummy.data[i] ~ dCell(pi.cell=pi.cell[s.cell[i]], InSS=InSS[s.cell[i]])
    # Detection model:
    for(j in 1:J){
      log(lam0[i,j]) <- lam0.beta0 
    }
    log(sigma[i]) <- sigma.beta0 
    # Calculate spatial detection rate:
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J, sigma=sigma[i], lam0=lam0[i,1:J], z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lambda=lam[i,1:J]*K1D[1:J],z=z[i])  # Model for complete capture histories (vectorized obs mod)
    # Create trait level indicators
    coat.indicator[i,1:4] <- GetCoat(G.true[i,2])
    bi.indicator[i] <- GetBi(G.true[i,3])
    # Conditional probability of assigning and ID to each detection
    logit(theta.i[i]) <- thin.beta0 + inprod(coat.indicator[i,1:4],thin.beta.coat[1:4]) + bi.indicator[i]*thin.beta.bi
    y.ID[i,1:J] ~ dBinomialVector(prob=theta.i[i], size=y.true[i,1:J], capcounts=capcounts[i])  # Model for ID process
  }
  # Number of times each individual was captured
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J])
  # Number of individuals sampled 
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],G.latent=G.latent[1:M,1:n.cat])
  # Sum for total abundance
  N <- sum(z[1:M])
})# end model
