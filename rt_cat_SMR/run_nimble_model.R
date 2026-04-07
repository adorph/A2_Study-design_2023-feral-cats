run_nimble_model <- function(
    data_file,
    nimbuild_file,
    combo_name,
    chains = 5,
    iterations = 50000,
    burnin = 15000,
    model_script,
    output_dir = "./outputs/",
    sbc = FALSE,
    verbose = TRUE
) {
  
  require(nimble)
  require(coda)
  require(snow)
  require(doSNOW)
  require(foreach)
  
  tryCatch({
    
    # Load data
    if (!file.exists(data_file)) stop("Data file not found:", data_file)
    if (!file.exists(nimbuild_file)) stop("Nimbuild file not found:", nimbuild_file)
    
    load(data_file)
    load(nimbuild_file)
    
    if (verbose) cat("✓ Loaded data and nimbuild for", combo_name, "\n")
    
    # Setup parameters
    n.ID <- data$detections$n.ID
    gamma <- data$ID.covariates$gamma
    n.cat <- data$ID.covariates$n.cat
    n.levels <- data$ID.covariates$n.levels
    M1_factor <- 1.3 # more room for undetected individuals
    M2_factor <- 1
    seeds <- round(runif(chains, 1, 10000000))
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Parallel processing
    cl.tmp <- makeCluster(rep("localhost", chains), type = "SOCK")
    registerDoSNOW(cl.tmp)
    clusterExport(cl = cl.tmp, 
                  list("data", "n.ID", "gamma", "n.cat", "n.levels", "M1_factor", "M2_factor", 
                       "nimbuild", "model_script", "burnin"),
                  envir = environment())
    
    if (verbose) cat("✓ Starting parallel MCMC with", chains, "chains\n")
    
    start.time <- Sys.time()
    
    out <- foreach(c = 1:chains) %dopar% {
      library(nimble)
      source(paste0("code/rt_cat_SMR/models/", model_script))
      source("code/rt_cat_SMR/nimble_functions.R")
      source('code/rt_cat_SMR/nimble_restart_MCMC.R')
      
      # Augmentation
      M1 <- ceiling(n.ID / 2 * 16 * M1_factor)
      M2 <- ceiling(n.ID / 2 * 16 * M2_factor)
      M <- M1 + M2
      
      gammaMat <- nimbuild$gammaMat
      G.true.init <- nimbuild$G.true
      G.true.data <- G.true.init * NA
      
      # Build thin.beta parameters based on n.cat
      thin.beta.params <- list(thin.beta0 = -1.5, thin.beta.coat = c(0.25, 3, 2.5, 2.5), thin.beta.bi = 0.5)
      if (!grepl("Grid1000m", combo_name)) thin.beta.params$thin.beta.hair <- 0
      
      key.params <- list(lam0.beta0 = -3.5, sigma.beta0 = log(496))
      if (grepl("define_model_M2", model_script)) key.params$lam0.beta.sex <- 0
      if (grepl("define_model_M3", model_script)) key.params$sigma.beta.sex <- 0.2
      if (grepl("define_model_M4", model_script)) {
        key.params$lam0.beta.sex <- 0
        key.params$sigma.beta.sex <- 0.2
      }
        
      Niminits <- c(
        list(z = rep(1, M), s = nimbuild$s, G.true = G.true.init, ID = nimbuild$ID,
             capcounts = rowSums(nimbuild$y.true), y.true = nimbuild$y.true, gammaMat = gammaMat,
             G.latent = nimbuild$G.latent), #, psi = 0.9
        key.params, 
        thin.beta.params
      )
      
      constants <- list(M = M, J = data$trap.covariates$J, K1D = data$occasion.covariates$K1D,
                        n.samples = nimbuild$n.samples, n.cat = n.cat, n.levels = n.levels,
                        xlim = data$trap.covariates$xlim, ylim = data$trap.covariates$ylim,
                        n.cells = nimbuild$n.cells, n.cells.x = nimbuild$n.cells.x,
                        n.cells.y = nimbuild$n.cells.y, res = nimbuild$res)
      
      z.data <- c(rep(1, n.ID), rep(NA, M - n.ID))
      Nimdata <- list(y.true = matrix(NA, nrow = M, ncol = data$trap.covariates$J),
                      y.ID = nimbuild$y.ID, G.true = G.true.data, ID = rep(NA, nimbuild$n.samples),
                      z = z.data, X = as.matrix(data$trap.covariates$X), capcounts = rep(NA, M),
                      dummy.data = rep(0, M), cells = nimbuild$cells, InSS = nimbuild$InSS)
      
      # If running SCB use reduced parameters to speed up computation
      if (sbc) {
        base_pars <- c("lam0.beta0", "sigma.beta0", "psi", "N")
        extra_pars <- c(
          if (grepl("define_model_M2|define_model_M4", model_script)) "lam0.beta.sex",
          if (grepl("define_model_M3|define_model_M4", model_script)) "sigma.beta.sex"
        )
        parameters <- c(base_pars, extra_pars)
        parameters2 <- c("z") 
      } else {
        base_pars <- c('lam0.beta0', 'sigma.beta0', 'psi', 'N', 'n',
                        'gammaMat', 'thin.beta0', 'thin.beta.coat', 'thin.beta.bi')
        extra_pars <- c(
                  if (!grepl("Grid1000m", combo_name)) 'thin.beta.hair',
                  if (grepl("define_model_M2|define_model_M4", model_script)) "lam0.beta.sex",
                  if (grepl("define_model_M3|define_model_M4", model_script)) "sigma.beta.sex"
        )
        parameters <- c(base_pars, extra_pars)
        parameters2 <- c("s", "z", "s.cell", "lambda.cell") 
      }
      
      Rmodel <- nimbleModel(code = NimModel, constants = constants, data = Nimdata, inits = Niminits, check = FALSE)
      
      conf_nodes <- c("z", "gammaMat", "s", "psi", 'lam0.beta0', 'sigma.beta0', 'thin.beta0', 'thin.beta.coat', 'thin.beta.bi')
      if (!grepl("Grid1000m", combo_name)) conf_nodes <- c(conf_nodes, 'thin.beta.hair')
      if (grepl("define_model_M2", model_script)) conf_nodes <- c(conf_nodes, 'lam0.beta.sex')
      if (grepl("define_model_M3", model_script)) conf_nodes <- c(conf_nodes, 'sigma.beta.sex')
      if (grepl("define_model_M4", model_script)) {
        conf_nodes <- c(conf_nodes, 'lam0.beta.sex')
        conf_nodes <- c(conf_nodes, 'sigma.beta.sex')
      }
      
      conf <- configureMCMC(Rmodel, enableWAIC = TRUE,
                            monitors = parameters, thin = 1, useConjugacy = TRUE,
                            nodes = conf_nodes, monitors2 = parameters2, thin2 = 10)
      
      # Replace samplers
      conf$removeSampler("y.true")
      conf$addSampler(target = paste0("y.true[1:", M, ",1:", data$trap.covariates$J, "]"),
                      type = 'IDSampler',
                      control = list(M = M, J = data$trap.covariates$J, K1D = data$occasion.covariates$K1D,
                                     n.cat = n.cat, n.samples = nimbuild$n.samples, n.ID = n.ID,
                                     this.j = nimbuild$this.j, G.noID = nimbuild$G.noID,
                                     G.latent.ID = nimbuild$G.latent.ID), silent = TRUE)
      
      conf$removeSampler("G.true")
      for (i in 1:M) {
        for (m in 1:n.cat) {
          conf$addSampler(target = paste("G.true[", i, ",", m, "]", sep = ""),
                          type = 'GSampler',
                          control = list(i = i, m = m, M = M, n.cat = n.cat, n.samples = nimbuild$n.samples,
                                         n.levels = n.levels, G.noID = nimbuild$G.noID), silent = TRUE)
        }
      }
      
      conf$removeSampler(paste("s[1:", M, ", 1:2]", sep = ""))
      for (i in 1:M) {
        conf$addSampler(target = paste("s[", i, ", 1:2]", sep = ""),
                        type = 'AF_slice', control = list(adaptive = TRUE), silent = TRUE)
      }
      
      if (grepl("define_model_M1", model_script)) {
        conf$removeSampler(c("lam0.beta0", "sigma.beta0"))
        conf$addSampler(target = c("lam0.beta0", "sigma.beta0"),
                        type = 'AF_slice', control = list(adaptive = TRUE), silent = TRUE)
      }
      if (grepl("define_model_M2", model_script)) {
        conf$removeSampler(c("lam0.beta0", "sigma.beta0", "lam0.beta.sex"))
        conf$addSampler(target = c("lam0.beta0", "sigma.beta0", "lam0.beta.sex"),
                        type = 'AF_slice', control = list(adaptive = TRUE), silent = TRUE)
      }
      if (grepl("define_model_M3", model_script)) {
        conf$removeSampler(c("lam0.beta0", "sigma.beta0", "sigma.beta.sex"))
        conf$addSampler(target = c("lam0.beta0", "sigma.beta0", "sigma.beta.sex"),
                        type = 'AF_slice', control = list(adaptive = TRUE), silent = TRUE)
      }
      if (grepl("define_model_M4", model_script)) {
        conf$removeSampler(c("lam0.beta0", "sigma.beta0", "lam0.beta.sex", "sigma.beta.sex"))
        conf$addSampler(target = c("lam0.beta0", "sigma.beta0", "lam0.beta.sex", "sigma.beta.sex"),
                        type = 'AF_slice', control = list(adaptive = TRUE), silent = TRUE)     
      }
      
      # Build and compile ----
      Rmcmc <- buildMCMC(conf)
      Cmodel <- compileNimble(Rmodel)
      Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
      
      # Run the model ----
      run.start <- Sys.time()
      Cmcmc$run(niter=iterations, nburnin=burnin, reset = TRUE)
      run.end <- Sys.time()
      
      # Get waic and store any text returned ----
      waic_output_text <- capture.output(
        waic_chain <- Cmcmc$getWAIC(),
        type = "output"
      )
      
      # Get summary information ----
      mvSamples <- as.matrix(Cmcmc$mvSamples)
      mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
      
      area <- nimbuild$area
      summary_results <- data.frame(
        total_time = paste(round(run.end - run.start, 2), units(run.end - run.start)),
        area = paste(round(area), "sqkm"), n_traps = data$trap.covariates$J,
        M1 = M1, M2 = M2, M = M, n_ID = n.ID, n_noID = dim(nimbuild$G.noID)[1],
        N_total = mean(mvSamples[-c(1:burnin), "N"]),
        n_spatial_recaptures = sum(apply(apply(data$detections$y.marked, c(1,2), sum), 1, function(x) sum(x>0)) > 1),
        D = round(mean(mvSamples[-c(1:burnin), "N"]) / area, 2)
      )
      
      return(list(mvSamples = mvSamples, mvSamples2 = mvSamples2, time = run.end - run.start,
                  run.seed = seeds[c], stateList = list(modelState = getModelState(Cmodel),
                                                        mcmcState = getMCMCstate(conf, Cmcmc)),
                  summary_results = summary_results, gammaMat = gammaMat,
                  waic_chain = waic_chain, waic_output_text=waic_output_text))
    }
    
    stopCluster(cl.tmp)
    if (verbose) cat("✓ MCMC completed\n")
    
    # Gamma key ----
    gammaMat_translation <- rbind(
      t(
        sapply(
          data$ID.covariates$category_key, 
          function(x) paste(c(as.character(x$value), 
                              rep("irrelevant", ncol(out[[1]]$gammaMat)-length(x$value)))
          )
        )
      )
    )
    # Process chains ----
    n.chain <- length(out)
    stateList <- mvSamplesList <- mvSamples2List <- vector("list", n.chain)
    file_to_append_to <- paste0(output_dir, combo_name, "_summary_results.csv")
    
    for (i in 1:n.chain) {
      if (sbc) {
        mvSamplesList[[i]] <- mcmc(out[[i]]$mvSamples[-1, ])
        mvSamples2List[[i]] <- mcmc(out[[i]]$mvSamples2)
        stateList[[i]] <- out[[i]]$stateList
        mvSamplesList[[i]] <- as.mcmc(mvSamplesList[[i]][complete.cases(mvSamplesList[[i]]), ])
      } else {
        mvSamplesList[[i]] <- mcmc(out[[i]]$mvSamples[-1, ])
        mvSamples2List[[i]] <- mcmc(out[[i]]$mvSamples2)
        stateList[[i]] <- out[[i]]$stateList
        colnames(mvSamplesList[[i]]) <- gsub("gammaMat", "gammaMat_translation", colnames(mvSamplesList[[i]]))
        colnames(mvSamplesList[[i]])[grepl("gammaMat_translation", colnames(mvSamplesList[[i]]))] <- 
          sapply(colnames(mvSamplesList[[i]])[grepl("gammaMat_translation", colnames(mvSamplesList[[i]]))],
                 function(y) eval(parse(text = y)))
        colnames(mvSamplesList[[i]])[4] <- "Bicolor FALSE"
        colnames(mvSamplesList[[i]])[7] <- "Bicolor TRUE"
        if (!grepl("Grid1000m", combo_name)) {
          colnames(mvSamplesList[[i]])[5] <- "Long Hair FALSE"
          colnames(mvSamplesList[[i]])[8] <- "Long Hair TRUE"
        }
        mvSamplesList[[i]] <- as.mcmc(mvSamplesList[[i]][complete.cases(mvSamplesList[[i]]), ])
      }
      if (i == 1) {
        write.csv(cbind(out[[i]]$summary_results, chain = i), file = file_to_append_to, row.names = FALSE)
      } else {
        write.csv(unique(rbind(read.csv(file_to_append_to), cbind(out[[i]]$summary_results, chain = i))),
                  file = file_to_append_to, row.names = FALSE)
      }
    }
    
    # Diagnostics ----
    pdf(paste0(output_dir, combo_name, "_trace.pdf"))
    gelman.plot(mcmc.list(mvSamplesList)[, 1], main = "Gelman convergence - N", ylim = c(0, 2))
    abline(h = 1.1, lty = 2)
    plot(mcmc.list(mvSamplesList))
    dev.off()
    
    # Save outputs ----
    category_key <- data$ID.covariates$category_key
    save(list = c("out", "data", "category_key"), file = paste0(output_dir, combo_name, "_results.Rdata"))
    
    if (verbose) cat("✓ Results saved\n")
    
    return(list(out = out, mvSamplesList = mvSamplesList, mvSamples2List = mvSamples2List,
                summary_file = file_to_append_to, total_time = Sys.time() - start.time))
    
  }, error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}