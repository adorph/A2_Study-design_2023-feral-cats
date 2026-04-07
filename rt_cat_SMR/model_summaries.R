rm(list=ls())
gc()

# ---- Load libraries ---
require(tidyverse)
require(HDInterval)
require(coda)

# ---- Set constants ---
BURNIN = 0
# Removed within model process now

# ---- Create helper function ---
event_stats <- function(df, sex_code) {
  df1 = df %>%
    filter(samp.type == 1, sex == sex_code) %>%
    distinct(Camera, Occasion, ID.name) %>%
    count(ID.name) %>%
    summarise(
      nEvents_ID = n(),
      meanEvents = mean(n),
      sdEvents = sd(n)
    )
  
  df2 = df %>%
    filter(samp.type == 1, sex == sex_code) %>%
    distinct(Camera, ID.name) %>%
    count(ID.name) %>%
    summarise(
      meanCamera = mean(n),
      sdCamera = sd(n)
    )
  return(cbind(df1, df2))
}

# ---- Get model combinations ---
all_files = list.files('outputs/rt_cat_SMR', recursive=TRUE, pattern="Rdata")
all_files = sub("_results.Rdata", "", all_files)
exclude_folder = "old|sbc|test"
filtered_files = all_files[!grepl(exclude_folder, all_files)]

# subset = "INF"
# filtered_files = filtered_files[grepl(subset, filtered_files)]

# ---- Set up empty output vectors ---
gelmanDF    = NULL
estimates   = NULL
inputs      = NULL
ess         = NULL
ess_summary = NULL
cors        = NULL
waic        = NULL
hdi_summary = NULL
# ---- Summarise each model/ array combination ---
for (model in filtered_files) {
  # ---- Load input data ---
  array_name = sub("(.*)/", "", model)
  load(paste0("data_edited/rt_cat_SMR/", array_name, "_data.Rdata"))
  load(paste0("data_edited/rt_cat_SMR/", array_name, "_data_nimbuild.Rdata"))
  # ---- Calculate event statistics ---
  stats_female  = event_stats(data$raw, 1)
  stats_male    = event_stats(data$raw, 2)
  stats_unknown = event_stats(data$raw, 0)
  # ---- Get summary of model inputs ---
  inputs[[model]] = as.data.frame(cbind(
    Area                = nimbuild$area, 
    AugmentedN          = nimbuild$M, 
    nEvents             = sum(data$detections$y.marked.noID) + sum(data$detections$y.marked.noID) + sum(data$detections$y.unmarked) + sum(data$detections$y.unk), 
    nID                 = data$detections$n.ID, 
    nID_female          = length(data$ID.covariates$G.marked[,1][data$ID.covariates$G.marked[,1] == 1]), 
    nID_male            = length(data$ID.covariates$G.marked[,1][data$ID.covariates$G.marked[,1] == 2]), 
    nID_unknown         = length(data$ID.covariates$G.marked[,1][data$ID.covariates$G.marked[,1] == 0]),
    nEvents_IDfemale    = nrow(data$raw %>% filter(samp.type==1, sex==1) %>% select(Camera, Occasion, ID.name) %>% distinct()), 
    nEvents_perF_mean   = stats_female$meanEvents, 
    nEvents_perF_sd     = stats_female$sdEvents, 
    nCameras_perF_mean  = stats_female$meanCamera, 
    nCameras_perF_sd    = stats_female$sdCamera,
    nEvents_IDmale      = nrow(data$raw %>% filter(samp.type==1, sex==2) %>% select(Camera, Occasion, ID.name) %>% distinct()), 
    nEvents_perM_mean   = stats_male$meanEvents, 
    nEvents_perM_sd     = stats_male$sdEvents, 
    nCameras_perM_mean  = stats_male$meanCamera, 
    nCameras_perM_sd    = stats_male$sdCamera,
    nEvents_IDunknown   = nrow(data$raw %>% filter(samp.type==1, sex==0) %>% select(Camera, Occasion, ID.name) %>% distinct()), 
    nEvents_perU_mean   = stats_unknown$meanEvents, 
    nEvents_perU_sd     = stats_unknown$sdEvents, 
    nCameras_perU_mean  = stats_unknown$meanCamera, 
    nCameras_perU_sd    = stats_unknown$sdCamera,
    nEvents_noID        = sum(data$detections$y.marked.noID), 
    nEvents_unmarked    = sum(data$detections$y.unmarked), 
    nEvents_unknown     = sum(data$detections$y.unk)
    )
  )
  # ---- Load model output ---
  load(paste0("outputs/rt_cat_SMR/", model, "_results.Rdata"))  
  # ---- Model chains ---
  gammaMat_translation = t(
    sapply(data$ID.covariates$category_key, 
           function(x) paste(c(as.character(x$value), rep("irrelevant", ncol(out[[1]]$gammaMat)-length(x$value))))
           )
    )
  mvSamples  = as.matrix(out[[1]]$mvSamples)
  mvSamples2 = as.matrix(out[[1]]$mvSamples2)
  # ---- Set up empty vectors ---
  n_chain        = length(out)
  new_samples    = vector("list", n_chain)
  ests           = vector("list", n_chain)
  cor            = vector("list", n_chain)
  waic_stats     = vector("list", n_chain)
  stateList      = vector("list", n_chain)
  mvSamplesList  = vector("list", n_chain)
  mvSamples2List = vector("list", n_chain)
  HDIoverlap     = vector("list", n_chain)
  # ---- Summarise each chain ---
  for (i in 1:n_chain) {
  # ---- Extract samples and rename columns to covariate categories ---
    mvSamplesList[[i]]  = mcmc(out[[i]]$mvSamples[-1,]) #remove first iteration bc screws up plotting
    mvSamples2List[[i]] = mcmc(out[[i]]$mvSamples2)
    stateList[[i]]      = out[[i]]$stateList
    # ---- Rename columns: If 1000m grid array do not include long hair summary
    if (grepl("Grid1000m", array_name)) {
      colnames(mvSamplesList[[i]]) = gsub("gammaMat", "gammaMat_translation", colnames(mvSamplesList[[i]]))
      colnames(mvSamplesList[[i]])[grepl("gammaMat_translation", colnames(mvSamplesList[[i]]))] = sapply(colnames(mvSamplesList[[i]])[grepl("gammaMat_translation",colnames(mvSamplesList[[i]]))], function(y) eval(parse(text = y)))
      colnames(mvSamplesList[[i]])[4] = "Bicolor FALSE"
      colnames(mvSamplesList[[i]])[7] = "Bicolor TRUE"
    } else {
      colnames(mvSamplesList[[i]]) = gsub("gammaMat", "gammaMat_translation", colnames(mvSamplesList[[i]]))
      colnames(mvSamplesList[[i]])[grepl("gammaMat_translation",colnames(mvSamplesList[[i]]))] = sapply(colnames(mvSamplesList[[i]])[grepl("gammaMat_translation", colnames(mvSamplesList[[i]]))], function(y) eval(parse(text = y)))
      colnames(mvSamplesList[[i]])[4] = "Bicolor FALSE"
      colnames(mvSamplesList[[i]])[5] = "Long Hair FALSE"
      colnames(mvSamplesList[[i]])[8] = "Bicolor TRUE"
      colnames(mvSamplesList[[i]])[9] = "Long Hair TRUE"
    }
    mvSamplesList[[i]] = as.mcmc(mvSamplesList[[i]][complete.cases(mvSamplesList[[i]]),])
    new_samples[[i]]   = mvSamplesList[[i]][, !grepl("irrelevant", colnames(mvSamplesList[[i]]))]
    # ---- Extract parameter estimates based on model type and removing the burnin ---
    N   = mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'N']
    est = mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]), "N"]/(nimbuild$area)
    if (grepl("M1", model)){
      sigma    = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta0'])
      g0       = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta0'])
      est_samp = as.data.frame(cbind(N, est, sigma, g0))
    } else if (grepl("M2", model)) {
      sigma    = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta0'])
      g0_f     = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta0'])
      g0_m     = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta0'] + mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta.sex'] )
      est_samp = as.data.frame(cbind(N, est, sigma, g0_f, g0_m))
    } else if (grepl("M3", model)) {
      sigma_f  = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta0'])
      sigma_m  = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta0'] + mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta.sex'] )
      g0       = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta0'])
      est_samp = as.data.frame(cbind(N, est, sigma_f, sigma_m, g0))
    } else {
      sigma_f  = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta0'])
      g0_f     = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta0'])
      sigma_m  = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta0'] + mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta.sex'] )
      g0_m     = exp(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta0'] + mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta.sex'] )
      est_samp = as.data.frame(cbind(N, est, sigma_f, sigma_m, g0_f, g0_m))
    }
    est_temp = est_samp %>% summarise_all(.funs=c('mean', 'median', 'sd'), na.rm = TRUE)
    est_temp2 = est_samp %>%
      summarise_all(hdi) %>%
      mutate(names = c("lcl", "ucl")) %>%
      pivot_wider(id_cols=NULL, names_from = names, values_from = everything())
    ests[[i]] = as.data.frame(cbind(est_temp, est_temp2)) %>% 
      select(starts_with("N"), starts_with("est"), starts_with("sigma"), starts_with("g0")) %>%
      select(-starts_with("names"))
    # ---- Check 95% intervals of posterior sex effects overlap zero ---
    if (grepl("M4", model)) {
      sigma_overlaps_0 = if ( quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta.sex'], 0.025) < 0 & quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta.sex'], 0.975) > 0 ) {
        TRUE  
      } else {FALSE}
      g0_overlaps_0 = if ( quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta.sex'], 0.025) < 0 & quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta.sex'], 0.975) > 0 ) {
        TRUE  
      } else {FALSE}
      HDIoverlap[[i]] = cbind(sigma_overlaps_0, g0_overlaps_0)
    } else if (grepl("M3", model)) {
      sigma_overlaps_0 = if ( quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta.sex'], 0.025) < 0 & quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'sigma.beta.sex'], 0.975) > 0 ) {
        TRUE  
      } else {FALSE}
      HDIoverlap[[i]] = cbind(sigma_overlaps_0, NA)
    } else if (grepl("M2", model)) {
      g0_overlaps_0 = if ( quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta.sex'], 0.025) < 0 & quantile(mvSamplesList[[i]][BURNIN:nrow(mvSamplesList[[i]]),'lam0.beta.sex'], 0.975) > 0 ) {
        TRUE  
      } else {FALSE}
      HDIoverlap[[i]] = cbind(NA, g0_overlaps_0)
    } else {
      HDIoverlap[[i]] = cbind(NA, NA)
    }
    # ---- Calculate correlations between sigma and g0 for males and females ---
    if (grepl("M4", model)) {
      cor[[i]] <- cbind (female = cor(est_samp$sigma_f, est_samp$g0_f, method="pearson"),
                         male   = cor(est_samp$sigma_m, est_samp$g0_m, method="pearson") )
    }
    # ---- Extract WAIC values ---
    waic_stats[[i]] <- as.data.frame(cbind(out[[i]]$waic_chain$lppd, out[[i]]$waic_chain$pWAIC, out[[i]]$waic_chain$WAIC, out[[i]]$waic_output_text))
  }
  # ---- Calculate effective sample size ---
  total_iterations = nrow(new_samples[[1]]) * length(new_samples)
  ess_all = effectiveSize(mcmc.list(new_samples))
  ess_ratio = ess_all / total_iterations
  summary_ess = data.frame(
    parameter = names(ess_all),
    ESS = as.numeric(ess_all),
    total_samples = total_iterations,
    ESS_ratio = as.numeric(ess_all) / total_iterations
  )  
  # ---- Convergence ---
  gd2 <- data.frame(variable = NA_character_, est = NA_real_, uci = NA_real_)
  for(s in 1:ncol(new_samples[[1]])) {
    gd      = gelman.diag(mcmc.list(new_samples)[,s])
    gd2[s,] = cbind(variable = colnames(new_samples[[1]])[s], est = gd$psrf[,1], uci = gd$psrf[,2])
  }
  # ---- Bind together results ---
  estimates[[model]] = colSums(bind_rows(ests)) / 5
  cors[[model]] = cor
  waic[[model]] = bind_rows(waic_stats) 
  ess[[model]]  = summary_ess
  ess_summary[[model]] = data.frame(
    mean_ESS = mean(summary_ess$ESS),
    min_ESS =  min(summary_ess$ESS),
    mean_ESS_ratio = mean(summary_ess$ESS_ratio),
    min_ESS_ratio = min(summary_ess$ESS_ratio)
  )
  hdi_summary[[model]] = HDIoverlap
  gelmanDF[[model]] <- gd2
  # ---- Tidy up ---
  all_objects = ls()
  object_to_keep = c('BURNIN', 'event_stats', 'filtered_files', 'inputs', 'estimates', 'cors', 'waic', 'ess', 'ess_summary', 'hdi_summary', 'gelmanDF')
  rm(list = setdiff(all_objects, object_to_keep))
}

# ---- Get key estimates --- 
est_main = bind_rows(estimates, .id = "list_name") %>%
  select(list_name, est_median, est_sd, est_lcl, est_ucl, 
         sigma_median, sigma_sd, sigma_lcl, sigma_ucl, 
         g0_median, g0_sd, g0_ucl, g0_lcl,
         sigma_f_median, sigma_f_sd, sigma_f_lcl, sigma_f_ucl, 
         sigma_m_median, sigma_m_sd, sigma_m_lcl, sigma_m_ucl, 
         g0_f_median, g0_f_sd, g0_f_lcl, g0_f_ucl,
         g0_m_median, g0_m_sd, g0_m_lcl, g0_m_ucl)

write.csv(est_main, "./outputs/nimble_summary_estimates.csv")

# ---- Check out large Gelman values ---
bind_rows(gelmanDF, .id = "list_name") %>%
  filter(uci > 1.1) %>%
  arrange(list_name)
# ---- Get ESS ---
bind_rows(ess_summary, .id = "list_name")
# ---- Get Correlations ---
cor_means <- t(sapply(cors, function(outer) {
  female <- mean(sapply(outer, function(x) x[, "female"]), na.rm=T)
  male   <- mean(sapply(outer, function(x) x[, "male"]), na.rm=T)
  c(female = female, male = male)
}))
# ---- Get WAIC ---
waic_bind = bind_rows(waic, .id = "list_name") %>%
  group_by(list_name) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  select(list_name, V3) %>%
  mutate(array = sub("(.*)/", "", list_name),
         WAIC = as.numeric(V3)) %>%
  group_by(array) %>%
  mutate(deltaWAIC = WAIC - min(WAIC)) #experimental


bind_rows(lapply(hdi_summary, function(x) as.data.frame(t(x))), .id = "list_name")

ins = bind_rows(inputs, .id = "list_name")
