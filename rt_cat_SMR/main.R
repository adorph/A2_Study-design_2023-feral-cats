rm(list=ls())
gc()

require(tidyverse)
source("code/rt_cat_SMR/prepare_camera_data.R")
source("code/rt_cat_SMR/prepare_nimble_data.R")
source("code/rt_cat_SMR/run_nimble_model.R")

# ===== DEFINE ALL CAMERA LOCATION FILES ======================================
camera_files <- list(
  Grid250m = "./data_raw/PL_0250m_Grid_CameraLocations.csv",
  Grid500m = "./data_raw/PL_0500m_Grid_CameraLocations.csv",
  Grid1000m = "./data_raw/PL_1000m_Grid_CameraLocations.csv",
  Road = "./data_raw/PL_Predator_CameraLocations.csv"
)

# Add in subsetted cameras
sub_camera_files <- list(
  Grid250m_sparseA = "./data_raw/PL_0250m_Grid_sparseA_CameraLocations.csv",
  Grid250m_sparseB = "./data_raw/PL_0250m_Grid_sparseB_CameraLocations.csv",
  Grid500m_sparseA = "./data_raw/PL_0500m_Grid_sparseA_CameraLocations.csv",
  Grid500m_sparseB = "./data_raw/PL_0500m_Grid_sparseB_CameraLocations.csv",
  Grid1000m_areaA = "./data_raw/PL_1000m_Grid_areaA_CameraLocations.csv",
  Grid1000m_areaB = "./data_raw/PL_1000m_Grid_areaB_CameraLocations.csv"
)

# ===== LOAD ALL CAMERA DATA WITH STANDARDISATION =============================
cameras <- lapply(camera_files, function(f) {
  read.csv(f) %>%
    rename(
      CameraCode = if_else("CameraCode" %in% names(.), "CameraCode", "name"),
      Latitude = if_else("Latitude" %in% names(.), "Latitude", names(.)[grep("Lat", names(.), ignore.case = TRUE)][1]),
      Longitude = if_else("Longitude" %in% names(.), "Longitude", names(.)[grep("Lon", names(.), ignore.case = TRUE)][1])
    ) %>%
    dplyr::select(CameraCode, Latitude, Longitude)
})

sub_cameras  <- lapply(sub_camera_files, function(f) {
  read.csv(f) %>%
    rename(
      CameraCode = if_else("CameraCode" %in% names(.), "CameraCode", "name"),
      Latitude = if_else("Latitude" %in% names(.), "Latitude", names(.)[grep("Lat", names(.), ignore.case = TRUE)][1]),
      Longitude = if_else("Longitude" %in% names(.), "Longitude", names(.)[grep("Lon", names(.), ignore.case = TRUE)][1])
    ) %>%
    dplyr::select(CameraCode, Latitude, Longitude)
})

# ===== GENERATE ALL COMBINATIONS + SINGLE ARRAYS ==============================
array_names <- names(cameras)
all_combinations <- list()
combination_names <- list()
combination_files <- list()

# Add single arrays first
for (arr in array_names) {
  combo_name <- arr
  all_combinations[[combo_name]] <- cameras[[arr]]
  combination_names[[combo_name]] <- arr
  
  output_file <- paste0("./data_raw/PL_", combo_name, "_locations.csv")
  write.csv(all_combinations[[combo_name]], output_file, row.names = FALSE)
  combination_files[[combo_name]] <- output_file
  
  cat(combo_name, ":", nrow(all_combinations[[combo_name]]), "cameras\n")
}

sub_array_names <- names(sub_cameras)
sub_all_combinations <- list()
sub_combination_names <- list()
sub_combination_files <- list()

# Now subset arrays
for (arr in sub_array_names) {
  combo_name <- arr
  sub_all_combinations[[combo_name]] <- sub_cameras[[arr]]
  sub_combination_names[[combo_name]] <- arr
  
  output_file <- paste0("./data_raw/PL_", combo_name, "_locations.csv")
  write.csv(sub_all_combinations[[combo_name]], output_file, row.names = FALSE)
  sub_combination_files[[combo_name]] <- output_file
  
  cat(combo_name, ":", nrow(sub_all_combinations[[combo_name]]), "cameras\n")
}

# ===== PREPARE ALL CAMERA DATA ================================================
cat_files <- c(
  "./data_raw/2023_PL_Grid250m_CAT_cleaned.csv",
  "./data_raw/2023_PL_Grid500m_CAT_cleaned.csv",
  "./data_raw/2023_PL_Grid1000m_CAT_cleaned.csv",
  "./data_raw/2023_PL_Predator_CAT_cleaned.csv"
)

op_files <- c(
  "./data_raw/2023_PL_S01_Grid250m_OT.csv",
  "./data_raw/2023_PL_S02_Grid250m_OT.csv",
  "./data_raw/2023_PL_S01_Grid500m_OT.csv",
  "./data_raw/2023_PL_S02_Grid500m_OT.csv",
  "./data_raw/2023_PL_S01_Grid1km_OT.csv",
  "./data_raw/2023_PL_S02_Grid1km_OT.csv",
  "./data_raw/2023_PL_S03_Predator_OT.csv",
  "./data_raw/2023_PL_S04_Predator_OT.csv"
)

# Loop through informed data combinations
results <- mapply(function(combo_name, combo_arrays, combo_file) {
  cat("\nProcessing:", combo_name, 
      "(", paste(combo_arrays, collapse = " + "), ")\n")
  
  prepare_camera_data(
    array_type = combo_arrays,
    camera_locations_file = combo_file,
    cat_data_files = cat_files,
    cat_data_type = "Informed",
    operation_tracking_files = op_files,
    include_longer_hair = TRUE,
    output_path = "./data_edited/rt_cat_SMR",
    output_filename = paste0("INF_", combo_name, "_data.Rdata")
  )
}, 
names(all_combinations), 
combination_names,
combination_files,
SIMPLIFY = FALSE)

# Loop through uninformed data combinations
uninformed_names <- c("Grid250m", "Grid500m", "Grid1000m", "Predator")
results <- mapply(function(combo_name, combo_arrays, combo_file, orig_name) {
  cat("\nProcessing:", combo_name, 
      "(", paste(combo_arrays, collapse = " + "), ")\n")
  
  prepare_camera_data(
    array_type = combo_arrays,
    camera_locations_file = combo_file,
    cat_data_files = paste0("./data_raw/2023_PL_", orig_name, "_CAT_cleaned.csv"),
    cat_data_type = "Uninformed",
    operation_tracking_files = op_files,
    include_longer_hair = TRUE,
    output_path = "./data_edited/rt_cat_SMR",
    output_filename = paste0("UNI_", combo_name, "_data.Rdata")
  )
}, 
names(all_combinations), 
combination_names,
combination_files,
uninformed_names,
SIMPLIFY = FALSE)

# Loop through subsets of informed data combinations
results <- mapply(function(combo_name, combo_arrays, combo_file) {
  cat("\nProcessing:", combo_name, 
      "(", paste(combo_arrays, collapse = " + "), ")\n")
  
  prepare_camera_data(
    array_type = combo_arrays,
    camera_locations_file = combo_file,
    cat_data_files = cat_files,
    cat_data_type = "Informed",
    operation_tracking_files = op_files,
    include_longer_hair = TRUE,
    output_path = "./data_edited/rt_cat_SMR",
    output_filename = paste0("INF_", combo_name, "_data.Rdata")
  )
}, 
names(sub_all_combinations), 
sub_combination_names,
sub_combination_files,
SIMPLIFY = FALSE)

# single array
# prepare_camera_data(
#   array_type = "Road",
#   camera_locations_file = "./data_raw/PL_Road_locations.csv",
#   cat_data_files = cat_files,
#   cat_data_type = "Uninformed",
#   operation_tracking_files = op_files,
#   include_longer_hair = TRUE,
#   output_path = "./data_edited/rt_cat_SMR",
#   output_filename = paste0("UNI_Road_data.Rdata")
# )

# ===== PREPARE ALL COMBINATIONS FOR USE IN NIMBLE =============================
# Get all files names
file_list <- c(paste0("INF_", array_names), paste0("UNI_", array_names))
file_list <- paste0("INF_", sub_array_names)
nimble_results <- lapply(file_list, function(combo_name) {
  cat("\n\n===== PREPARING NIMBLE DATA:", combo_name, "=====\n")
  
  data_file <- paste0("./data_edited/rt_cat_SMR/", combo_name, "_data.Rdata")
  
  if (!file.exists(data_file)) {
    cat("✗ Data file not found:", data_file, "\n")
    return(NULL)
  }
  
  prepare_nimble_data(
    data_file = data_file,
    output_dir = "./data_edited/rt_cat_SMR",
    M1_factor = 1.3,
    M2_factor = 1,
    verbose = TRUE
  )
})

names(nimble_results) <- file_list

# Check results generated
for (name in names(nimble_results)) {
  if (!is.null(nimble_results[[name]])) {
    cat(name, "- M =", nimble_results[[name]]$M, 
        ", n.cells =", nimble_results[[name]]$n.cells, "\n")
  } else {
    cat(name, "- FAILED\n")
  }
}

# Single array
# prepare_nimble_data(
#   data_file = "./data_edited/rt_cat_SMR/UNI_Road_data.Rdata",
#   output_dir = "./data_edited/rt_cat_SMR",
#   M1_factor = 1.3,
#   M2_factor = 1,
#   verbose = TRUE
# )

# Check nimble data generated ok
# Check generation of co-variates in data - make sure all are represented
# If detections are incorrect there is likely an issue with the creation of your data (e.g. cameras ordered incorrectly)
for(file in file_list) {
  
  cat("\n\n===== CHECKING NIMBLE DATA:", file, "=====\n")
  
  data_file <- paste0("./data_edited/rt_cat_SMR/", file, "_data.Rdata")
  nimble_file <- paste0("./data_edited/rt_cat_SMR/", file, "_data_nimbuild.Rdata")
  
  if (!file.exists(data_file)) {
    cat("✗ Data file not found:", data_file, "\n")
    return(NULL)
  }
  
  load(data_file)
  load(nimble_file)

  covariates = list(
    sex = round(cbind(
      female   = table(nimbuild$G.true[,1])[1] / sum(table(nimbuild$G.true[,1])),
      male     = table(nimbuild$G.true[,1])[2] / sum(table(nimbuild$G.true[,1]))
    ), 4), 
    Primary_pattern = round(cbind(
      black    = table(nimbuild$G.true[,2])[1] / sum(table(nimbuild$G.true[,2])),
      ginger   = table(nimbuild$G.true[,2])[2] / sum(table(nimbuild$G.true[,2])),
      classic  = table(nimbuild$G.true[,2])[3] / sum(table(nimbuild$G.true[,2])),
      mackerel = table(nimbuild$G.true[,2])[4] / sum(table(nimbuild$G.true[,2])),
      tortoise = table(nimbuild$G.true[,2])[5] / sum(table(nimbuild$G.true[,2]))
    ), 4),
    Bicolor = round(cbind(
      biFALSE  = table(nimbuild$G.true[,3])[1] / sum(table(nimbuild$G.true[,3])),
      biTRUE   = table(nimbuild$G.true[,3])[2] / sum(table(nimbuild$G.true[,3]))
    ), 4),
    Longer_hair = round(cbind(
      longFALSE = table(nimbuild$G.true[,4])[1] / sum(table(nimbuild$G.true[,4])),
      longTRUE  = table(nimbuild$G.true[,4])[2] / sum(table(nimbuild$G.true[,4]))
    ), 4)
  )  
  print(covariates)
  
  X <- data$trap.covariates$X 
  tmp2 <- as.data.frame(data$trap.covariates$X)
  
  for(i in 1:data$detections$n.ID) {
    idtmp <- (X[nimbuild$y.true[i,]>0,1:2])
    if (length(idtmp) <= 2 ) {
      idtmp <- t(as.data.frame(idtmp))
    }
    gg <- ggplot(tmp2, aes(x=X, y=Y)) +
      geom_point() +
      geom_point(idtmp, mapping=aes(x=X, y=Y), color="red", size=3) + 
      ggtitle(paste("Data:", file, "\nIndividual:", i))
    print(gg)
  }
}

# To visualise a specific mask area:
# image(matrix(nimble_results$Grid1000m$InSS,nimble_results$Grid1000m$n.cells.x,nimble_results$Grid1000m$n.cells.y),main="Mask Area")


# ===== PROCESS ALL COMBINATIONS THROUGH NIMBLE ================================

# Uninformed single arrays:
# Model mapping 
uninformed_mapping <- list(
  UNI_Grid250m_M1 = "define_model_M1A.R",
  UNI_Grid500m_M1 = "define_model_M1A.R",
  UNI_Grid1000m_M1 = "define_model_M1B.R",
  UNI_Road_M1 = "define_model_M1A.R",
  UNI_Grid250m_M2 = "define_model_M2A.R",
  UNI_Grid500m_M2 = "define_model_M2A.R",
  UNI_Grid1000m_M2 = "define_model_M2B.R",
  UNI_Road_M2 = "define_model_M2A.R",
  UNI_Grid250m_M3 = "define_model_M3A.R",
  UNI_Grid500m_M3 = "define_model_M3A.R",
  UNI_Grid1000m_M3 = "define_model_M3B.R",
  UNI_Road_M3 = "define_model_M3A.R",
  UNI_Grid250m_M4 = "define_model_M4A.R",
  UNI_Grid500m_M4 = "define_model_M4A.R",
  UNI_Grid1000m_M4 = "define_model_M4B.R",
  UNI_Road_M4 = "define_model_M4A.R"
)

# Run all
uninformed_results <- lapply(names(uninformed_mapping), function(combo_name) {
  cat("\n===== PROCESSING:", combo_name, "=====\n")

  file_name <- sub("_M\\d+", "", combo_name)
  
  data_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data.Rdata")
  nimbuild_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data_nimbuild.Rdata")
  
  if (!file.exists(data_file) || !file.exists(nimbuild_file)) {
    cat("✗ Files not found\n"); return(NULL)
  }
  
  run_nimble_model(
    data_file = data_file, nimbuild_file = nimbuild_file, combo_name = file_name,
    chains = 5, iterations = 30000, burnin = 10000,
    model_script = uninformed_mapping[[combo_name]], 
    output_dir = paste0("./outputs/rt_cat_SMR/", sub(".*_(M\\d+).*", "\\1", uninformed_mapping[[combo_name]]), "/"), verbose = TRUE
  )
})

names(uninformed_results) <- names(uninformed_mapping)


# Informed single arrays:
# Model mapping
informed_mapping <- list(
  INF_Grid250m_M1 = "define_model_M1A.R",
  INF_Grid500m_M1 = "define_model_M1A.R",
  INF_Grid1000m_M1 = "define_model_M1B.R",
  INF_Road_M1 = "define_model_M1A.R",
  INF_Grid250m_M2 = "define_model_M2A.R",
  INF_Grid500m_M2 = "define_model_M2A.R",
  INF_Grid1000m_M2 = "define_model_M2B.R",
  INF_Road_M2 = "define_model_M2A.R",
  INF_Grid250m_M3 = "define_model_M3A.R",
  INF_Grid500m_M3 = "define_model_M3A.R",
  INF_Grid1000m_M3 = "define_model_M3B.R",
  INF_Road_M3 = "define_model_M3A.R",
  INF_Grid250m_M4 = "define_model_M4A.R",
  INF_Grid500m_M4 = "define_model_M4A.R",
  INF_Grid1000m_M4 = "define_model_M4B.R",
  INF_Road_M4 = "define_model_M4A.R"
)

# Run all
informed_results <- lapply(names(informed_mapping), function(combo_name) {
  cat("\n===== PROCESSING:", combo_name, "=====\n")
  
  file_name <- sub("_M\\d+", "", combo_name)
  
  data_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data.Rdata")
  nimbuild_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data_nimbuild.Rdata")
  
  if (!file.exists(data_file) || !file.exists(nimbuild_file)) {
    cat("✗ Files not found\n"); return(NULL)
  }
  
  run_nimble_model(
    data_file = data_file, nimbuild_file = nimbuild_file, combo_name = file_name,
    chains = 5, iterations = 30000, burnin = 10000,
    model_script = informed_mapping[[combo_name]], 
    output_dir = paste0("./outputs/rt_cat_SMR/", sub(".*_(M\\d+).*", "\\1", informed_mapping[[combo_name]]), "/"), verbose = TRUE
  )
})

names(informed_results) <- names(informed_mapping)

# Subset Informed single arrays:
# Model mapping
subset_mapping <- list(
  INF_Grid250m_sparseA_M1 = "define_model_M1A.R",
  INF_Grid250m_sparseB_M1 = "define_model_M1A.R",
  INF_Grid500m_sparseA_M1 = "define_model_M1A.R",
  INF_Grid500m_sparseB_M1 = "define_model_M1A.R",
  INF_Grid1000m_areaA_M1 = "define_model_M1B.R",
  INF_Grid1000m_areaB_M1 = "define_model_M1B.R",
  INF_Grid250m_sparseA_M2 = "define_model_M2A.R",
  INF_Grid250m_sparseB_M2 = "define_model_M2A.R",
  INF_Grid500m_sparseA_M2 = "define_model_M2A.R",
  INF_Grid500m_sparseB_M2 = "define_model_M2A.R",
  INF_Grid1000m_areaA_M2 = "define_model_M2B.R",
  INF_Grid1000m_areaB_M2 = "define_model_M2B.R",
  INF_Grid250m_sparseA_M3 = "define_model_M3A.R",
  INF_Grid250m_sparseB_M3 = "define_model_M3A.R",
  INF_Grid500m_sparseA_M3 = "define_model_M3A.R",
  INF_Grid500m_sparseB_M3 = "define_model_M3A.R",
  INF_Grid1000m_areaA_M3 = "define_model_M3B.R",
  INF_Grid1000m_areaB_M3 = "define_model_M3B.R",
  INF_Grid250m_sparseA_M4 = "define_model_M4A.R",
  INF_Grid250m_sparseB_M4 = "define_model_M4A.R",
  INF_Grid500m_sparseA_M4 = "define_model_M4A.R",
  INF_Grid500m_sparseB_M4 = "define_model_M4A.R",
  INF_Grid1000m_areaA_M4 = "define_model_M4B.R",
  INF_Grid1000m_areaB_M4 = "define_model_M4B.R"
)

# Run all
subset_results <- lapply(names(subset_mapping), function(combo_name) {
  cat("\n===== PROCESSING:", combo_name, "=====\n")
  
  file_name <- sub("_M\\d+", "", combo_name)
  
  data_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data.Rdata")
  nimbuild_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data_nimbuild.Rdata")
  
  if (!file.exists(data_file) || !file.exists(nimbuild_file)) {
    cat("✗ Files not found\n"); return(NULL)
  }
  
  run_nimble_model(
    data_file = data_file, nimbuild_file = nimbuild_file, combo_name = file_name,
    chains = 5, iterations = 30000, burnin = 10000,
    model_script = subset_mapping[[combo_name]], 
    output_dir = paste0("./outputs/rt_cat_SMR/", sub(".*_(M\\d+).*", "\\1", subset_mapping[[combo_name]]), "/"), verbose = TRUE
  )
})

names(subset_results) <- names(subset_mapping)


# Single array
# single_map <- list(
#   UNI_Road_M1 = "define_model_M1A.R",
#   UNI_Road_M2 = "define_model_M2A.R",
#   UNI_Road_M3 = "define_model_M3A.R",
#   UNI_Road_M4 = "define_model_M4A.R"
# )
# 
# subset_results <- lapply(names(single_map), function(combo_name) {
#   cat("\n===== PROCESSING:", combo_name, "=====\n")
#   
#   file_name <- sub("_M\\d+", "", combo_name)
#   
#   data_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data.Rdata")
#   nimbuild_file <- paste0("data_edited/rt_cat_SMR/", file_name, "_data_nimbuild.Rdata")
#   
#   if (!file.exists(data_file) || !file.exists(nimbuild_file)) {
#     cat("✗ Files not found\n"); return(NULL)
#   }
#   
#   run_nimble_model(
#     data_file = data_file, 
#     nimbuild_file = nimbuild_file, 
#     combo_name = file_name,
#     chains = 5, iterations = 30000, burnin = 10000,
#     model_script = single_map[[combo_name]], 
#     output_dir = paste0("./outputs/rt_cat_SMR/", sub(".*_(M\\d+).*", "\\1", single_map[[combo_name]]), "/"), 
#     verbose = TRUE
#   )
# })
