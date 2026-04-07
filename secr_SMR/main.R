rm(list=ls())
gc()

require(tidyverse)
require(secr)
source("code/secr_SMR/prepare_camera_data.R")
source("code/secr_SMR/run_secr_model.R")

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
    output_path = "./data_edited/secr_SMR/",
    output_filename = paste0("INF_", combo_name)
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
    output_path = "./data_edited/secr_SMR/",
    output_filename = paste0("UNI_", combo_name)
  )
}, 
names(all_combinations)[1:4], 
combination_names[1:4],
combination_files[1:4],
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
    output_path = "./data_edited/secr_SMR/",
    output_filename = paste0("INF_", combo_name)
  )
}, 
names(sub_all_combinations), 
sub_combination_names,
sub_combination_files,
SIMPLIFY = FALSE)


# ===== PROCESS ALL COMBINATIONS THROUGH SECR ================================

file_list <- c(paste0("INF_", array_names), paste0("UNI_", array_names), paste0("INF_", sub_array_names))

results <- lapply(file_list, function(file_name) {
  cat("\n===== PROCESSING:", file_name, "=====\n")
  
  captpath <- paste0('./data_edited/secr_SMR/', file_name, '_marked_cats.txt')
  Tupath  <- paste0('./data_edited/secr_SMR/', file_name, '_unmarked_cats.txt')
  Tmpath  <- paste0('./data_edited/secr_SMR/', file_name, '_uncertain.txt')
  trappath <- paste0('./data_edited/secr_SMR/', file_name, '_detectors.txt')
  oppath <-  paste0('./data_edited/secr_SMR/', file_name, '_operation.txt')
  
  if (!file.exists(captpath) || !file.exists(Tupath) || !file.exists(Tmpath) || !file.exists(trappath) || !file.exists(oppath)) {
    cat("✗ Files not found\n"); return(NULL)
  }
  
  run_secr_model(
    data_file = captpath, 
    Tu_file = Tupath, 
    Tm_file = Tmpath, 
    trap_file = trappath, 
    operation_file = oppath,
    combo_name = file_name,
    output_dir = "./outputs/secr_SMR/", 
    verbose = TRUE
  )
})

# Run single model 
file_name <- "UNI_Road"
captpath <- paste0('./data_edited/secr_SMR/', file_name, '_marked_cats.txt')
Tupath  <- paste0('./data_edited/secr_SMR/', file_name, '_unmarked_cats.txt')
Tmpath  <- paste0('./data_edited/secr_SMR/', file_name, '_uncertain.txt')
trappath <- paste0('./data_edited/secr_SMR/', file_name, '_detectors.txt')
oppath <-  paste0('./data_edited/secr_SMR/', file_name, '_operation.txt')

if (!file.exists(captpath) || !file.exists(Tupath) || !file.exists(Tmpath) || !file.exists(trappath) || !file.exists(oppath)) {
  cat("✗ Files not found\n"); return(NULL)
}
run_secr_model(
  data_file = captpath, 
  Tu_file = Tupath, 
  Tm_file = Tmpath, 
  trap_file = trappath, 
  operation_file = oppath,
  combo_name = file_name,
  output_dir = "./outputs/secr_SMR/", 
  verbose = TRUE
)


# ===== CHECK MODELS ================================
model_file <- list.files("./outputs/secr_SMR/", full.names =T)

for (i in 1:length(model_file)) {
  file_name <- model_file[i]
  load(file_name)
  
  ids <- rownames(catCH2)
  cols <- Polychrome::createPalette(N = 19, seedcolors = c("#E69F00", "#56B4E9"))
  
  # plot using only secr functions (tracks = TRUE draws spatial connections)
  plot(catCH2,
       tracks   = TRUE,    # join consecutive captures with lines
       varycol  = TRUE,    # use per-individual colours
       icolours = cols,    # supply our colours
       rad      = 0.05,    # small radial jitter to separate overlapping points
       hidetraps = FALSE)  # show traps (set TRUE to hide)
  legend("topright",
         legend = ids,
         col    = cols,
         pch    = 19,
         cex    = 0.7,
         title  = "Individuals",
         bg     = "white",
         box.lty = 0)
  
  secr::esaPlot(fit_M3)
  # i=i+1
}



# ===== EXTRACT ESTIMATES FROM MODELS ================================

# Define all dataset combinations to run (must match previous loop output names)
rm(densityest)
for (file_name in model_file) {
  load(file_name)
  
  array_name = sub("\\.Rdata$", "", sub(".*/[0-9]+_", "", file_name))
  
  if(!exists("densityest")) {  
    densityest <- cbind(array=array_name, 
                        print(summary(fit_M3)$predicted[1, ]),
                        print(summary(fit_M3)$predicted[2, ]),
                        print(summary(fit_M3)$predicted[3, ]))
  } else {
    densityest <- rbind(densityest, 
                        cbind(array=array_name, 
                              print(summary(fit_M3)$predicted[1, ]),
                              print(summary(fit_M3)$predicted[2, ]),
                              print(summary(fit_M3)$predicted[3, ])))
  }
}
colnames(densityest) <- c("array", "link", "estimate", "SE.estimate", "lcl", "ucl",
                          "g0.link", "g0.estimate", "g0.SE.estimate", "g0.lcl", "g0.ucl",
                          "sigma.link", "sigma.estimate", "sigma.SE.estimate", "sigma.lcl", "sigma.ucl" )
densityest$estimate <- densityest$estimate*100
densityest$SE.estimate <- densityest$SE.estimate*100
densityest$lcl <- densityest$lcl * 100
densityest$ucl <- densityest$ucl * 100

# Remove first road estimate (wrong buffer size)
densityest <- densityest[-14,]

write.csv(densityest, "./outputs/secr_summary_estimates.csv")


# ===== TABLE CAPTURE SUMMARIES ================================

sumres <- NULL
for(file_name in model_file) {
  # Load model - informed
  load(file_name)
  # Get capthist
  caps <- apply(catCH2, 1, sum)
  # Capture matrix: 1 if animal i was detected at trap j at any occasion
  ID_matrix <- apply(catCH2, c(3, 2), function(x) any(x > 0))
  Tn_matrix <- apply(attr(catCH2, "Tm"), c(1, 2), function(x) any(x > 0))
  Tu_matrix <- apply(attr(catCH2, "Tu"), c(1, 2), function(x) any(x > 0))
  combined <- ID_matrix | Tn_matrix | Tu_matrix # combine all detections
  detected_per_trap <- apply(combined, 1, any)# collapse across occasions
  # Events
  EventsID <- sum(caps)
  EventsTm <- sum(attr(catCH2, "Tm")) #marked no ID
  EventsTu <- sum(attr(catCH2, "Tu")) #unmarked
  EventsTn <- sum(attr(catCH2, "Tn")) #uncertain
  TotalEvents <- EventsID + EventsTm + EventsTu + EventsTn
  # Effort
  camCTN <- rowSums(attr(attr(catCH2, "traps"), "usage"))
  # Results
  res <- data.frame(
    Array = sub("\\.Rdata$", "", sub(".*/[0-9]+_", "", file_name)),
    `Area` = round((attr(mask, "area")*nrow(mask))/100, 0),
    `Camera_density` = round(nrow(attr(catCH2, "traps")) / ((attr(mask, "area")*nrow(mask))/100),1),
    `Cameras_deployed` = as.integer(nrow(attr(catCH2, "traps"))),
    `Cameras_detections` = as.integer(sum(detected_per_trap)), # number of unique traps with ≥1 detection
    `CTN` = as.integer(sum(attr(attr(catCH2, "traps"), "usage"))),
    `Min_CTN` = as.integer(min(camCTN)),
    `Second_Min_CTN` = as.integer(min(camCTN[camCTN != min(camCTN)])),
    `Mean_CTN` = as.integer(mean(camCTN)),
    `Median_CTN` =  as.integer(median(camCTN)),
    `Max_CTN` = as.integer(max(camCTN)),
    `No_Identified` = as.integer(length(caps)),
    `No_detections` = as.integer(TotalEvents),
    `No_unid_marked` = as.integer(EventsTm),
    `No_unmarked` = as.integer(EventsTu + EventsTn),
    `Recaptures` = as.integer(sum(pmax(0, caps - 1))),
    `Spatial_recaptures` = sum(unlist(secr::moves(catCH2)) > 0),
    `Mean_distance` = round(mean(unlist(secr::moves(catCH2))),0),
    `Max_distance` = round(max(unlist(secr::moves(catCH2))),0)
  )
  
  sumres[[file_name]] <- res
}

summary_captures = t(bind_rows(sumres) %>%
  mutate(model_id = case_when(  
                               Array == "UNI_Grid500m" ~ "A1", 
                               Array == "INF_Grid500m" ~ "A2",
                               Array == "UNI_Grid250m" ~ "B1", 
                               Array == "INF_Grid250m" ~ "B2",
                               Array == "UNI_Grid1000m" ~ "C1", 
                               Array == "INF_Grid1000m" ~ "C2",
                               Array == "UNI_Road" ~ "D1", 
                               Array == "INF_Road" ~ "D2",
                               Array == "INF_Grid500m_sparseB" ~ "SA1", 
                               Array == "INF_Grid500m_sparseA" ~ "SA2", 
                               Array == "INF_Grid1000m_areaA" ~ "SA3", 
                               Array == "INF_Grid1000m_areaB" ~ "SA4", 
                               Array == "INF_Grid250m_sparseB" ~ "SB1", 
                               Array == "INF_Grid250m_sparseA" ~ "SB2")) %>%
  mutate(model_id = factor(model_id, levels=c("A1",  "B1", "C1", "D1", 
                                              "A2", "B2", "C2", "D2", 
                                              "SA1", "SA2", "SA3", "SA4", "SB1", "SB2"))) %>%
  arrange(model_id))

