prepare_camera_data <- function(
    array_type,
    camera_locations_file,
    cat_data_files,
    cat_data_type,
    operation_tracking_files,
    study_start_date = "2023-06-30",
    study_end_date = "2023-09-10",
    output_path = NULL,
    output_filename = NULL
) {
  
  require(tidyverse)
  require(sf)
  
  # ===== LOAD TRAP INFORMATION ====================================================
  camlist <- read.csv(camera_locations_file) %>%
    rename(CameraCode = if_else(
      "CameraCode" %in% names(.), 
      "CameraCode", 
      if_else("name" %in% names(.), "name", names(.)[1])
    )) %>%
    arrange(CameraCode) # Don't forget this! It will ruin your day(s) if you do
  
  # Transform camera coordinates to study CRS (3308)
  initX <- st_as_sf(
    camlist %>% dplyr::select(CameraCode, Longitude, Latitude),
    coords = c("Longitude", "Latitude"),
    crs = 4326
  )
  X <- st_coordinates(st_transform(initX, crs = 3308, crs_from = 4326))
  rownames(X) <- initX$CameraCode
  
  Xdf <- as.data.frame(X)
  Xdf$Camera <- initX$CameraCode
  
  # ===== LOAD AND PROCESS CAT DETECTION DATA ===================================
  if (cat_data_type == "Informed") {
    catimages <- do.call(rbind, lapply(cat_data_files, read.csv)) %>%
      filter(!grepl("BADCAM", Keywords_raw)) %>%
      mutate(
        DateTimeOriginal = as.POSIXct(DateTimeOriginal,format = "%Y-%m-%d %H:%M:%S",tz = "Australia/Brisbane")
      ) %>%
      arrange(Camera, DateTimeOriginal) %>%
      group_by(Camera, FinalID) %>%
      mutate(
        time_diff = as.numeric(difftime(lead(DateTimeOriginal, 1), DateTimeOriginal, "secs")),
        time_diff = replace_na(time_diff, 1000)
      ) %>%
      filter(time_diff > 60) %>%
      dplyr::select(-time_diff) %>%
      ungroup() %>%
      mutate(
        hour = as.numeric(round(difftime(DateTimeOriginal,as.POSIXct(paste(study_start_date, "00:00:01"),format = "%Y-%m-%d %H:%M:%S",tz = "Australia/Brisbane"), unit = "hours"), 0)),
        Bicolor = if_else(Bicolor == "Y", 2, 1),
        sex = case_when(
          Sex == "MALE" ~ 2,
          Sex == "FEMALE" ~ 1,
          Coat == "GING_D" ~ 2,
          Coat == "TORT" ~ 1,
          .default = 0
        ),
        Primary_Pattern = case_when(
          Coat == "BLACK" ~ 1,
          Coat %in% c("WHITE_IR", "GING_D") ~ 2,
          Coat == "CLAS_TABBY" ~ 3,
          Coat == "MACK_TABBY" ~ 4,
          Coat == "TORT" ~ 5,
          .default = 0
        ),
        Longer_hair = if_else(Long_Hair == "Y", 2, 1),
        samp.type = case_when(
          FinalID != "NO_FINAL_ID" ~ 1,
          FinalID == "NO_FINAL_ID" & Primary_Pattern > 1 ~ 2,
          FinalID == "NO_FINAL_ID" & Primary_Pattern == 1 ~ 3,
          FinalID == "NO_FINAL_ID" & Primary_Pattern == 0 ~ 4
        ),
        adjustDateTime = DateTimeOriginal + 12 * 60 * 60,
        Dayinit = as.Date(adjustDateTime),
        Day = as.numeric(Dayinit - as.Date(study_start_date)) + 1,
        FinalID = factor(FinalID)
      ) %>%
      arrange(samp.type, FinalID) %>%
      filter(Day <= as.numeric(as.Date(study_end_date) - as.Date(study_start_date))) %>%
      dplyr::select(
        Camera, DateTimeOriginal, hour, Occasion = Day, ID.name = FinalID,
        samp.type, Primary_Pattern, Bicolor, sex, Longer_hair
      ) %>%
      filter(Camera %in% camlist$CameraCode)
  } else {
    catimages <- do.call(rbind, lapply(cat_data_files, read.csv)) %>%
      filter(!grepl("BADCAM", Keywords_raw)) %>%
      mutate(
        DateTimeOriginal = as.POSIXct(DateTimeOriginal,format = "%Y-%m-%d %H:%M:%S",tz = "Australia/Brisbane")
      ) %>%
      arrange(Camera, DateTimeOriginal) %>%
      group_by(Camera, VerifiedID) %>%
      mutate(
        time_diff = as.numeric(difftime(lead(DateTimeOriginal, 1), DateTimeOriginal, "secs")),
        time_diff = replace_na(time_diff, 1000)
      ) %>%
      filter(time_diff > 60) %>%
      dplyr::select(-time_diff) %>%
      ungroup() %>%
      mutate(
        hour = as.numeric(round(difftime(DateTimeOriginal,as.POSIXct(paste(study_start_date, "00:00:01"),format = "%Y-%m-%d %H:%M:%S",tz = "Australia/Brisbane"), unit = "hours"), 0)),
        Bicolor = if_else(Bicolor == "Y", 2, 1),
        sex = case_when(
          Sex == "MALE" ~ 2,
          Sex == "FEMALE" ~ 1,
          Coat == "GING_D" ~ 2,
          Coat == "TORT" ~ 1,
          .default = 0
        ),
        Primary_Pattern = case_when(
          Coat == "BLACK" ~ 1,
          Coat %in% c("WHITE_IR", "GING_D") ~ 2,
          Coat == "CLAS_TABBY" ~ 3,
          Coat == "MACK_TABBY" ~ 4,
          Coat == "TORT" ~ 5,
          .default = 0
        ),
        Longer_hair = if_else(Long_Hair == "Y", 2, 1),
        samp.type = case_when(
          VerifiedID != "NO_ID" ~ 1,
          VerifiedID == "NO_ID" & Primary_Pattern > 1 ~ 2,
          VerifiedID == "NO_ID" & Primary_Pattern == 1 ~ 3,
          VerifiedID == "NO_ID" & Primary_Pattern == 0 ~ 4
        ),
        adjustDateTime = DateTimeOriginal + 12 * 60 * 60,
        Dayinit = as.Date(adjustDateTime),
        Day = as.numeric(Dayinit - as.Date(study_start_date)) + 1,
        VerifiedID = factor(VerifiedID)
      ) %>%
      arrange(samp.type, VerifiedID) %>%
      filter(Day <= as.numeric(as.Date(study_end_date) - as.Date(study_start_date))) %>%
      dplyr::select(
        Camera, DateTimeOriginal, hour, Occasion = Day, ID.name = VerifiedID,
        samp.type, Primary_Pattern, Bicolor, sex, Longer_hair
      ) %>%
      filter(Camera %in% camlist$CameraCode)
  }
  
  # Check sex assigned within individual uninformed datasets 
  if(cat_data_type == "Uninformed") {
    if (array_type == "Grid250m"){
      catimages <- catimages %>%
        mutate(sex = case_when(
          ID.name %in% c("PL_2023_CAT02","PL_2023_CAT03","PL_2023_CAT05") ~ 2,
          ID.name == "PL_2023_CAT04" ~ 1,
          .default = sex))
    } else if (array_type == "Grid500m") {
      catimages <- catimages %>%
        mutate(sex = case_when(
          ID.name %in% c("PL_2023_CAT05","PL_2023_CAT07","PL_2023_CAT09","PL_2023_CAT20") ~ 2,
          ID.name == "PL_2023_CAT02" ~ 1,
          .default = sex))  
    } else if (array_type == "Grid1000m") {
      catimages <- catimages %>%
        mutate(sex = case_when(
          ID.name %in% c("PL_2023_CAT06","PL_2023_CAT08","PL_2023_CAT019") ~ 2,
          .default = sex))  
    } else if (array_type == "Road") {
      catimages <- catimages %>%
        mutate(sex = case_when(
          ID.name %in% c("PL_2023_CAT07","PL_2023_CAT18") ~ 2,
          ID.name == c("PL_2023_CAT19") ~ 1,
          .default = sex))  
    }
  }
  
  # ===== LOAD CAMERA OPERATION DATA ============================================
  operation_data <- do.call(rbind, lapply(operation_tracking_files, read.csv)) %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d", tz = "Australia/Brisbane")) %>%
    dplyr::select(Camera, Date, Status)
  
  Kinit <- operation_data %>%
    filter(
      Date >= as.Date(study_start_date),
      Date < as.Date(study_end_date),
      Camera %in% camlist$CameraCode,
      Status == "working"
    ) %>%
    mutate(bin = 1) %>%
    distinct(Camera, Date, .keep_all = TRUE) %>%
    arrange(Date) %>%
    pivot_wider(id_cols = Camera, names_from = Date, values_from = bin, values_fill = 0) %>%
    arrange(Camera) %>%
    column_to_rownames("Camera") %>%
    as.matrix()
  write.table(Kinit, paste0(output_path, "/", output_filename, "_operation.txt"), col.names = F, row.names = F, sep = " ", quote = F)
  
  # ===== CREATE DETECTION TABLES ============================================
  # MARKED
  marked <- catimages %>%
    filter(samp.type == 1) %>%
    mutate(Session = 1) %>%
    select(Session, ID.name, Occasion, Camera) %>%
    distinct() %>%
    arrange(ID.name)
  write.table(marked, paste0(output_path, "/", output_filename, "_marked_cats.txt"), col.names = F, row.names = F, sep = " ", quote = F)
  
  # COMBOS
  combos <- expand.grid(Camera = Xdf$Camera, Occasion = 1:72)
  
  # UNMARKED
  unmarked <- catimages %>%
    filter(samp.type %in% c(3,4)) %>%
    select(Camera, Occasion) %>%
    mutate(bin = 1) %>%
    distinct() %>%
    full_join(combos, by = c("Camera", "Occasion")) %>%
    arrange(Occasion, Camera) %>%
    mutate(Session = 1, bin = replace_na(bin, 0)) %>%
    pivot_wider(id_cols = c(Session, Camera), names_from = Occasion, values_from = bin) %>%
    select(-Camera)
  write.table(unmarked, paste0(output_path, "/", output_filename, "_unmarked_cats.txt"), col.names = F, row.names = F, sep = " ", quote = F)
  
  # MARKED UNIDENTIFIED
  noid <- catimages %>%
    filter(samp.type %in% c(2)) %>%
    select(Camera, Occasion) %>%
    mutate(bin = 1) %>%
    distinct() %>%
    full_join(combos, by = c("Camera", "Occasion")) %>%
    arrange(Occasion, Camera) %>%
    mutate(Session = 1, bin = replace_na(bin, 0)) %>%
    pivot_wider(id_cols = c(Session, Camera), names_from = Occasion, values_from = bin) %>%
    select(-Camera)
  write.table(noid, paste0(output_path, "/", output_filename, "_uncertain.txt"), col.names = F, row.names = F, sep = " ", quote = F)
  
  # DETECTORS
  detectors <- Xdf %>%
    select(Detector = Camera, X, Y) %>%
    arrange(Detector)
  write.table(detectors, paste0(output_path, "/", output_filename,  "_detectors.txt"), col.names = F, row.names = F, sep = " ", quote = F)
}