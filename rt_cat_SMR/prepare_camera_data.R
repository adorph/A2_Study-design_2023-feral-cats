prepare_camera_data <- function(
    array_type,
    camera_locations_file,
    cat_data_files,
    cat_data_type,
    operation_tracking_files,
    include_longer_hair = TRUE,
    study_start_date = "2023-06-30",
    study_end_date = "2023-09-10",
    output_path = NULL,
    output_filename = NULL
) {
  
  require(tidyverse)
  require(sf)
  require(abind)
  
  # array_type <- match.arg(array_type)
  
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
  rm(initX)
  
  J <- nrow(X)
  if(grepl("Road", array_type) & cat_data_type == "Uninformed") {
    buff <- 8000
  }else {
    buff <- 3000
  }

  xlim <- c(min(X[, 1]) - buff, max(X[, 1]) + buff)
  ylim <- c(min(X[, 2]) - buff, max(X[, 2]) + buff)
  
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

  
  # Join camera locations
  catimages <- catimages %>%
    left_join(Xdf, by = "Camera")
  
  # ===== EXTRACT ID INFORMATION ================================================
  if (cat_data_type == "Informed") {
    ID.marked <- unique(
      catimages %>%
        filter(!grepl("NO_FINAL_ID", ID.name)) %>%
        dplyr::select(ID.name) %>%
        droplevels()
    )$ID.name
  } else {
    ID.marked <- unique(
      catimages %>%
        filter(!grepl("NO_ID", ID.name)) %>%
        dplyr::select(ID.name) %>%
        droplevels()
    )$ID.name
  }
  
  n.ID <- length(ID.marked)
  K <- as.numeric(as.Date(study_end_date) - as.Date(study_start_date))
  
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
    group_by(Camera) %>%
    summarise(Freq = sum(bin)) %>%
    mutate(Freq = if_else(Freq > K, K, Freq))
  
  K1D <- as.vector(Kinit$Freq)
  names(K1D) <- Kinit$Camera
  
  # ===== BUILD DETECTION MATRICES ==============================================
  # Helper function to build detection array for a given sampling type
  build_detection_array <- function(samp_type_val, n_individuals, J, K, cameras, catimages_df) {
    initimages <- catimages_df %>%
      filter(samp.type == samp_type_val) %>%
      dplyr::select(Camera, Occasion, ID.name)
    
    n_ind <- nrow(initimages)
    if (n_ind == 0) {
      return(array(numeric(), c(0, J, K)))
    }
    
    expand_grid <- as.data.frame(
      expand.grid(
        Camera = cameras,
        Occasion = seq(1, K),
        ID.name = seq(1:n_ind)
      )
    ) %>%
      mutate(
        Camera = as.character(Camera),
        Occasion = as.integer(Occasion),
        ID.name = as.character(ID.name)
      )
    
    init <- initimages %>%
      mutate(bin = 1) %>%
      mutate(ID.name = as.character(row_number())) %>%
      full_join(expand_grid, by = c("Camera", "Occasion", "ID.name")) %>%
      arrange(Occasion) %>%
      pivot_wider(
        id_cols = c("ID.name", "Camera"),
        names_from = "Occasion",
        values_from = bin
      ) %>%
      arrange(ID.name, Camera) %>%
      replace(is.na(.), 0)
    
    y_array <- array(numeric(), c(n_ind, J, K))
    for (i in 1:n_ind) {
      init2 <- as.matrix(
        init %>%
          filter(ID.name == i) %>%
          dplyr::select(-ID.name) %>%
          arrange(Camera) %>%
          column_to_rownames(var = "Camera")
      )
      y_array[i, , ] <- init2
    }
    
    return(y_array)
  }
  
  # Event 1 - marked identified individuals
  expand1 <- as.data.frame(
    expand.grid(
      Camera = camlist$CameraCode,
      Occasion = seq(1, K),
      ID.name = ID.marked
    )
  ) %>%
    mutate(
      Camera = as.character(Camera),
      Occasion = as.integer(Occasion),
      ID.name = as.character(ID.name)
    )
  
  y.IDinit <- catimages %>%
    filter(samp.type == 1) %>%
    dplyr::select(Camera, Occasion, ID.name) %>%
    distinct() %>%
    mutate(bin = 1) %>%
    full_join(expand1, by = c("Camera", "Occasion", "ID.name")) %>%
    arrange(Occasion) %>%
    pivot_wider(
      id_cols = c("ID.name", "Camera"),
      names_from = "Occasion",
      values_from = bin
    ) %>%
    arrange(ID.name, Camera) %>%
    replace(is.na(.), 0)
  
  y.ID <- array(numeric(), c(n.ID, J, K))
  for (i in 1:length(ID.marked)) {
    y.IDinits <- as.matrix(
      y.IDinit %>%
        filter(ID.name == ID.marked[i]) %>%
        dplyr::select(-ID.name) %>%
        arrange(Camera) %>%
        column_to_rownames(var = "Camera")
    )
    y.ID[i, , ] <- y.IDinits
  }
  rm(y.IDinit, y.IDinits, expand1)
  y.marked <- y.ID
  
  # Events 2-4 using helper function
  y.marked.noID <- build_detection_array(2, n.ID, J, K, camlist$CameraCode, catimages)
  y.marked.noID[y.marked.noID > 1] <- 1
  
  y.unmarked <- build_detection_array(3, n.ID, J, K, camlist$CameraCode, catimages)
  
  y.unk <- build_detection_array(4, n.ID, J, K, camlist$CameraCode, catimages)
  
  # ===== BUILD ID COVARIATES ===================================================
  # Define gamma distributions
  sex_gamma <- as.table(c(0.5, 0.5))
  names(sex_gamma) <- c(1, 2)
  
  pat_gamma <- as.table(c(0.4, 0.15, 0.07, 0.21, 0.17))
  names(pat_gamma) <- c(1, 2, 3, 4, 5)
  
  bicol_gamma <- as.table(c(0.95, 0.05))
  names(bicol_gamma) <- c(1, 2)
  
  hair_gamma <- as.table(c(0.99, 0.01))
  names(hair_gamma) <- c(1, 2)
  
  # Select covariates based on include_longer_hair
  if (include_longer_hair) {
    gamma <- list(sex_gamma, pat_gamma, bicol_gamma, hair_gamma)
    names(gamma) <- c("sex", "Primary_Pattern", "Bicolor", "Longer_hair")
    covariate_cols <- c("sex", "Primary_Pattern", "Bicolor", "Longer_hair")
    n.levels <- as.integer(c(2, 5, 2, 2))
  } else {
    gamma <- list(sex_gamma, pat_gamma, bicol_gamma)
    names(gamma) <- c("sex", "Primary_Pattern", "Bicolor")
    covariate_cols <- c("sex", "Primary_Pattern", "Bicolor")
    n.levels <- as.integer(c(2, 5, 2))
  }
  
  names(n.levels) <- names(gamma)
  
  IDcovs <- lapply(
    list(
      sex = c(1, 2),
      Primary_Pattern = c(1, 2, 3, 4, 5),
      Bicolor = c(1, 2),
      Longer_hair = c(1, 2)
    ),
    as.integer
  )
  
  if (!include_longer_hair) {
    IDcovs$Longer_hair <- NULL
  }
  
  category_key <- list(
    sex = data.frame(IDcovs = as.integer(c(1, 2)), value = c("Female", "Male")),
    Primary_Pattern = data.frame(
      IDcovs = as.integer(c(1, 2, 3, 4, 5)),
      value = c("BLACK", "GING/WHITE", "CLAS_TABBY", "MACK_TABBY", "TORT")
    ),
    Bicolor = data.frame(IDcovs = as.integer(c(1, 2)), value = c("FALSE", "TRUE")),
    Longer_hair = data.frame(IDcovs = as.integer(c(1, 2)), value = c("FALSE", "TRUE"))
  )
  
  if (!include_longer_hair) {
    category_key$Longer_hair <- NULL
  }
  
  n.cat <- as.integer(length(gamma))
  
  # Build covariate matrices
  G.marked <- as.matrix(
    catimages %>%
      filter(samp.type == 1) %>%
      dplyr::select(ID.name, all_of(covariate_cols)) %>%
      distinct() %>%
      column_to_rownames(var = "ID.name")
  )
  
  G.marked.noID <- as.matrix(
    catimages %>%
      filter(samp.type == 2) %>%
      dplyr::select(all_of(covariate_cols))
  )
  
  G.unmarked <- as.matrix(
    catimages %>%
      filter(samp.type == 3) %>%
      dplyr::select(all_of(covariate_cols))
  )
  
  G.unk <- as.matrix(
    catimages %>%
      filter(samp.type == 4) %>%
      dplyr::select(all_of(covariate_cols))
  )
  
  G.obs <- rbind(G.marked, G.marked.noID, G.unmarked, G.unk)
  
  # ===== CONSTRUCT SAMPLING TYPE VECTOR ========================================
  samp.type <- c(rep(1, n.ID), catimages[catimages$samp.type != 1, ]$samp.type)

  # ===== COMPILE FINAL DATA LIST ===============================================
  data <- list(
    raw = catimages,
    detections = list(
      y.marked = y.marked,
      y.marked.noID = y.marked.noID,
      y.unmarked = y.unmarked,
      y.unk = y.unk,
      ID.marked = ID.marked,
      n.ID = n.ID,
      samp.type = samp.type
    ),
    occasion.covariates = list(K = K, K1D = K1D),
    trap.covariates = list(
      X = X,
      xlim = xlim,
      ylim = ylim,
      buff = buff,
      J = J
    ),
    ID.covariates = list(
      gamma = gamma,
      IDcovs = IDcovs,
      n.cat = n.cat,
      n.levels = n.levels,
      category_key = category_key,
      G.marked = G.marked,
      G.marked.noID = G.marked.noID,
      G.unmarked = G.unmarked,
      G.unk = G.unk,
      G.obs = G.obs
    ),
    metadata = list(
      array_type = array_type,
      include_longer_hair = include_longer_hair,
      study_start = study_start_date,
      study_end = study_end_date,
      processed_date = Sys.Date()
    )
  )
  
  # ===== SAVE OUTPUT ===========================================================
  if (!is.null(output_path)) {
    if (is.null(output_filename)) {
      output_filename <- paste0("2023_PL_", array_type, "_data.Rdata")
    }
    save(data, file = file.path(output_path, output_filename))
    cat("Data saved to:", file.path(output_path, output_filename), "\n")
  }
  
  return(data)
}