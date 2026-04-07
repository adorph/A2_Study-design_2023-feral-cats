prepare_nimble_data <- function(data_file, output_dir = NULL, M1_factor = 1.3, M2_factor = 1, verbose = TRUE) {
  
  require(tidyverse)
  require(abind)
  require(sf)
  require(terra)
  
  tryCatch({
    
    # ===== LOAD DATA ===========================================================
    if (!file.exists(data_file)) stop("Data file not found:", data_file)
    load(data_file)
    
    if (verbose) cat("✓ Loaded data from", basename(data_file), "\n")
    
    # ===== EXTRACT PARAMETERS ==================================================
    category_key <- data$ID.covariates$category_key
    gamma <- data$ID.covariates$gamma
    n.cat <- data$ID.covariates$n.cat
    n.levels <- unlist(lapply(data$ID.covariates$IDcovs, length))
    
    y.ID <- data$detections$y.marked
    y.noID <- abind(data$detections$y.marked.noID, data$detections$y.unmarked, 
                    data$detections$y.unk, along = 1)
    n.ID <- data$detections$n.ID
    
    X <- data$trap.covariates$X
    J <- data$trap.covariates$J
    K <- data$occasion.covariates$K
    K1D <- data$occasion.covariates$K1D
    buff <- data$trap.covariates$buff
    xlim <- data$trap.covariates$xlim
    ylim <- data$trap.covariates$ylim
    
    IDcovs <- data$ID.covariates$IDcovs
    G.ID <- data$ID.covariates$G.marked
    G.noID <- abind(data$ID.covariates$G.marked.noID, data$ID.covariates$G.unmarked, 
                    data$ID.covariates$G.unk, along = 1)
    
    # ===== VALIDATE DATA STRUCTURES ============================================
    if (!is.matrix(G.ID)) G.ID <- matrix(G.ID)
    if (!is.matrix(G.noID)) G.noID <- matrix(G.noID)
    if (!is.list(IDcovs)) stop("IDcovs must be a list")
    if (ncol(G.ID) != n.cat) stop("G.ID needs n.cat number of columns")
    if (ncol(G.noID) != n.cat) stop("G.noID needs n.cat number of columns")
    if (length(dim(y.ID)) != 3) stop("dim(y.ID) must be 3")
    if (length(dim(y.noID)) != 3) stop("dim(y.noID) must be 3")
    if (!is.list(gamma)) stop("gamma must be a list")
    
    if (verbose) cat("✓ Data structures validated\n")
    
    # ===== POPULATION AUGMENTATION =============================================
    # Initial augmentation
    # M1 <- ceiling(n.ID * 16 * M1_factor)
    # M2 <- ceiling(M1 / 2 * M2_factor)
    # M <- M1 + M2
    # More conservative augmentation
    M1 <- ceiling(n.ID / 2 * 16 * M1_factor) # augment female population
    M2 <- ceiling(n.ID / 2 * 16 * M2_factor) # augment male population
    M <- M1 + M2
    
    if (verbose) cat("✓ Augmentation: M1 =", M1, ", M2 =", M2, ", M =", M, "\n")
    
    # ===== INITIALIZE IDS ======================================================
    G.true <- matrix(0, nrow = M, ncol = n.cat)
    G.true[1:n.ID, ] <- G.ID
    
    n.samples <- nrow(y.noID)
    ID <- rep(NA, n.samples)
    nextID <- n.ID + 1
    
    y.noID2D <- apply(y.noID, c(1, 2), sum)
    this.j <- apply(y.noID2D, 1, function(x) which(x > 0))
    y.ID2D <- apply(y.ID, c(1, 2), sum)
    y.ID2D <- rbind(y.ID2D, matrix(0, nrow = M - n.ID, ncol = J))
    y.true2D <- y.ID2D
    
    # Match samples with consistent covariates
    for (l in 1:n.samples) {
      obs.id <- which(G.noID[l, ] != 0)
      matches <- which(apply(G.true, 1, function(row) {
        obs.id.x2 <- which(row != 0)
        sameobsidx <- intersect(obs.id, obs.id.x2)
        if (length(sameobsidx) == 0) return(FALSE)
        all(row[sameobsidx] == G.noID[l, sameobsidx])
      }))
      
      sametrap <- if (length(matches) > 0) y.true2D[matches, this.j[[l]]] > 0 else logical(0)
      
      if (any(sametrap)) {
        ID[l] <- matches[which(sametrap)[1]]
        y.true2D[ID[l], this.j[[l]]] <- y.true2D[ID[l], this.j[[l]]] + 1
        notobsidx <- which(G.true[ID[l], ] == 0)
        G.true[ID[l], notobsidx] <- G.noID[l, notobsidx]
      } else {
        if (nextID > M) stop("Need to raise M to initialize data")
        ID[l] <- nextID
        G.true[nextID, obs.id] <- G.noID[l, obs.id]
        y.true2D[nextID, this.j[[l]]] <- y.true2D[nextID, this.j[[l]]] + 1
        nextID <- nextID + 1
      }
    }
    
    # Verify initialization
    for (l in 1:n.samples) {
      obs.id <- which(G.noID[l, ] != 0)
      if (!all(G.true[ID[l], obs.id] == G.noID[l, obs.id])) 
        stop("Error in initialization at sample", l)
    }
    
    # Fill remaining covariates
    inits_gamma <- gamma
    for (i in 1:M) {
      for (c in 1:n.cat) {
        if (G.true[i, c] == 0) {
          G.true[i, c] <- sample(IDcovs[[c]], 1, prob = inits_gamma[[c]])
        }
      }
    }
    
    # Sex-dependent coat color constraints
    for (i in 1:M) {
      if (G.true[i, 2] == 2) G.true[i, 1] <- 2
      if (G.true[i, 2] == 5) G.true[i, 1] <- 1
    }
    
    z <- 1 * (rowSums(y.true2D) > 0)
    
    if (verbose) cat("✓ IDs initialized for", sum(z), "individuals\n")
    
    # ===== CREATE HABITAT MASK ==================================================
    initX <- st_as_sf(as.data.frame(X), coords = c('X', 'Y'), crs = 3308)
    buffzone <- st_make_valid(st_union(st_buffer(initX, dist = buff)))
    area <- as.numeric(st_area(buffzone) / 1000000)
    
    vbuff <- vect(buffzone)
    mask_r <- rast(ext(vbuff), res = 150)
    InSS <- rasterize(vect(buffzone), mask_r, background = 0)
    
    res <- res(InSS)
    x.vals <- xFromCol(InSS)
    y.vals <- rev(yFromRow(InSS))
    dSS <- as.matrix(cbind(expand.grid(x.vals, y.vals)))
    cells <- matrix(1:nrow(dSS), nrow = length(x.vals), ncol = length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    
    InSS2 <- flip(InSS, direction = "vertical")
    InSS.tmp <- as.vector(terra::values(InSS2))
    
    # Validate mask
    if (nrow(dSS) != n.cells) stop("'dSS' must have 'n.cells' rows")
    if (ncol(dSS) != 2) stop("'dSS' must have 2 columns")
    if (dim(cells)[1] != n.cells.x) stop("'cells' should be of length 'n.cells.x'")
    if (dim(cells)[2] != n.cells.y) stop("'cells' should be of length 'n.cells.y'")
    
    xlim <- range(x.vals) + c(-res[1], res[1])
    ylim <- range(y.vals) + c(-res[2], res[2])
    
    
    if (!all(range(x.vals) + c(-res[1], res[1]) == xlim)) 
      stop("'xlim' does not match 'xvals'")
    if (!all(range(y.vals) + c(-res[2], res[2]) == ylim)) 
      stop("'ylim' does not match 'yvals'")    
    
    if (verbose) cat("✓ Habitat mask created:", n.cells, "cells\n")

    
    # ===== COMPREHENSIVE MASK VALIDATION ========================================
    if (verbose) cat("Running comprehensive mask validation...\n")
    
    for (i in 1:n.cells) {
      s.cell.x <- i %% n.cells.x
      s.cell.y <- floor(i / n.cells.x) + 1
      if (s.cell.x == 0) {
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y - 1
      }
      
      match <- which(cells == i, arr.ind = TRUE)
      if (!all(c(s.cell.x, s.cell.y) == match)) 
        stop("Mask validation error 1 at cell", i, ": cell indexing mismatch")
      
      xlim.cell <- c(s.cell.x - 1, s.cell.x) * res[1]
      ylim.cell <- c(s.cell.y - 1, s.cell.y) * res[2]
      
      if (max(x.vals[s.cell.x] + c(-res[1]/2, res[1]/2) - xlim.cell) > max(xlim))
        stop("Mask validation error 2 at cell", i, ": x-limits mismatch")
      if (max(y.vals[s.cell.y] + c(-res[2]/2, res[2]/2) - ylim.cell) > max(ylim))
        stop("Mask validation error 3 at cell", i, ": y-limits mismatch")
    }
    
    if (verbose) cat("✓ Mask validation passed (", n.cells, "cells checked)\n")
    
    # ===== INITIALIZE ACTIVITY CENTERS =========================================
    s <- st_coordinates(st_sample(buffzone, M))
    idx <- which(rowSums(y.true2D) > 0)
    
    for (i in idx) {
      trps <- matrix(X[y.true2D[i, ] > 0, 1:2], ncol = 2, byrow = FALSE)
      s[i, ] <- if (nrow(trps) > 1) c(mean(trps[, 1]), mean(trps[, 2])) else trps
    }
    
    # Adjust points outside habitat to nearest habitat cell
    e2dist <- function(x, y) {
      i <- sort(rep(1:nrow(y), nrow(x)))
      dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
      matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = FALSE)
    }
    
    alldists <- e2dist(s, dSS)
    alldists[, InSS.tmp == 0] <- Inf
    
    # BUG FIX 20260304: Add explicit bounds checking
    out_of_bounds_individuals <- c()
    for (i in 1:M) {
      # Calculate cell indices
      cell_idx_x <- trunc(n.cells.x - (xlim[2] - s[i, 1]) / res[1]) + 1
      cell_idx_y <- trunc(n.cells.y - (ylim[2] - s[i, 2]) / res[2]) + 1
      
      # Check bounds explicitly
      if (cell_idx_x < 1 || cell_idx_x > n.cells.x || 
          cell_idx_y < 1 || cell_idx_y > n.cells.y) {
        out_of_bounds_individuals <- c(out_of_bounds_individuals, i)
        
        if (verbose) {
          cat("WARNING: Individual", i, "- coords (", s[i, 1], ",", s[i, 2], 
              ") map to invalid cell (", cell_idx_x, ",", cell_idx_y, 
              ") - snapping to nearest habitat\n")
        }
        
        # Force snap to nearest habitat cell
        new.cell <- which.min(alldists[i, ])
        s[i, ] <- dSS[new.cell, ]
        
      } else {
        # Safe to index
        this.cell <- cells[cell_idx_x, cell_idx_y]
        
        if (InSS.tmp[this.cell] == 0) {
          new.cell <- which.min(alldists[i, ])
          s[i, ] <- dSS[new.cell, ]
          
          if (verbose && i <= 5) {  # Log first 5
            cat("Individual", i, "adjusted to nearest habitat\n")
          }
        }
      }
    }
    if (length(out_of_bounds_individuals) > 0) {
      cat("\n⚠ WARNING:", length(out_of_bounds_individuals), 
          "individuals were outside bounds and snapped to habitat\n")
      cat("Individual IDs:", out_of_bounds_individuals, "\n")
    }
    if (verbose) cat("✓ Activity centers initialized\n")
    
    # ===== VALIDATE ALL ACTIVITY CENTERS ARE IN HABITAT ===========================
    cat("Validating all activity centers...\n")
    
    invalid_individuals <- c()
    
    for (i in 1:M) {
      cell_idx_x <- trunc(n.cells.x - (xlim[2] - s[i, 1]) / res[1]) + 1
      cell_idx_y <- trunc(n.cells.y - (ylim[2] - s[i, 2]) / res[2]) + 1
      
      # Bounds check
      if (cell_idx_x < 1 || cell_idx_x > n.cells.x || 
          cell_idx_y < 1 || cell_idx_y > n.cells.y) {
        invalid_individuals <- c(invalid_individuals, i)
        cat("✗ Individual", i, "STILL out of bounds:", s[i, ], "\n")
        next
      }
      
      this.cell <- cells[cell_idx_x, cell_idx_y]
      
      if (InSS.tmp[this.cell] == 0) {
        invalid_individuals <- c(invalid_individuals, i)
        cat("✗ Individual", i, "in non-habitat cell\n")
      }
    }
    
    if (length(invalid_individuals) > 0) {
      stop("ERROR: ", length(invalid_individuals), 
           " individuals still invalid after adjustment:\n",
           paste(invalid_individuals, collapse = ", "))
    }
    
    cat("✓ All", M, "activity centers verified within habitat\n")
    
    # ===== CALCULATE DETECTION RATES ===========================================
    lam0 <- c(exp(-3.5), exp(-3.5))
    sigma <- c(exp(5.7), exp(5.7) * 1.67)
    
    D <- e2dist(s, X)
    lamd <- lam0[1] * exp(-D^2 / (2 * sigma[1]^2))
    
    for (i in seq_len(nrow(lamd))) {
      lamd[i, ] <- lam0[which(G.true[i, 1] == 1:length(lam0))] * 
        exp(-D[i, ]^2 / (2 * sigma[which(G.true[i, 1] == 1:length(sigma))]^2))
    }
    
    # Validate likelihood
    ll.y <- dpois(y.true2D, K1D * lamd * z, log = TRUE)
    if (!is.finite(sum(ll.y))) 
      stop("Starting likelihood not finite. Check K1D or initialization.")
    
    if (verbose) cat("✓ Detection rates calculated\n")
    
    # ===== PREPARE CATEGORICAL PARAMETERS ======================================
    gammaMat <- matrix(0, nrow = n.cat, ncol = max(n.levels))
    for (l in 1:n.cat) gammaMat[l, 1:n.levels[l]] <- gamma[[l]]
    
    G.latent <- matrix(1, nrow = M, ncol = n.cat)
    # BUG FIX: Previously set all known cats to zero for all covariates - but sex is not always known for identified individuals
    # G.latent.ID <- 0*(G.ID>0)
    # Adjusted to set only known values to 0 for known cats - otherwise allow the model to sample the values
    G.latent[1:n.ID, ] <- ifelse(G.ID > 0, 0, 1)
    
    for (l in 1:n.samples) {
      for (m in 1:n.cat) {
        if (G.noID[l, m] != 0) G.latent[ID[l], m] <- 0
      }
    }
    
    if (verbose) cat("✓ Categorical parameters prepared\n")
    
    # ===== COMPILE NIMBLE BUILD LIST ===========================================
    nimbuild <- list(
      M = M, s = s, z = z, G.true = G.true, ID = ID, G.noID = G.noID,
      y.ID = y.ID2D, y.true = y.true2D, K1D = K1D, n.samples = n.samples,
      gammaMat = gammaMat, this.j = this.j, G.latent = G.latent,
      G.latent.ID = 0 * (G.ID > 0), area = area, cells = cells, res = res,
      InSS = InSS.tmp, n.cells = n.cells, n.cells.x = n.cells.x, n.cells.y = n.cells.y
    )
    
    # ===== SAVE OUTPUT =========================================================
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      # Update and save data
      data$trap.covariates$xlim <- xlim
      data$trap.covariates$ylim <- ylim
      data_out <- file.path(output_dir, sub(".Rdata$", ".Rdata", basename(data_file)))
      save(data, file = data_out)
      
      # Save nimbuild
      nimbuild_out <- file.path(output_dir, sub(".Rdata$", "_nimbuild.Rdata", basename(data_file)))
      save(nimbuild, file = nimbuild_out)
      
      if (verbose) {
        cat("✓ Saved updated data to:", basename(data_out), "\n")
        cat("✓ Saved nimbuild to:", basename(nimbuild_out), "\n")
      }
    }
    
    if (verbose) cat("\n✓ NIMBLE data preparation complete\n")
    return(nimbuild)
    
  }, error = function(e) {
    cat("\n✗ ERROR in prepare_nimble_data():\n", conditionMessage(e), "\n")
    return(NULL)
  })
}