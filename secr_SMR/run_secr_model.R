run_secr_model <- function(
    data_file, 
    Tu_file,
    Tm_file, 
    trap_file, 
    operation_file,
    combo_name,
    output_dir = "./outputs/",
    verbose = TRUE
) {
  
  require(secr)
  require(tidyverse)
  require(sf)
  
  tryCatch({
    
    # Load data
    if (!file.exists(data_file)) stop("Marked file not found:", data_file)
    if (!file.exists(Tu_file)) stop("Tu_file file not found:", Tu_file)
    if (!file.exists(Tm_file)) stop("Tm_file file not found:", Tm_file)
    if (!file.exists(trap_file)) stop("Detector file not found:", trap_file)
    if (!file.exists(operation_file)) stop("Operation file not found:", operation_file)
    
    captures <- read.table(data_file, col.names = c("Session", "ID", "occassion", "trap"))
    trapdata <- read.table(trap_file, col.names = c("trap", "x", "y")) %>% arrange(trap)
    operdata <- read.table(operation_file)
    
    # Create traps object
    traps <- secr::read.traps(data = trapdata,
                        detector = "proximity",
                        trapID = "trap",
                        markocc = rep(0, 72))  # number of marking occasions specified
    camlist <- unique(trapdata$trap)
    
    # Bind usage to traps data
    usage(traps) <- operdata
    
    # Create capture history (mark-recapture)
    catCH <- make.capthist(captures = captures, traps = traps,
                           fmt = "trapID", noccasions = c(72, 72))
    
    # Add unmarked sightings
    catCH2 <- addSightings(catCH, unmarked=Tu_file, nonID=Tm_file)
    
    # Make mask
    if(grepl("UNI_Road", combo_name)) {
      mask <- make.mask(traps(catCH2), buffer = 8000, spacing = 150, type = 'trapbuffer')
    } else {
      mask <- make.mask(traps(catCH2), buffer = 3000, spacing = 150, type = 'trapbuffer')
    }
    
    # Fit base model
    fit_M1 <- secr.fit(catCH2, mask = mask,
                       model = list(D ~ 1, g0 ~ 1, sigma ~ 1), 
                       detectfn = 'HN',
                       trace = FALSE,
                       details = list(knownmarks = FALSE))
    
    # Estimate overdispersion
    fit_M2 <- secr.fit(catCH2, mask = mask, 
                       model = list(D ~ 1, g0 ~ 1, sigma ~ 1), 
                       detectfn = 'HN',
                       trace = F, 
                       details = list(knownmarks = FALSE, nsim = 10000),
                       start = fit_M1)
    
    # Fit final model adjusting for overdispersion
    fit_M3 <- secr.fit(catCH2, mask = mask, 
                       model = list(D ~ 1, g0 ~ 1, sigma ~ 1), 
                       detectfn = 'HN', 
                       trace = F,
                       details = list(knownmarks = FALSE, chat = fit_M2$details$chat),
                       start = fit_M2) 
    
    # Save all
    save(list = c("fit_M1", "fit_M2", "fit_M3", "catCH2", "mask"),
         file = paste0(output_dir, "/", format(Sys.Date(), "%Y%m%d"), "_", combo_name, ".Rdata"))
    
    # Report
    cat("Finished:", combo_name, "\n\n")
    
  }, error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
  
}