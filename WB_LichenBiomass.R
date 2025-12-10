defineModule(sim, list(
  name = "WB_LichenBiomass",
  description = paste("Compute lichen biomass according to the Greuel & Degre-Timmons model (2021)"),
  keywords = c("lichen", "biomass", "western boreal"),
  authors =  c(
    person("Pierre", "Racine", email= "pierre.racine@sbf.ulaval.ca", role = "aut"),
    person("Andres", "Caseiro Guilhem", email= "andres.caseiro-guilhem.1@ulaval.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(WB_LichenBiomass = "0.0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  # citation = list("citation.bib"),
  # documentation = list("NEWS.md", "README.md", "WB_LichenBiomass.Rmd"),
  reqdPkgs = list("reproducible"),
  loadOrder = list(after = c("WB_HartJohnstoneForestClasses")),
  parameters = rbind(
    defineParameter("WB_LichenBiomassTimeStep", "numeric", 1, NA, NA,
                    "Simulation time at which the drainage map is regenerated.")
  ),
  inputObjects = rbind(
    expectsInput(objectName = "WB_HartJohnstoneForestClassesMap",
                 objectClass = "SpatRast",
                 desc = paste("WB_HartJohnstoneForestClassesMap from the ",
                              "WB_HartJohnstoneForestClasses module used as ",
                              "standtype"),
                 sourceURL = NA),
    expectsInput(objectName = "EcoProvincesMap", 
                 objectClass = "SpatVector", 
                 desc = "", 
                 sourceURL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip"),
    expectsInput(objectName = "ageMap", 
                 objectClass = "SpatVector", 
                 desc = "", 
                 sourceURL = NA)
  ),
  outputObjects = rbind(
    createsOutput(objectName = "WB_LichenBiomassMap", 
                  objectClass = "SpatRast", 
                  desc = "Lichen biomass map predicted from the model")
  )
))

doEvent.WB_LichenBiomass = function(sim, eventTime, eventType) {
  switch(
    eventType,
    
    init = {
      sim <- scheduleEvent(sim, time(sim), "WB_LichenBiomass", "reComputeLichenBiomassMap", 2)
    },
    
    reComputeLichenBiomassMap = {
      sim <- reComputeLichenBiomassMap(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$WB_LichenBiomassTimeStep, "WB_LichenBiomass", "reComputeLichenBiomassMap")
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}
# ---------------------------------------------------
# 1. Model coefficients from Greuel et al. (2021)
# ---------------------------------------------------
# Logistic regression for presence/absence (g2)
g2_coefs <- list(
  Intercept = 0,                       # baseline = "conimix" (level 3)
  standtypeDeciduous = -0.801061,      # level 1 
  standtypeDeciduous_conifer_mix = -0.380172, # level 2
  standtypeJack_pine = 0.953791,       # level 4
  standtypeLarch = -15.068047,         # level 5
  standtypePoorly_drained_spruce = 0.395372, # level 7
  standtypeWell_drained_spruce = -0.587152,  # level 6
  ecoprovinceCoppermine = -0.398703,
  ecoprovinceGreat_Bear_Plains = -0.602398,
  ecoprovinceHay_Slave_River = -0.724421,
  TSSRF = 0.014524
)

# Non-linear biomass model (m4)
m4_alpha <- list(
  Intercept = 912.7731,              # baseline = "conimix" (level 3)
  standtypeJack_pine = 519.2378,     # level 4
  standtypePoorly_drained_spruce = 1379.2158, # level 7
  standtypeWell_drained_spruce = 178.3520     # level 6
)
m4_xmid <- list(
  Intercept = 28.0413,
  ecoprovinceCoppermine = 17.2021,
  ecoprovinceGreat_Bear_Plains = 44.9796,
  ecoprovinceHay_Slave_River = 41.7653
)
m4_scal <- 3.0647

# Function predicting lichen biomass (kg/ha).
# TSSRF (Time Since Stand-Replacing Fire) is approximated by sim$standAgeMap
predict_lichen_biomass <- function(ecoprov, standtype, TSSRF) {
  # --- Logistic regression (presence probability) ---
  logit <- g2_coefs$Intercept +
    ifelse(standtype == 1, g2_coefs$standtypeDeciduous, 0) +
    ifelse(standtype == 2, g2_coefs$standtypeDeciduous_conifer_mix, 0) +
    ifelse(standtype == 4, g2_coefs$standtypeJack_pine, 0) +
    ifelse(standtype == 5, g2_coefs$standtypeLarch, 0) +
    ifelse(standtype == 6, g2_coefs$standtypeWell_drained_spruce, 0) +
    ifelse(standtype == 7, g2_coefs$standtypePoorly_drained_spruce, 0) +
    ifelse(ecoprov == 0, 0,
           ifelse(ecoprov == 1, g2_coefs$ecoprovinceCoppermine,
                  ifelse(ecoprov == 2, g2_coefs$ecoprovinceGreat_Bear_Plains,
                         g2_coefs$ecoprovinceHay_Slave_River))) +
    g2_coefs$TSSRF * TSSRF
  
  p_presence <- 1 / (1 + exp(-logit))
  
  # --- Nonlinear biomass (only for valid stands) ---
  valid <- standtype %in% c(3, 4, 6, 7)  # conimix, jackpine, wd_spruce, pd_spruce
  
  alpha <- rep(NA, length(standtype))
  xmid  <- rep(NA, length(standtype))
  
  alpha[valid] <- m4_alpha$Intercept +
    ifelse(standtype[valid] == 4, m4_alpha$standtypeJack_pine, 0) +
    ifelse(standtype[valid] == 6, m4_alpha$standtypeWell_drained_spruce, 0) +
    ifelse(standtype[valid] == 7, m4_alpha$standtypePoorly_drained_spruce, 0)
  
  xmid[valid] <- m4_xmid$Intercept +
    ifelse(ecoprov[valid] == 1, m4_xmid$ecoprovinceCoppermine,
           ifelse(ecoprov[valid] == 2, m4_xmid$ecoprovinceGreat_Bear_Plains,
                  m4_xmid$ecoprovinceHay_Slave_River))
  
  biomass <- rep(NA, length(standtype))
  biomass[valid] <- alpha[valid] / (1 + exp((xmid[valid] - TSSRF[valid]) / m4_scal))
  
  # --- Combine (expected biomass = probability × biomass) ---
  expected_biomass <- p_presence * biomass
  
  return(expected_biomass)
}

reComputeLichenBiomassMap <- function(sim) {
  browser()
  message("Recomputing sim$WB_LichenBiomassMap for ", 
          format(ncell(sim$WB_HartJohnstoneForestClassesMap), scientific = FALSE), " pixels..")

  # Get a current age map from cohortData 
  currentCohortAgeMap <- LandR::standAgeMapGenerator(
    sim$cohortData,
    sim$pixelGroupMap
  )
  names(currentCohortAgeMap) <- "current_age"
  varnames(currentCohortAgeMap) <- "current_age"
  
  sim$WB_LichenBiomassMap <- app(c(
    sim$EcoProvincesMap, 
    sim$WB_HartJohnstoneForestClassesMap, 
    currentCohortAgeMap
    ),
    fun = function(x) predict_lichen_biomass(x[1], x[2], x[3])
  )

  # Assign it a name 
  names(sim$WB_LichenBiomassMap) <- "lichenBiomass"
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  userTags <- c(currentModule(sim), "function:.inputObjects")
  ##############################################################################
  # Generate a fake WB_HartJohnstoneForestClassesMap if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("WB_HartJohnstoneForestClassesMap", sim)){
    rastWidth <- 1000
    message("##############################################################################")   
    message("WB_HartJohnstoneForestClassesMap not supplied.")   
    message("Please couple with the WB_HartJohnstoneForestClasses module. ")   
    message("Creating random map ", rastWidth, " pixels by ", rastWidth, " pixels for 6 forest ")   
    message("classes (\"deci (1)\", \"mixed (2)\", \"conimix (3)\", \"jackpine (4)\", ")   
    message("\"larch (5)\" and \"spruce (6)\")...")

    sim$WB_HartJohnstoneForestClassesMap <- Cache(
      getRandomCategoricalMap,
      origin = c(-667296, 1758502),
      ncol = rastWidth,
      nrow = rastWidth,
      crs = "ESRI:102002",
      nbregion = 2000,
      valuevect = 1:6,
      seed = 100,
      userTags = c(userTags, "WB_HartJohnstoneForestClassesMap")
    )
    
    # Convert to factor and add proper labels
    sim$WB_HartJohnstoneForestClassesMap <- terra::as.factor(sim$WB_HartJohnstoneForestClassesMap)
    levels(sim$WB_HartJohnstoneForestClassesMap) <- data.frame(
      value = c(1L, 2L, 3L, 4L, 5L, 6L),
      class = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce")
    )
    names(sim$WB_HartJohnstoneForestClassesMap) <- "standtype"
  }

  if (!is.null(sim$WB_HartJohnstoneForestClassesMap)){
    baseRast <- sim$WB_HartJohnstoneForestClassesMap
  }
  else if (!is.null(sim$pixelGroupMap)){
    baseRast <- sim$pixelGroupMap
  }
  else if (!is.null(sim$rasterToMatch)){
    baseRast <- sim$rasterToMatch
  }
  else {
    stop(paste("At least one of WB_HartJohnstoneForestClassesMap, pixelGroupMap or ",
               "rasterToMatch must be defined in sim before WB_LichenBiomass can be initialized..."))
  }
  # baseExtent <- ext(baseRast)
  # baseCRS <- crs(baseRast)
  
  ##############################################################################
  # Create a dummy standAgeMap
  ##############################################################################
  if(!suppliedElsewhere("standAgeMap", sim)){
    message("##############################################################################")   
    message("standAgeMap not supplied.")   
    message("Please couple with the Biomass_core module which simulate stand age...")  
    message("Creating random map ", rastWidth, " pixels by ", rastWidth, " pixels ")   
    message("with ages from 0 to 150...")
    sim$standAgeMap <- Cache(
      getRandomCategoricalMap,
      # origin = c(-667296, 1758502),
      origin = c(xmin(baseRast), ymin(baseRast)),
      ncol = ncol(baseRast),
      nrow = nrow(baseRast),
      crs = crs(baseRast),
      nbregion = 2000,
      valuevect = 0:150,
      seed = 100,
      userTags = c(userTags, "standAgeMap")
    )
    names(sim$standAgeMap) <- "standage"
  }
  
  ##############################################################################
  # Download, process and cache ecoprovince if it is not supplied
  # https://www.epa.gov/eco-research/ecoregions-north-america
  ##############################################################################
  if(!suppliedElsewhere("EcoProvincesMap", sim)){
    message("##############################################################################")   
    message("EcoProvincesMap not supplied.")   
    message("Downloading and projecting raster NAs with SoilGrids values...")
    ecoProv <- Cache(
      prepInputs,
      url = extractURL("EcoProvincesMap", sim),
      targetFile = "NA_CEC_Eco_Level3.shp",
      destinationPath = getPaths()$cache,
      projectTo = baseRast,
      cropTo = vect(ext(baseRast), crs(baseRast)),
      writeTo = file.path(getPaths()$cache, "NA_CEC_Eco_Level3_postProcessed.shp"),
      fun = terra::vect,
      userTags = c(userTags, "NA_CEC_Eco_Level3_postProcessed.shp"),
      overwrite = TRUE
    )

    ecoProv <- ecoProv[, c("NA_L3NAME")]
    ecoProv$NA_L3NAME <- as.factor(ecoProv$NA_L3NAME)
    names(ecoProv)[names(ecoProv) == "NA_L3NAME"] <- "ecoprov"
    sim$EcoProvincesMap <- rasterize(ecoProv, baseRast, field = "ecoprov")
    
    # Alternative dataset from Canada Open
    # # https://open.canada.ca/data/en/dataset/98fa7335-fbfe-4289-9a0e-d6bf3874b424
    # sim$EcoProvincesMap <- Cache(
    #   prepInputs,
    #   url = "https://agriculture.canada.ca/atlas/data_donnees/nationalEcologicalFramework/data_donnees/geoJSON/ep/nef_ca_ter_ecoprovince_v2_2.geojson",
    #   targetFile = "nef_ca_ter_ecoprovince_v2_2.geojson",
    #   destinationPath = getPaths()$cache,
    #   projectTo = plotAndPixelGroupAreaRast,
    #   cropTo = plotAndPixelGroupArea,
    #   writeTo = file.path(getPaths()$cache, paste0("nef_ca_ter_ecoprovince_v2_2_postProcessed.geojson")),
    #   fun = terra::vect
    # )
    # sim$EcoProvincesMap <- sim$EcoProvincesMap[, c("ECOPROVINCE_NAME_EN")]
    # sim$EcoProvincesMap$ECOPROVINCE_NAME_EN <- as.factor(sim$EcoProvincesMap$ECOPROVINCE_NAME_EN)
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
