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
    expectsInput("cohortData",  "data.table",
                 desc = paste("Initial community table, created from available biomass (g/m2)",
                              "age and species cover data, as well as ecozonation information",
                              "Columns: B, pixelGroup, speciesCode")),
    expectsInput("pixelGroupMap", "SpatRast",
                 desc = "Initial community map that has mapcodes match initial community table")
    
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
      # We initiate a fake cohortData here because it is dependent on sim$pixelGroupMap
      # We do that only if it is not supplied by another module (like Biomass_core)
      # We have also set loadOrder to "after Biomass_core" so pixelGroupMap should be initialized
      if (!suppliedElsewhere(sim$cohortData) && !is.null(sim$pixelGroupMap)) {
        nbGroup <- length(unique(values(sim$pixelGroupMap)))
        message("##############################################################################")   
        message("cohortData not supplied.")   
        message("Please provide one. Generating random cohort data for ", nbGroup, " pixel groups...")
        sim$cohortData <- getRandomCohortData(nbPixelGroup = nbGroup, 
                                              pixelSize = res(sim$pixelGroupMap)[1])
      }
      
      if(!suppliedElsewhere("WB_HartJohnstoneForestClassesMap", sim) && !is.null(sim$cohortData) && !is.null(sim$pixelGroupMap)){
        source("https://raw.githubusercontent.com/pedrogit/WB_HartJohnstoneForestClasses/refs/heads/main/R/WB_HartJohnstoneForestClasses.r")
        sim$WB_HartJohnstoneForestClassesMap <- classifyStand(
          cohortData = sim$cohortData, 
          pixelGroupMap = sim$pixelGroupMap,
          time = time(sim)
        )
      }
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
# TSSRF (Time Since Stand-Replacing Fire) is approximated by the sim$cohortData age column
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
  message("Recomputing sim$WB_LichenBiomassMap for ", 
          format(ncell(sim$WB_HartJohnstoneForestClassesMap), scientific = FALSE), " pixels..")

  # Get a current age map from cohortData 
  currentCohortAgeMap <- LandR::standAgeMapGenerator(
    sim$cohortData,
    sim$pixelGroupMap
  )
  names(currentCohortAgeMap) <- "current_age"
  # varnames(currentCohortAgeMap) <- "current_age"
  
  
  # This is very slow. Try instead to isolate the unique values, compute the 
  # biomass for those values and use rasterizeReduced to create the map.
    currentCohortAgeMap

  combined_rast <- currentCohortAgeMap * 10000 + sim$EcoProvincesMap * 10 + sim$WB_HartJohnstoneForestClassesMap
  names(combined_rast) <- "combined_val"
  combined_rast_unique_values <- unique(values(combined_rast))
  age <- combined_rast_unique_values %/% 10000
  eco_prov <- (combined_rast_unique_values - age * 10000) %/% 10
  for_class <- combined_rast_unique_values - age* 10000 - eco_prov * 10
  biomass <- data.table(
    combined_val = combined_rast_unique_values, 
    biomass = predict_lichen_biomass(eco_prov, for_class, age)
  )
  names(biomass) <- c("combined_val", "biomass")
  
  # Biomass is in kg/ha. We have to divide by 10000 and multiply by the pixel size
  biomass$biomass <- biomass$biomass / 10000 * prod(res(sim$pixelGroupMap))
  
  sim$WB_LichenBiomassMap <- SpaDES.tools::rasterizeReduced(
    reduced = biomass,
    fullRaster = combined_rast,
    mapcode = "combined_val", 
    newRasterCols = "biomass"
  )
  
  # sim$WB_LichenBiomassMap <- app(c(
  #   sim$EcoProvincesMap, 
  #   sim$WB_HartJohnstoneForestClassesMap, 
  #   currentCohortAgeMap
  #   ),
  #   fun = function(x) predict_lichen_biomass(x[1], x[2], x[3])
  # )
  
  # Compute the biomass for the non forested area using reclass
  # We must convert mean_biomass per ha to the area of a pixel.
  # e.g. if the table says 2kg/ha and one pixel is 0.25 ha then we must divide by 4
  # mean_biomass * pixel
  
  # Assign it a name 
  names(sim$WB_LichenBiomassMap) <- "lichenBiomass"
  # varnames(sim$WB_LichenBiomassMap) <- "lichenBiomass"
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  userTags <- c(currentModule(sim), "function:.inputObjects")
  if(!suppliedElsewhere("pixelGroupMap", sim)){
    nbGroup <- 200
    pixelGroupRastWidth <- 1000
    message("##############################################################################")   
    message("pixelGrouMap not supplied.")   
    message("Please provide one. Creating random map ", pixelGroupRastWidth, " pixels by ", 
            pixelGroupRastWidth, " pixels with ", nbGroup, " groups...")
    
    sim$pixelGroupMap <- Cache(
      getRandomCategoricalMap,
      origin = c(-667296, 1758502),
      ncol = pixelGroupRastWidth,
      nrow = pixelGroupRastWidth,
      crs = "ESRI:102002",
      nbregion = nbGroup,
      seed = 100,
      userTags = c(userTags, "pixelGroupMap"), 
      omitArgs = c("userTags")
    )
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
