defineModule(sim, list(
  name = "WB_LichenBiomass",
  description = paste("Compute lichen biomass according to the Greuel, Degre-Timmons model (2021)"),
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
                              "standtype."),
                 sourceURL = NA),
    expectsInput(objectName = "EcoProvincesMap", 
                 objectClass = "SpatVector", 
                 desc = "SpatRaster of north american ecoprovinces.", 
                 sourceURL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip"),
    expectsInput("cohortData",  "data.table",
                 desc = paste("Initial community table, created from available ",
                              "biomass (g/m2) age and species cover data, as well ",
                              "as ecozonation information",
                              "Columns: B, pixelGroup, speciesCode.")),
    expectsInput("pixelGroupMap", "SpatRast",
                 desc = "Initial community map that has mapcodes match initial ",
                        "community table"),
    expectsInput("WB_MeanBiomassPerVegClasses", "data.table",
                 desc = "CSV table of static values for mean biomass per MVI and ",
                        "some HartJohnstone forest classes")
    
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
      sim <- Init(sim)
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

Init <- function(sim){
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
  return(invisible(sim))
}

reComputeLichenBiomassMap <- function(sim) {
  message("Recomputing sim$WB_LichenBiomassMap for ", 
          format(ncell(sim$WB_HartJohnstoneForestClassesMap), scientific = FALSE), " pixels..")

  sim$WB_LichenBiomassMap <- computeLichenBiomassMap(
    sim$cohortData,
    sim$pixelGroupMap,
    sim$EcoProvincesMap,
    sim$WB_HartJohnstoneForestClassesMap,
    sim$WB_NonForestedVegClassesMap,
    sim$WB_MeanBiomassPerVegClasses
  )
  
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
      userTags = c(userTags, "WB_pixelGroupMap"), 
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

  ##############################################################################
  # Load the supplied biomass table if an alternative one is not provided by the 
  # user. 
  ##############################################################################
  if(!suppliedElsewhere("WB_MeanBiomassPerVegClasses")){
    sim$WB_MeanBiomassPerVegClasses <- fread(file.path(dataPath(sim), "meanBiomassPerVegClasses.csv"))
  }
  return(invisible(sim))
}
