# ---------------------------------------------------
# 1. Model coefficients from Greuel et al. (2021)
# ---------------------------------------------------
# Logistic regression for presence/absence (lichen_prob)
lichen_prob_coefs <- list(
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

# Non-linear biomass model (predicted_biomass)
predicted_biomass_alpha <- list(
  Intercept = 912.7731,              # baseline = "conimix" (level 3)
  standtypeJack_pine = 519.2378,     # level 4
  standtypePoorly_drained_spruce = 1379.2158, # level 7
  standtypeWell_drained_spruce = 178.3520     # level 6
)
predicted_biomass_xmid <- list(
  Intercept = 28.0413,
  ecoprovinceCoppermine = 17.2021,
  ecoprovinceGreat_Bear_Plains = 44.9796,
  ecoprovinceHay_Slave_River = 41.7653
)
predicted_biomass_scal <- 3.0647

# Function predicting lichen biomass (kg/ha) from ecoprovince, stand type and TSSRF vectors.
# TSSRF (Time Since Stand-Replacing Fire) is approximated by the cohortData age column.
predict_lichen_biomass <- function(ecoprov, standtype, TSSRF) {
  # --- Logistic regression (presence probability) ---
  logit <- lichen_prob_coefs$Intercept +
    ifelse(standtype == 1, lichen_prob_coefs$standtypeDeciduous, 0) +
    ifelse(standtype == 2, lichen_prob_coefs$standtypeDeciduous_conifer_mix, 0) +
    ifelse(standtype == 4, lichen_prob_coefs$standtypeJack_pine, 0) +
    ifelse(standtype == 5, lichen_prob_coefs$standtypeLarch, 0) +
    ifelse(standtype == 6, lichen_prob_coefs$standtypeWell_drained_spruce, 0) +
    ifelse(standtype == 7, lichen_prob_coefs$standtypePoorly_drained_spruce, 0) +
    ifelse(ecoprov == 0, 0,
           ifelse(ecoprov == 1, lichen_prob_coefs$ecoprovinceCoppermine,
                  ifelse(ecoprov == 2, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains,
                         lichen_prob_coefs$ecoprovinceHay_Slave_River))) +
    lichen_prob_coefs$TSSRF * TSSRF
  
  p_presence <- 1 / (1 + exp(-logit))
  
  # --- Nonlinear biomass (only for valid stands) ---
  valid <- standtype %in% c(3, 4, 6, 7)  # conimix, jackpine, wd_spruce, pd_spruce
  
  alpha <- rep(NA, length(standtype))
  xmid  <- rep(NA, length(standtype))
  
  alpha[valid] <- predicted_biomass_alpha$Intercept +
    ifelse(standtype[valid] == 4, predicted_biomass_alpha$standtypeJack_pine, 0) +
    ifelse(standtype[valid] == 6, predicted_biomass_alpha$standtypeWell_drained_spruce, 0) +
    ifelse(standtype[valid] == 7, predicted_biomass_alpha$standtypePoorly_drained_spruce, 0)
  
  xmid[valid] <- predicted_biomass_xmid$Intercept +
    ifelse(ecoprov[valid] == 1, predicted_biomass_xmid$ecoprovinceCoppermine,
           ifelse(ecoprov[valid] == 2, predicted_biomass_xmid$ecoprovinceGreat_Bear_Plains,
                  predicted_biomass_xmid$ecoprovinceHay_Slave_River))
  
  biomass <- rep(NA, length(standtype))
  biomass[valid] <- alpha[valid] / (1 + exp((xmid[valid] - TSSRF[valid]) / predicted_biomass_scal))
  
  # --- Combine (expected biomass = probability × biomass) ---
  expected_biomass <- p_presence * biomass
  
  return(expected_biomass)
}

computeLichenBiomassMap <- function(
    cohortData,
    pixelGroupMap,
    ecoProvincesMap,
    HJForestclassesMap,
    nonForestedVegClassesMap,
    biomassMeansPerVegClassesDT
){
  # Get a current age map from cohortData 
  currentCohortAgeMap <- LandR::standAgeMapGenerator(
    cohortData,
    pixelGroupMap
  )
  names(currentCohortAgeMap) <- "current_age"
  
  # Merge the HJForestclassesMap and the nonForestedVegClassesMap together 
  # so we have to deal with only one raster. HJForestclassesMap have precedence
  # over nonForestedVegClassesMap.
  mergedVegMap <- cover(HJForestclassesMap, nonForestedVegClassesMap)

  # Isolate the rasters stack unique values
  stackedMap <- c(currentCohortAgeMap, ecoProvincesMap, mergedVegMap)
  
  # Create a unique ID for each combination of unique values
  allVals <- as.data.table(values(stackedMap))
  allVals[, id := .GRP, by = .(current_age, ecoprov, standtype)]
  
  # Assign the IDs back to a new raster
  IDMap <- setValues(currentCohortAgeMap, allVals$id)
  names(IDMap) <- c("id")
  
  # Add the Id raster to the stack
  stackedMap <- c(stackedMap, IDMap)
  # Determine the unique values for which to compute the biomass values 
  # (this is much faster than computing the biomass for every raster pixels)
  uniqueValsDT <- as.data.table(unique(values(stackedMap)))
  names(uniqueValsDT) <- c("age", "eco_prov", "veg_class", "id")
  
  # Compute the biomass for forested areas
  uniqueValsDT[, biomass := predict_lichen_biomass(eco_prov, veg_class, age)]

  # The model does not predict values for some forest 1 (deci), 2 (mixed) and 
  # 5 (larch) and vegetabion classes, we take them from the provided table.
  uniqueValsDT[biomassMeansPerVegClassesDT, on = "veg_class", biomass := i.mean_lichen_biomass]

# browser()
  # uniqueValsDT[order(veg_class, age), .(veg_class, age, biomass = format(biomass, scientific = FALSE))]

  # Biomass is in kg/ha. We have to divide by 10000 and multiply by the pixel size
  uniqueValsDT$biomass <- uniqueValsDT$biomass / 10000 * prod(res(pixelGroupMap))
  
  # Use rasterizeReduced to create the map.
  lichenBiomassMap <- SpaDES.tools::rasterizeReduced(
    reduced = uniqueValsDT,
    fullRaster = IDMap,
    mapcode = "id", 
    newRasterCols = "biomass"
  )

  # Assign it a name 
  names(lichenBiomassMap) <- "lichenBiomass"

  return(lichenBiomassMap)
}