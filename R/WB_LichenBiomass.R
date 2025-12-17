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

computeLichenBiomassMap <- function(
    cohortData,
    pixelGroupMap,
    ecoProvincesMap,
    HJForestclassesMap
){
  # Get a current age map from cohortData 
  currentCohortAgeMap <- LandR::standAgeMapGenerator(
    cohortData,
    pixelGroupMap
  )
  names(currentCohortAgeMap) <- "current_age"
  # varnames(currentCohortAgeMap) <- "current_age"
  
  # Isolate the unique values, compute the biomass for those values
  # and use rasterizeReduced to create the map.
  combined_rast <- currentCohortAgeMap * 10000 + ecoProvincesMap * 10 + HJForestclassesMap
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
  biomass$biomass <- biomass$biomass / 10000 * prod(res(pixelGroupMap))
  
  lichenBiomassMap <- SpaDES.tools::rasterizeReduced(
    reduced = biomass,
    fullRaster = combined_rast,
    mapcode = "combined_val", 
    newRasterCols = "biomass"
  )
  
  # slower code
  # lichenBiomassMap <- app(c(
  #   ecoProvincesMap, 
  #   HJForestClassesMap, 
  #   currentCohortAgeMap
  #   ),
  #   fun = function(x) predict_lichen_biomass(x[1], x[2], x[3])
  # )
  
  # Compute the biomass for the non forested area using reclass
  # We must convert mean_biomass per ha to the area of a pixel.
  # e.g. if the table says 2kg/ha and one pixel is 0.25 ha then we must divide by 4
  # mean_biomass * pixel
  
  # Assign it a name 
  names(lichenBiomassMap) <- "lichenBiomass"
  # varnames(sim$WB_LichenBiomassMap) <- "lichenBiomass"
  return(lichenBiomassMap)
}