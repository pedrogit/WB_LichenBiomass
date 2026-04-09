# ---------------------------------------------------
# 1. Model coefficients from Greuel et al. (2021)
# ---------------------------------------------------

# Logistic regression for presence/absence (lichen_prob)
lichen_prob_coefs <- list(
  Intercept = -0.500169,

  # standtype main effects
  standtypeDeciduous = -0.801061,
  standtypeDeciduous_conifer_mix = -0.380172,
  standtypeJack_pine = 0.953791,
  standtypeLarch = -15.068047,
  standtypePoorly_drained_spruce = 0.395372,
  standtypeWell_drained_spruce = -0.587152,

  # ecoprovince main effects
  ecoprovinceCoppermine = -0.398703,
  ecoprovinceGreat_Bear_Plains = -0.602398,
  ecoprovinceHay_Slave_River = -0.724421,

  # interactions
  ecoprovinceCoppermine_standtypeDeciduous = 0.969629,
  ecoprovinceGreat_Bear_Plains_standtypeDeciduous = -13.243396,
  ecoprovinceHay_Slave_River_standtypeDeciduous = 0.478505,

  ecoprovinceCoppermine_standtypeDeciduous_conifer_mix = 0.434797,
  ecoprovinceGreat_Bear_Plains_standtypeDeciduous_conifer_mix = -13.664284,
  ecoprovinceHay_Slave_River_standtypeDeciduous_conifer_mix = -0.533582,

  ecoprovinceCoppermine_standtypeJack_pine = 1.595913,
  ecoprovinceGreat_Bear_Plains_standtypeJack_pine = -0.074130,
  ecoprovinceHay_Slave_River_standtypeJack_pine = -0.591750,

  ecoprovinceHay_Slave_River_standtypeLarch = 13.223912,

  ecoprovinceCoppermine_standtypePoorly_drained_spruce = -0.533398,
  ecoprovinceGreat_Bear_Plains_standtypePoorly_drained_spruce = 0.818873,
  ecoprovinceHay_Slave_River_standtypePoorly_drained_spruce = -0.318033,

  ecoprovinceCoppermine_standtypeWell_drained_spruce = 0.929588,
  ecoprovinceGreat_Bear_Plains_standtypeWell_drained_spruce = 1.393756,
  ecoprovinceHay_Slave_River_standtypeWell_drained_spruce = 0.624810,

  TSSRF = 0.014524
)

# ---------------------------------------------------
# Nonlinear biomass model
# ---------------------------------------------------

predicted_biomass_alpha <- list(

  Intercept = 912.773075,

  standtypeJack_pine = 519.237824,
  standtypePoorly_drained_spruce = 1379.215795,
  standtypeWell_drained_spruce = 178.352000,

  ecoprovinceCoppermine = 72.311354,
  ecoprovinceGreat_Bear_Plains = -382.594163,
  ecoprovinceHay_Slave_River = -390.314399,

  ecoprovinceCoppermine_standtypeJack_pine = -484.423470,
  ecoprovinceCoppermine_standtypePoorly_drained_spruce = -1531.744131,
  ecoprovinceCoppermine_standtypeWell_drained_spruce = -461.555551,

  ecoprovinceGreat_Bear_Plains_standtypeJack_pine = 2201.925623,
  ecoprovinceGreat_Bear_Plains_standtypePoorly_drained_spruce = 693.948820,
  ecoprovinceGreat_Bear_Plains_standtypeWell_drained_spruce = 639.592595,

  ecoprovinceHay_Slave_River_standtypeJack_pine = -118.815728,
  ecoprovinceHay_Slave_River_standtypePoorly_drained_spruce = -543.170075,
  ecoprovinceHay_Slave_River_standtypeWell_drained_spruce = 232.210759
)

predicted_biomass_xmid <- list(
  Intercept = 28.041319,
  ecoprovinceCoppermine = 17.202108,
  ecoprovinceGreat_Bear_Plains = 44.979595,
  ecoprovinceHay_Slave_River = 41.765281
)

predicted_biomass_scal <- 3.0647


# ---------------------------------------------------
# 2. Prediction function
# ---------------------------------------------------

predict_lichen_biomass <- function(ecoprov, standtype, TSSRF){

  # -----------------------------
  # Logistic regression (presence)
  # -----------------------------
  logit <- lichen_prob_coefs$Intercept +

    ifelse(standtype==1, lichen_prob_coefs$standtypeDeciduous, 0) +
    ifelse(standtype==2, lichen_prob_coefs$standtypeDeciduous_conifer_mix, 0) +
    ifelse(standtype==4, lichen_prob_coefs$standtypeJack_pine, 0) +
    ifelse(standtype==5, lichen_prob_coefs$standtypeLarch, 0) +
    ifelse(standtype==6, lichen_prob_coefs$standtypeWell_drained_spruce, 0) +
    ifelse(standtype==7, lichen_prob_coefs$standtypePoorly_drained_spruce, 0) +

    ifelse(ecoprov==10, lichen_prob_coefs$ecoprovinceCoppermine, 0) +
    ifelse(ecoprov==11, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains, 0) +
    ifelse(ecoprov==12, lichen_prob_coefs$ecoprovinceHay_Slave_River, 0) +

    ifelse(ecoprov==10 & standtype==1, lichen_prob_coefs$ecoprovinceCoppermine_standtypeDeciduous,0) +
    ifelse(ecoprov==11 & standtype==1, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains_standtypeDeciduous,0) +
    ifelse(ecoprov==12 & standtype==1, lichen_prob_coefs$ecoprovinceHay_Slave_River_standtypeDeciduous,0) +

    ifelse(ecoprov==10 & standtype==2, lichen_prob_coefs$ecoprovinceCoppermine_standtypeDeciduous_conifer_mix,0) +
    ifelse(ecoprov==11 & standtype==2, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains_standtypeDeciduous_conifer_mix,0) +
    ifelse(ecoprov==12 & standtype==2, lichen_prob_coefs$ecoprovinceHay_Slave_River_standtypeDeciduous_conifer_mix,0) +

    ifelse(ecoprov==10 & standtype==4, lichen_prob_coefs$ecoprovinceCoppermine_standtypeJack_pine,0) +
    ifelse(ecoprov==11 & standtype==4, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains_standtypeJack_pine,0) +
    ifelse(ecoprov==12 & standtype==4, lichen_prob_coefs$ecoprovinceHay_Slave_River_standtypeJack_pine,0) +

    ifelse(ecoprov==12 & standtype==5, lichen_prob_coefs$ecoprovinceHay_Slave_River_standtypeLarch,0) +

    ifelse(ecoprov==10 & standtype==7, lichen_prob_coefs$ecoprovinceCoppermine_standtypePoorly_drained_spruce,0) +
    ifelse(ecoprov==11 & standtype==7, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains_standtypePoorly_drained_spruce,0) +
    ifelse(ecoprov==12 & standtype==7, lichen_prob_coefs$ecoprovinceHay_Slave_River_standtypePoorly_drained_spruce,0) +

    ifelse(ecoprov==10 & standtype==6, lichen_prob_coefs$ecoprovinceCoppermine_standtypeWell_drained_spruce,0) +
    ifelse(ecoprov==11 & standtype==6, lichen_prob_coefs$ecoprovinceGreat_Bear_Plains_standtypeWell_drained_spruce,0) +
    ifelse(ecoprov==12 & standtype==6, lichen_prob_coefs$ecoprovinceHay_Slave_River_standtypeWell_drained_spruce,0) +

    lichen_prob_coefs$TSSRF * TSSRF

  p_presence <- 1/(1+exp(-logit))


  # -----------------------------
  # Nonlinear biomass
  # -----------------------------

  valid <- standtype %in% c(3,4,6,7)

  alpha <- rep(NA,length(standtype))
  xmid  <- rep(NA,length(standtype))

  alpha[valid] <-

    predicted_biomass_alpha$Intercept +

    ifelse(standtype[valid]==4, predicted_biomass_alpha$standtypeJack_pine,0) +
    ifelse(standtype[valid]==6, predicted_biomass_alpha$standtypeWell_drained_spruce,0) +
    ifelse(standtype[valid]==7, predicted_biomass_alpha$standtypePoorly_drained_spruce,0) +

    ifelse(ecoprov[valid]==10, predicted_biomass_alpha$ecoprovinceCoppermine,0) +
    ifelse(ecoprov[valid]==11, predicted_biomass_alpha$ecoprovinceGreat_Bear_Plains,0) +
    ifelse(ecoprov[valid]==12, predicted_biomass_alpha$ecoprovinceHay_Slave_River,0) +

    ifelse(ecoprov[valid]==10 & standtype[valid]==4,
           predicted_biomass_alpha$ecoprovinceCoppermine_standtypeJack_pine,0) +

    ifelse(ecoprov[valid]==10 & standtype[valid]==7,
           predicted_biomass_alpha$ecoprovinceCoppermine_standtypePoorly_drained_spruce,0) +

    ifelse(ecoprov[valid]==10 & standtype[valid]==6,
           predicted_biomass_alpha$ecoprovinceCoppermine_standtypeWell_drained_spruce,0) +

    ifelse(ecoprov[valid]==11 & standtype[valid]==4,
           predicted_biomass_alpha$ecoprovinceGreat_Bear_Plains_standtypeJack_pine,0) +

    ifelse(ecoprov[valid]==11 & standtype[valid]==7,
           predicted_biomass_alpha$ecoprovinceGreat_Bear_Plains_standtypePoorly_drained_spruce,0) +

    ifelse(ecoprov[valid]==11 & standtype[valid]==6,
           predicted_biomass_alpha$ecoprovinceGreat_Bear_Plains_standtypeWell_drained_spruce,0) +

    ifelse(ecoprov[valid]==12 & standtype[valid]==4,
           predicted_biomass_alpha$ecoprovinceHay_Slave_River_standtypeJack_pine,0) +

    ifelse(ecoprov[valid]==12 & standtype[valid]==7,
           predicted_biomass_alpha$ecoprovinceHay_Slave_River_standtypePoorly_drained_spruce,0) +

    ifelse(ecoprov[valid]==12 & standtype[valid]==6,
           predicted_biomass_alpha$ecoprovinceHay_Slave_River_standtypeWell_drained_spruce,0)


  xmid[valid] <- predicted_biomass_xmid$Intercept +
    ifelse(ecoprov[valid]==10, predicted_biomass_xmid$ecoprovinceCoppermine,0) +
    ifelse(ecoprov[valid]==11, predicted_biomass_xmid$ecoprovinceGreat_Bear_Plains,0) +
    ifelse(ecoprov[valid]==12, predicted_biomass_xmid$ecoprovinceHay_Slave_River,0)


  biomass <- rep(NA,length(standtype))

  biomass[valid] <- alpha[valid] /
    (1 + exp((xmid[valid] - TSSRF[valid]) / predicted_biomass_scal))


  # -----------------------------
  # Expected biomass
  # -----------------------------

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
  # uniqueValsDT$biomass <- uniqueValsDT$biomass / 10000 * prod(res(pixelGroupMap))
  
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