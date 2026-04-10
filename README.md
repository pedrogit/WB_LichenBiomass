# WB_LichenBiomass

## Introduction

WB_LichenBiomass is a [SpaDES](https://spades.predictiveecology.org/) module complementing the [LandR](https://landr-manual.predictiveecology.org/) ecosystem of modules for forest biomass and succession simulation. It is part of an ensemble of modules that provide to LandR the statistical prediction of terrestrial lichen biomass from stand type, time-since fire, and terrestrial ecoprovince.

These modules are an implementatiom of [Greuel and Degré-Timmons et al (2021)](https://esajournals-onlinelibrary-wiley-com.acces.bibl.ulaval.ca/doi/full/10.1002/ecs2.3481), developed to support lichen biomass modelling for woodland caribou conservation in the Northwest Territories. The geographical area wherein the model may reasonably be applied should be assessed from Figure 2 of the cited paper. 

The components of the module ensemble are:

- [WB_HartJohnstoneForestClasses](https://github.com/pedrogit/WB_HartJohnstoneForestClasses) - Generates a map classifying LandR forested pixels to 6 (or 7) classes.
- [WB_VegBasedDrainage](https://github.com/pedrogit/WB_VegBasedDrainage) - Generates a map of two drainage classes.
- [WB_NonForestedVegClasses](https://github.com/pedrogit/WB_NonForestedVegClasses) - Generates a map of land cover classes for areas LandR considers to be non-forested.
- [WB_LichenBiomass](https://github.com/pedrogit/WB_LichenBiomass) - This module. Generates a wall-to-wall map of predicted lichen biomass density for forested and non-forested pixels.

These modules are derived from extensive empirical research in the northwest boreal of North America, as described in [Greuel and Degré-Timmons et al (2021)](https://esajournals-onlinelibrary-wiley-com.acces.bibl.ulaval.ca/doi/full/10.1002/ecs2.3481), Casheiro-Guilhem et. al (in prep.) and foundational papers by [Hart and Johnstone et al. (2018)](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.14550).

## WB_LichenBiomass Module Overview

At each simulation step, WB_LichenBiomass generate a raster map of lichen density (kg/ha) using a statistical model for forested areas and a reference table of observed lichen biomass for different classes of non-forested areas. The statistical model takes ecoprovince, forest type (from the WB_HartJohnstoneForestClasses module) and age (provided by the biomass_core module) into account. The reference table used for the non-forested area takes only the vegetation classes (provided by the WB_NonForestedVegClasses module) into account.

Areas not covered by Biomass_core and the land cover raster are set to NA.

The module time step would normally be the same 10 year period used by LandR.

WB_LichenBiomass is dynamic because it depends on inputs by the biomass_core module, the WB_HartJohnstoneForestClasses module and the WB_NonForestedVegClasses module which are themself dynamic.

### Authors and Citation

* Pierre Racine <pierre.racine@sbf.ulaval.ca> [aut, cre]
* Andres Caseiro Guilhem <andres.caseiro-guilhem.1@ulaval.ca> [aut]
* Steven G. Cumming <stevec.boreal@gmail.com> [aut]
* Ruth J. Greuel
* Geneviève É. Degré-Timmons

Racine, P., Caseiro Guilhem, A., Cumming, S.G., Greuel, Ruth J., Degré-Timmons, Geneviève É. (2026) *WB_LichenBiomass: A SpaDES module to compute lichen biomass density in western boreal forests of Canada.* SpaDES Module.

### Module Parameters

| Parameter | Class | Default | Description |
| --- | --- | --- | --- |
| WB_LichenBiomassTimeStep | integer | 10 | Module return interval, at which the SpatRaster is regenerated. |


### Expected Module Inputs

| Input Object | Class | Description |
| --- | --- | --- |
| WB_HartJohnstoneForestClassesMap | SpatRast | Forest classification map produced by the WB_HartJohnstoneForestClasses modules or equivalent. |
| EcoProvincesMap | SpatVector | SpatVector of the North American ecoprovinces level III or equivalent. |
| cohortData | data.table | LandD main state of the vegetation table or equivalent. |
| pixelGroupMap | SpatRast | LandD raster of active forested pixel groups. |
| WB_MeanBiomassPerVegClasses | data.table | Table of lichen biomass by vegetation classes used for non-forested areas. |

### Module Outputs

| Output Object | Class | Description |
| --- | --- | --- |
| WB_LichenBiomassMap | SpatRast | Raster map lichen biomass densities (kg/ha). |

### Code

The code is available here: https://github.com/pedrogit/WB_LichenBiomass

### Minimal Self Contained Workflow Example

Soon...



