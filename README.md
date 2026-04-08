Code and data associated with the article "Rivers of salt: the extent of
salinization in European running waters" by Crabot et al.

The code/ folder contains three files:
- 1_calculating_salinization.R: Loading data on salinity and environment, defining
reference sites, using these sites to predict baseline salinity with a random
forest model, calculating salinization as the difference between observed and
predicted baseline salinity.
- 2_finding_salinization_drivers.R: computing linear models and spatial 
autoregressive model of salinization in function of the environmental predictors
- 3_Gibbs_plot.R: plotting salinity in function of the Na / (Na + Ca) ratio

The data/ folder contains 8 files and two subfolders:
- coastal_catchments.csv: vector of ID of catchments that are along the coast
- environment.csv: combines all environmental variables that stem from
RiverATLAS, CHELSA, FutureStreams, LUCAS LUC, GIRES.
- geology_upstream_12catchments_percentage.csv: total % of the upstream area
covered by each geological class (as in IGME5000) for all level-12 HydroBASINS
catchments
- hyriv_hybas_equivalence.csv: equivalence between ID of HydroRIVERS
reaches and ID of catchments of HydroBASINS
- Ratios_EC_Cl_Na_Ca_time.csv: table with ion concentrations in Cl, Na, Ca
- salty_temporal_EC.csv: salinity data, stemming from the GlobSalt dataset
- seawater_intrusion_basin12.csv: total area of seawater intrusion in level-12
BasinATLAS catchments
- to_the_sea.csv: distance to the sea in kilometers for all river reaches
- hydroatlas/ folder: files of HydroRIVERS and HydroBASINS (8th-level catchment),
   downloaded from https://www.hydrosheds.org/hydroatlas
- output/ folder: folder where the output of the scripts will be written
