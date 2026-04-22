#### Setting anthropogenic pressures to zero ####

# parameters for the RF
set.seed(17)
num.trees = 1000
# mtry = number of variables to split at each node
mtry = 7 # default = square roots of nb of predictors, would be 11 here
# min node size = number of individuals at which the RF is not splitting the node anymore
min.node.size = 5 # 5 = default for regression 

# RF model of conductivity on all sites with observations
model_all <- ranger(EU_salt$EC ~ .,
                   data = EU_salt %>% select(alkali_basalt:predprob30), # without random factor = necessary to do predictions
                   num.trees = num.trees,
                   mtry= mtry,
                   replace = T,
                   oob.error = T,
                   keep.inbag = T,
                   num.threads = 4,
                   importance ='impurity',
                   probability = FALSE)

# Characteristics of the model
var_imp_all <- importance(model_all)[order(importance(model_all),decreasing=TRUE)]
r.squared_all <- model_all$r.squared
to_the_sea_index <- which(names(var_imp_all) == "to_the_sea")
seawater_index <- which(names(var_imp_all) == "seawater_area")

var_imp_all # Variable importance rank
r.squared_all # R squared
to_the_sea_index # Variable importance of distance to the sea
seawater_index # Variable importance of seawater intrusion

# Plotting variable importance
plot(var_imp_all, type = "n", ylab= "Variable importance", main = "Variable importance in predicting conductivity of allerence sites")
text(seq(1, 114, 1), var_imp_all, labels=names(var_imp_all), cex= 0.6, pos=3)


#### Building salinity map ####

if (! ("rivers" %in% ls(envir = .GlobalEnv))) { # checks if the file is not already loaded to save time
  
  # River network
  rivers <-terra::vect(paste0(wd, "/hydroatlas/HydroRIVERS_v10_eu.shp"),
                       extent= c(-10.5,32.702, 34.856, 71.31)) %>% # already cropped to the right extent
    select(HYRIV_ID) # only keeping the ID column
  
  # Catchments
  hybas <-terra::vect(paste0(wd, "/hydroatlas/hybas_eu_lev08_v1c.shp"),
                      extent=c(-10.5,32.702, 34.856, 71.31))
  
}

## Prediction on all river reaches across Europe
current_EC_all <- predict(model_all, data = EU_salt %>%
                            dplyr::select(any_of(names(var_imp_all)))) #prediction on all HYRIV reaches
# Building table with predictions
EC_pred_all <- data.table(HYRIV_ID = EU_salt$HYRIV_ID,
                            pred_EC = as.integer(current_EC_all$predictions))

# merge with shapefile of river network
EC_pred_all_shp <- merge(rivers, EC_pred_all, by = "HYRIV_ID")


#### Export maps ####
#### Recoding environment with 0 as human impact value ####

env_modif <- EU_salt  %>%
  select(any_of(c("HYRIV_ID", names(var_imp_all)))) %>%
  mutate(hft_ix_c09 = 0 ,
         impacted_lc = 0 ,
         pop_ct_csu = 0
         )

##  Prediction with artificially modified environement
baseline_EC <- predict(model_all, data = env_modif %>%
                      dplyr::select(any_of(names(var_imp_all)))) #prediction on all HYRIV reaches
# Tables with predictions
baseline_EC_df <- data.table(HYRIV_ID = env_modif$HYRIV_ID,
                          pred_EC = as.integer(baseline_EC$predictions))

# merge with shapefile of river network
baseline_EC_shp <- merge(rivers, baseline_EC_df, by = "HYRIV_ID")

# salinization = only on catchments with observation
salinization <- EU_salt %>%
  select(HYRIV_ID, EC) %>%
    left_join(baseline_EC_df %>% rename(baseline_EC = pred_EC), by = "HYRIV_ID") %>%
  mutate(salinization = EC - baseline_EC) %>%
  left_join(hyriv_hybas, by = "HYRIV_ID")

# averaging at catchment scale
salinization_avg08 <- salinization %>%
  select(HYBAS_L08, salinization, baseline_EC) %>%
  rename(HYBAS_ID = HYBAS_L08) %>%
  group_by(HYBAS_ID) %>%
  summarize(salinization = mean(salinization, na.rm = T),
            baseline_EC = mean(baseline_EC, na.rm = T)) %>%
  ungroup() %>%
  mutate(salinization_pc = salinization/ (baseline_EC+0.1)*100) # salinization percentage

# # % of negative values
# sum(salinization_avg08$salinization < 0)/nrow(salinization_avg08)*100

# Merging with shapefile of catchments
salinization_hybas_shp <- merge(hybas, salinization_avg08, by = "HYBAS_ID")

# export shapefile baseline EC
writeVector(baseline_EC_shp,
            filename = paste0(wd, "/output/EC_pred_ref_setting0s.shp"),
            overwrite=TRUE)

# export shapefile salinization
writeVector(salinization_hybas_shp,
            filename = paste0(wd, "/output/salinization_setting0s.shp"),
            overwrite=TRUE)
