# ------ Working environment ####
## Loading packages
pacman::p_load(here, # for the working directory
               data.table, dplyr, janitor, # data wrangling
               usdm, corrplot, # variable correlation
               terra, tidyterra, #gis
               ranger, # random forest
               ggplot2) # graphics 

## Working directory
wd <- here::here()
setwd(paste0(wd, "/data"))

# ------ Loading data ####

## Salinity data, stemming from the GlobSalt dataset
EU_salt <- read.csv("salty_temporal_EC.csv") 

## Environment (CHELSA, FutureStreams, LUCAS LUC, RiverATLAS, GIRES) for all rivers
env <- read.csv("environment.csv", sep=",")

## Geology - Total % of the upstream area covered by each geological class
## for all level-12 BasinATLAS catchments
geology <- read.csv("geology_upstream_12catchments_percentage.csv", sep=",")

## Distance to the sea in km for each river reach
to_the_sea  <- read.csv("to_the_sea.csv", header=T, sep=",",
                        stringsAsFactors = FALSE)

## Total area of seawater intrusion in level-12 BasinATLAS catchments
seawater <- read.csv("seawater_intrusion_basin12.csv", header=T, sep=",",
                     stringsAsFactors = FALSE)

## List of coastal level-12 BasinATLAS catchments
coastal <- read.csv("coastal_catchments.csv", header=T, sep=",",
                    stringsAsFactors = FALSE)

## Equivalence between the ID of river reaches (RiverATLAS) and ID of level-12
hyriv_hybas <- read.csv("hyriv_hybas_equivalence.csv",header=T,sep=",")


# ------ Processing data before running models ####
# Joining environmental data ####

env <- env %>%
  # join the geology to the other environmental variables
  left_join(geology  %>% select(-c(HYBAS_ID, UP_AREA, area_tot)),
            by = "HYRIV_ID")  %>%
  select_if(function(.) n_distinct(.) >1) %>% #remove columns with unique values
  janitor::clean_names() %>%  # to have clean column names  
  rename(HYRIV_ID = hyriv_id) %>%   # uppercase for the ID again after using janitor
  # computing categories of natural and impacted land covers from LUCAS LUC
  mutate(natural_lc = lc3_mean + lc4_mean + lc5_mean + lc7_mean + lc8_mean + lc9_mean,
         impacted_lc = lc13_mean + lc14_mean + lc15_mean)%>%
  # joining the distance to the sea
  left_join(to_the_sea %>% select(-outlet), by = "HYRIV_ID") %>%
  # joining the data on seawater intrusion
  left_join(seawater %>% select(-HYBAS_ID), by = "HYRIV_ID")

# clean global environment
rm(list=c("geology", "seawater", "to_the_sea"))
gc()

# Removing correlated variables ####

# The "env" file already encompassed variables with low collinearity
# Checking for other possible collinearities
var_final = c("HYRIV_ID", "HYBAS_ID",
              # conductivity
              "EC",
              # geology
              names(env %>% select(alkali_basalt:volcanic_alkaline_group)),
              # distance to the sea
              "to_the_sea", 
              # seawater intrusion
              "seawater_area", 
              # water intrusion
              "wt_min_mean", "wt_wm_mean", "wt_range_mean", "wt_ztw_mean",
              # discharge
              "q_dq_mean", "q_zfw_mean", "q_si_mean",
              # bioclimatic variables
              "scd_mean","bio4_mean","bio12_mean","bio15_mean","bio9_mean",
              # land cover
              "impacted_lc", "natural_lc",
              # other variables from RiverATLAS
              "pop_ct_csu", "hft_ix_c09", "ari_ix_uav", "gla_pc_use",
              "ria_ha_csu", "riv_tc_usu", "gwt_cm_cav", "slt_pc_cav",
              "snd_pc_uav", "ele_mt_cmn", "aet_mm_cyr", "ero_kh_cav",
              # probability of drying
              "predprob30")

# cleaning NA and infinite values
env <- env %>%
  select(any_of(var_final)) %>%
  select_if(function(.) n_distinct(.) >1) %>% #remove columns with unique values
  filter_all(all_vars(!is.infinite(.))) %>% # remove rows with Inf values
  na.omit() # remove all rows with NA

# Considering the remaining variables...
var_remain <- as.data.frame(env %>% select(-c(HYRIV_ID)))
# ... and running a VIF procedure on these variables with a 0.7 threshold
vif_remain <- usdm::vifcor(var_remain, th = 0.7)

# Selecting variables retained by the VIF procedure
env <- env %>%
  select(-any_of(vif_remain@excluded))

# Joining the environment to the salinity table
EU_salt  <- EU_salt %>%
  select(-c(year, HYRIV_ID_yr)) %>%
  group_by(HYBAS_ID, HYRIV_ID) %>%
  # Averaging conductivity across years when there are several years
  summarize(EC = mean(EC, na.rm = T)) %>%
  ungroup() %>%
  left_join(env, by = "HYRIV_ID") %>%
  na.omit() %>%
  distinct

# Splitting reference and impacted sites ####
thr_hfi = 50         # threshold for the human footprint index
thr_impact_lc = 0.10 # threshold for the impacted land cover
adding_seawater = TRUE

# coastal model
EU_salt_ref_coast <- EU_salt %>% 
  filter(hft_ix_c09 < thr_hfi) %>%
  filter(impacted_lc < thr_impact_lc) %>%
  filter(HYRIV_ID %in% coastal$HYRIV_ID) 

if (adding_seawater) {
  EU_salt_ref_coast <- EU_salt_ref_coast %>%
    #adding other coastal rivers exposed to seawater intrusion     
    bind_rows(EU_salt %>%
    filter(HYRIV_ID %in% coastal$HYRIV_ID) %>%  
    filter(seawater_area > 0))
}

# inland model
EU_salt_ref_inland <- EU_salt %>% 
  filter(hft_ix_c09 < thr_hfi) %>%
  filter(impacted_lc < thr_impact_lc) %>%
  filter(!(HYRIV_ID %in% coastal$HYRIV_ID))  

# ------ Model of salinity prediction and salinization calculation  ####
for (area in c("coastal", "inland")){ # running separately coastal and inland models
  for (random in c(TRUE, FALSE)){ # whether to include a random variable in the RF
  # Random Forest - predicting baseline salinity ####
  # parameters for the RF
  set.seed(17)
  num.trees = 1000
  # mtry = number of variables to split at each node
  mtry = 7 # default = square roots of nb of predictors, would be 11 here
  # min node size = number of individuals at which the RF is not splitting the node anymore
  min.node.size = 5 # 5 = default for regression 
  
  # Adapting the set of reference sites depending on the coastal/ inland model
  if (area == "coast"){
    EU_salt_ref <- EU_salt_ref_coast
  } else {
    EU_salt_ref <- EU_salt_ref_inland
  }
  
  if (random == TRUE){
  #  generating a random variable to use as threshold on variable importance
  EU_salt_ref$random <- runif(dim(EU_salt_ref)[1], min=0, max=100)
  
  # RF model of conductivity on reference sites
  model_ref<- ranger(EU_salt_ref$EC ~ .,
                      data = EU_salt_ref %>%
                       select(alkali_basalt:predprob30, random), # with random factor
                      num.trees = num.trees,
                      mtry= mtry,
                      min.node.size = min.node.size,
                      replace = T,
                      oob.error = T,
                      keep.inbag = T,
                      num.threads = 4,
                      importance ='impurity',
                      probability = FALSE)
  } else {
    model_ref<- ranger(EU_salt_ref$EC ~ .,
                       data = EU_salt_ref %>%
                         select(alkali_basalt:predprob30), # without random factor
                       num.trees = num.trees,
                       mtry= mtry,
                       min.node.size = min.node.size,
                       replace = T,
                       oob.error = T,
                       keep.inbag = T,
                       num.threads = 4,
                       importance ='impurity',
                       probability = FALSE)
  }
  
  var_imp_ref <- ranger::importance(model_ref)[order(importance(model_ref),
                                                     decreasing=TRUE)]
  r.squared_ref <- model_ref$r.squared
  
  ## Plotting variable importance
  # plot(var_imp_ref, type = "n", ylab= "Variable importance",
  #      main = "Variable importance in predicting conductivity of reference sites")
  # text(seq(1, 114, 1), var_imp_ref, labels=names(var_imp_ref), cex= 0.6, pos=3)
  
  if (random == FALSE){
    # Building salinity map ####
    
    if (! ("rivers" %in% ls(envir = .GlobalEnv))) { # checks if the file is not already loaded to save time
      
      # River network
      rivers <-terra::vect(paste0(wd, "/hydroatlas/HydroRIVERS_v10_eu.shp"),
                           extent= c(-10.5,32.702, 34.856, 71.31)) %>% # already cropped to the right extent
        select(HYRIV_ID) # only keeping the ID column
      
      # Catchments
      hybas <-terra::vect(paste0(wd, "/hydroatlas/hybas_eu_lev08_v1c.shp"),
                          extent=c(-10.5,32.702, 34.856, 71.31))
      
    }
    
    ## Predicting salinity levels with the RF model
    predicted_EC_ref <- predict(model_ref, data = env %>%
                                dplyr::select(any_of(names(var_imp_ref))))
    
    # Building table with predictions
    EC_pred_ref <- data.table(HYRIV_ID = env$HYRIV_ID,
                              pred_EC = as.integer(predicted_EC_ref$predictions)) %>%
      # adding the column with the 8th-level HydroBASIN catchment ID
      left_join(hyriv_hybas %>%
                  select(HYRIV_ID, HYBAS_L08), by = "HYRIV_ID")
    
    # Selecting only relevant catchments according to the coastal/ inland model
    if (area == "coast"){
      EC_pred_ref <- EC_pred_ref %>%
        filter(HYBAS_L08 %in% coastal$HYBAS_L08)
    } else {
      EC_pred_ref <- EC_pred_ref %>%
        filter(!(HYBAS_L08 %in% coastal$HYBAS_L08))
    }
    
    # Merging table with shapefile of river network
    EC_pred_ref_shp <- merge(rivers, EC_pred_ref, by = "HYRIV_ID")
    
    # Building salinization map ####
    salinization <- EU_salt %>%
      select(HYRIV_ID, EC) %>%
      left_join(EC_pred_ref %>% distinct(), by = "HYRIV_ID") %>%
      filter(!is.na(pred_EC)) %>%
      # calculating salinization as the difference between
      # observed EC and predicted baseline EC
      mutate(salinization = EC - pred_EC)
    
    # avering salinization at the 8th-level HydroBASIN catchment level
    salinization_avg08 <- salinization %>%
      select(HYBAS_L08, salinization, pred_EC) %>%
      rename(HYBAS_ID = HYBAS_L08) %>%
      group_by(HYBAS_ID) %>%
      summarize(salinization = mean(salinization, na.rm = T),
                pred_EC = mean(pred_EC, na.rm = T)) %>%
      ungroup() %>%
      mutate(salinization_pc = salinization/ (pred_EC+0.1)*100) # salinization percentage
    
    if (area == "coast"){
      salinization_avg08 <- salinization_avg08 %>%
        filter(HYBAS_ID %in% coastal$HYBAS_L08)
    } else {
      salinization_avg08 <- salinization_avg08 %>%
        filter(!(HYBAS_ID %in% coastal$HYBAS_L08))
    }
    
    sum(salinization_avg08$salinization < 0) / nrow(salinization_avg08) * 100
    
    # Merging with shapefile of river network
    salinization_shp <- merge(rivers, salinization, by = "HYRIV_ID")
    # Merging with shapefile of catchments
    salinization_hybas_shp <- merge(hybas, salinization_avg08, by = "HYBAS_ID")
  
    # Export maps and tables ####
    
    # export shapefile baseline EC
    writeVector(EC_pred_ref_shp,
                filename = paste0(wd, "/output/EC_pred_ref_", area, ".shp"),
                overwrite=TRUE)
    
    # export shapefiles salinization
    writeVector(salinization_hybas_shp,
                filename = paste0(wd, "/output/salinization_", area ,".shp"),
                overwrite=TRUE)
    
    # Prepare export of variable importance
    var_imp_ref = c(r.squared_ref, var_imp_ref)
    names(var_imp_ref)[1] = "r.squared"
    
    # Variable importance in predicting baseline EC
    write.csv(var_imp_ref,
              paste0(wd, "/output/var_importance_", area,".csv"), row.names = TRUE)
    # Final salinization table
    write.csv(salinization,
              paste0(wd, "/output/pred_salinization_", area,".csv"), row.names = FALSE)
    
  } else {
    # Export variable importance ####

    # Prepare export of variable importance
    var_imp_ref = c(r.squared_ref, var_imp_ref)
    names(var_imp_ref)[1] = "r.squared"
    
    # Variable importance including random factor
    write.csv(var_imp_ref,
              paste0(wd, "/output/var_importance_WITH_RANDOM_", area,".csv"), row.names = TRUE)

  }
  }
}

# Group coastal and inland salinization results ####

saliniz_inland  <- read.csv(paste0(wd, "/output/pred_salinization_inland.csv"),
                            header=T, sep=",", stringsAsFactors = FALSE)

saliniz_coast <- read.csv(paste0(wd, "/output/pred_salinization_coastal.csv"),
                          header=T, sep=",", stringsAsFactors = FALSE)

saliniz <- as.data.frame(rbind(saliniz_inland, saliniz_coast)) %>%
  select(HYRIV_ID, HYBAS_L08, EC, pred_EC, salinization) %>%
  rename(HYBAS_ID = HYBAS_L08) %>%
  group_by(HYBAS_ID) %>%
  mutate(salinization_pc = salinization/ (pred_EC+0.1)*100)

write.csv(saliniz,
          paste0(wd, "/output/salinization.csv"), row.names = F)

# ------ Summary salinization length across Europe #########
rivers_df <-terra::vect(paste0(wd, "/hydroatlas/HydroRIVERS_v10_eu.shp"),
                     extent= c(-10.5,32.702, 34.856, 71.31)) %>% # cropped to the right extent
  as.data.frame() %>%
  select(HYRIV_ID, LENGTH_KM)

salinization  <- read.csv(paste0(getwd(), "/output/salinization.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE)

salinization_km <- salinization %>%
  left_join(rivers_df, by = "HYRIV_ID")

tot_length <- sum(salinization_km$LENGTH_KM)
length_salinized <- with(salinization_km, sum(LENGTH_KM[salinization > 0]))
pc_salinized <- length_salinized/tot_length * 100

salinization_pc_km_cat <- salinization_km %>%
  mutate(sal_pc_cat = factor(case_when(salinization_pc <= -50~ "-100 - -50",
                                       salinization_pc <= -10~ "-50 - -10",
                                       salinization_pc <= 10 ~ "-10 - 10",
                                       salinization_pc <= 50 ~ "10-50",
                                       salinization_pc <= 100 ~ "50-100",
                                       salinization_pc <= 500 ~ "100-500",
                                       salinization_pc > 500 ~ ">500"),
         levels= c("-100 - -50", "-50 - -10", "-10 - 10", "10-50","50-100","100-500", ">500"))) %>%
  select(sal_pc_cat, LENGTH_KM) %>%
  group_by(sal_pc_cat) %>%
  summarize(nb_reach = length(LENGTH_KM),
            pc_reach = round(nb_reach/nrow(salinization_km)*100,0),
            LENGTH_KM = round(sum(LENGTH_KM),0)) %>%
  mutate(LENGTH_KM_PC = round(LENGTH_KM / tot_length*100, 0))

# Supplementary Material 2a
salinization_km %>%
  ggplot(aes(x=salinization)) + theme_classic() + xlab("Salinization percentage") + ylab("Count") +
  geom_histogram(binwidth =100, color = "black", fill = "#009999") +
  geom_vline(aes(xintercept=median(salinization)), col = "red", linewidth = 1)

# Supplementary Material 2b
pc99 <- quantile(salinization_km$salinization_pc, 0.99)

salinization_km %>%
#   filter(abs(salinization_pc)<=pc99) %>%
  filter(abs(salinization_pc)<=500) %>%
  ggplot(aes(x=salinization_pc)) + theme_classic() + xlab("Salinization percentage") + ylab("Count") +
  geom_histogram(binwidth =10, color = "black", fill = "#009999") +
  # geom_histogram(binwidth =50, color = "black", fill = "#009999") +
  geom_vline(aes(xintercept=median(salinization_pc)), col = "red", linewidth = 1)

