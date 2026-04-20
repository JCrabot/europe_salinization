# ------------------ Working environment ####
## Loading packages
pacman::p_load(here, # for the working directory
               dplyr, tidyr, # data wrangling
               lme4, lmerTest,# for glm
               sf, spdep, gstat, 
               spatialreg # for SAR
               ) 

## Working directory
wd <- here::here()
setwd(paste0(wd, "/data"))

# ------------------ Loading data ####

# Shapefile of the catchments of HydroBASINS
basins <- st_read(paste0(wd, "/hydroatlas/hybas_eu_lev08_v1c.shp"))

# Salinization file: output of the first script
saliniz  <- read.csv(paste0(getwd(), "/output/salinization.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE)

# Environmental variables as compiled in the first script
preds  <- read.csv(paste0(getwd(), "/output/preds.csv"),
                   header=T, sep=",", stringsAsFactors = FALSE)

# Small function to rescale variables
rescale <- function(x){(x-min(x))/(max(x)-min(x))}
rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  span <- diff(rng)
  
  # Handle exceptions
  # -> span is 0 (all values identical)
  # -> span is not finite (Inf, -Inf, NaN)
  if (!is.finite(span) || span == 0) {
    return(rep(0, length(x)))
  }
  
  # Rescale: shift by min, then divide by range
  # Result maps min -> 0 and max -> 1
  (x - rng[1]) / span
}


# Joining salinization and environment tables
preds <- preds  %>%
  select(-EC) %>% # because already present in salinization table
  left_join(saliniz %>% select(-HYBAS_ID), by = "HYRIV_ID") %>%
  select(-c(alkali_basalt:volcanic_alkaline_group)) %>% # removing geology
  mutate(across(where(is.numeric), ~ ifelse(. <= -999, NA, .))) %>% # recoding NAs
  drop_na() %>%
  mutate(across(c(to_the_sea:predprob30), ~rescale(.))) %>% # rescaling predictors
  mutate(HYBAS_ID = as.factor(HYBAS_ID))

# ------------------ Averaging values at catchment scale ####### 
salinization_avg08 <- preds %>%
  select(-HYRIV_ID) %>%
  group_by(HYBAS_ID) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  filter(!is.na(salinization))

# ------------------ lm on raw salinization #####
full_lm_mean_mod <-lm(salinization~.,
                 data=salinization_avg08 %>%
                   select(to_the_sea:predprob30, salinization))

summary(full_lm_mean_mod)

# ------------------ lm on percentage salinization #####
full_lm_mean_mod_pc <-lm(salinization_pc ~.,
                    data=salinization_avg08 %>%
                      select(to_the_sea:predprob30, salinization_pc))

summary(full_lm_mean_mod_pc)
AIC(full_lm_mean_mod_pc)

# ------------------ SAR on salinization ####
basins$HYBAS_ID <- as.character(basins$HYBAS_ID)
salinization_avg08$HYBAS_ID <- as.character(salinization_avg08$HYBAS_ID)

basins_sal <- basins %>%
  inner_join(salinization_avg08, by = "HYBAS_ID")%>%
  select(to_the_sea:predprob30, salinization)

# Convert in spatial object
basins_sal_sp <- as(basins_sal, "Spatial")

# Neighboring with k nearest neighbors
coords <- coordinates(basins_sal_sp)
knn <- knearneigh(coords, k = 5)
nb <- knn2nb(knn)
listw <- nb2listw(nb, style = "W")

# If the following Moran test is significant, it means that
# there is remaining spatial autocorrelation
moran_mc_result <- moran.mc(basins_sal$salinization, listw, nsim = 999)
print(moran_mc_result)

### Spatial AutoRegressive (SAR) model

data_model <- basins_sal_sp@data %>%
  select(to_the_sea:predprob30, salinization)

modele_sar <- lagsarlm(salinization~.,
                       data=data_model,
                       listw = listw)
summary(modele_sar)

# ------------------ SAR on salinization percentage ####

basins$HYBAS_ID <- as.character(basins$HYBAS_ID)
salinization_avg08$HYBAS_ID <- as.character(salinization_avg08$HYBAS_ID)

basins_sal <- basins %>%
  inner_join(salinization_avg08, by = "HYBAS_ID")%>%
  select(to_the_sea:predprob30, salinization_pc)

# Convert in spatial object
basins_sal_sp <- as(basins_sal, "Spatial")

# Neighboring by k nearest neighbors
coords <- coordinates(basins_sal_sp)
knn <- knearneigh(coords, k = 5)  
nb <- knn2nb(knn)
listw <- nb2listw(nb, style = "W")

# # If the following Moran test is significant, it means that
# # there is spatial autocorrelation in salinization in the input data
# moran_mc_result <- moran.mc(basins_sal$salinization_pc, listw, nsim = 999)
# print(moran_mc_result)

### Spatial AutoRegressive (SAR) model

data_model <- basins_sal_sp@data %>%
  select(to_the_sea:predprob30, salinization_pc)

modele_sar <- lagsarlm(salinization_pc~.,
                       data=data_model,
                       listw = listw)
summary(modele_sar)
AIC(modele_sar)

impacts_result <- impacts(modele_sar, listw = listw, R = 1000)  # R for bootstrapping
summary(impacts_result, zstats = TRUE)
