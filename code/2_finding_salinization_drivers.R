#### Working environment ####
setwd("C:/Users/user2380/Documents/Julie/articles/salinization_jrc/data")
# setwd("/home/julie/articles/salinization_jrc/data/") 

pacman::p_load(dplyr, tidyr, # data wrangling
               lme4, lmerTest,# for glm, not used in the end
               sf, spdep, gstat, sp, ggplot2,
               parallel, GWmodel,
               spatialreg
               ) 

#### Loading data ####
rescale <- function(x){(x-min(x))/(max(x)-min(x))}

coastal <- read.csv("selec_coastal_catchments_v2.csv", header=T, sep=",", stringsAsFactors = FALSE)

# saliniz  <- read.csv(paste0(getwd(), "/output/salinization_june_v2.csv"),
saliniz  <- read.csv(paste0(getwd(), "/output/salinization_sept_v2.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE) #%>%
  # filter(!(HYRIV_ID %in% coastal$HYRIV_ID))

# mean((saliniz %>% select(HYRIV_ID, HYBAS_ID) %>% distinct %>%group_by(HYBAS_ID) %>% summarize(HYRIV_ID = length(HYRIV_ID)))$HYRIV_ID)
# sd((saliniz %>% select(HYRIV_ID, HYBAS_ID) %>% distinct %>%group_by(HYBAS_ID) %>% summarize(HYRIV_ID = length(HYRIV_ID)))$HYRIV_ID)


preds  <- read.csv(paste0(getwd(), "/output/preds.csv"),
                   header=T, sep=",", stringsAsFactors = FALSE) %>%
  select(-EC) %>%
  left_join(saliniz %>% select(-HYBAS_ID), by = "HYRIV_ID") %>%
  select(-c(alkali_basalt:volcanic_alkaline_group)) %>%
  # rowwise() %>% 
  mutate(across(where(is.numeric), ~ ifelse(. <= -999, NA, .))) %>% 
  drop_na() %>%
  mutate(across(c(to_the_sea:predprob30), ~rescale(.))) %>%
  mutate(HYBAS_ID = as.factor(HYBAS_ID))

# ------------------  using averaged value at catchment scale - KEPT ####### 
salinization_avg08 <- preds %>%
  select(-HYRIV_ID) %>%
  group_by(HYBAS_ID) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  filter(!is.na(salinization))

# sum(salinization_avg08$salinization < 0)/nrow(salinization_avg08)*100
# sum(salinization_avg08$salinization_pc < 0)/nrow(salinization_avg08)*100

# hist((salinization_avg08 %>% filter(salinization_pc<=500))$salinization_pc, breaks= 100)
# nrow(salinization_avg08 %>% filter(salinization_pc>500))/nrow(salinization_avg08)*100
# hist((salinization_avg08 %>% filter(abs(salinization)<=2000))$salinization, breaks= 100)
# nrow(salinization_avg08 %>% filter(abs(salinization)>2000))/nrow(salinization_avg08)*100
# 
# salinization_avg08 %>%
#   filter(abs(salinization)<=2000) %>%
#   ggplot(aes(x=salinization)) + theme_classic() + xlab("Salinization (µS/cm)") + ylab("Count") +
#   geom_histogram(binwidth =50, color = "black", fill = "#009999") + 
#   geom_vline(aes(xintercept=mean(salinization)), col = "red")
# 
# salinization_avg08 %>%
#   # filter(abs(salinization_pc)<=500) %>%
#   ggplot(aes(x=salinization_pc)) + theme_classic() + xlab("Salinization (µS/cm)") + ylab("Count") +
#   geom_histogram(binwidth =100, color = "black", fill = "#009999") + 
#   geom_vline(aes(xintercept=mean(salinization_pc)), col = "red")
#
# boxplot(salinization_avg08$salinization)
#
# length(unique(saliniz$HYBAS_ID))

##### lm on raw salinization #####
full_lm_mean_mod <-lm(salinization~.,
                 data=salinization_avg08 %>% select(to_the_sea:predprob30, salinization))
# data=salinization_avg08 %>% select(alkali_basalt:predprob30, salinization))

summary(full_lm_mean_mod)

##### lm on percentage salinization #####
full_lm_mean_mod_pc <-lm(salinization_pc ~.,
                    # data=salinization_avg08 %>% select(alkali_basalt:predprob30, salinization_pc))
                    data=salinization_avg08 %>% select(to_the_sea:predprob30, salinization_pc))

summary(full_lm_mean_mod_pc)
AIC(full_lm_mean_mod_pc)

# ##### glm on raw salinization #####
# full_glm_mean_mod <-glm(salinization~.,
#                         # data=salinization_avg08 %>% select(alkali_basalt:predprob30, salinization))
#                         data=salinization_avg08 %>% select(to_the_sea:predprob30, salinization))
# 
# summary(full_glm_mean_mod)
# # null_glm_mean_mod <-lm(salinization~1,data=salinization_avg08%>% select(alkali_basalt:predprob30, salinization))
# 
# ##### glm on percentage salinization #####
# full_glm_mean_pc_mod <-stats::glm(salinization_pc~.,
#                            # data=salinization_avg08 %>% select(alkali_basalt:predprob30, salinization_pc))
#                            data=salinization_avg08 %>%
#                              select(to_the_sea:predprob30, salinization_pc) #%>%
#                            #   filter(salinization_pc>0),
#                            # family = gaussian(link="log")
#                            )
# 
# summary(full_glm_mean_pc_mod)
# # library(pscl)
# pR2(full_glm_mean_pc_mod)
# # null_glm_pc_mod <-lm(salinization_pc~1,data=salinization_avg08%>% select(alkali_basalt:predprob30, salinization_pc))
# # anova(full_glm_pc_mod,null_glm_pc_mod)

# #### Model dredging ####
# library(MuMIn)
# fmla_full <- as.formula(paste("salinization_pc~ ", paste(colnames(salinization_avg08[,2:21]), collapse= "+")))
# pfullmod <- lm(fmla_full, data = salinization_avg08)
# options(na.action = "na.fail")
# options(digits=2)
# selectionmodeles <- MuMIn::dredge(pfullmod, rank="AICc", m.lim=c(0,15))
# # getmodels(selectionmodeles)
# selectionmodeles_inf2 <- subset(selectionmodeles, delta<=2)
# selectionmodeles_inf2
# avgmod.95 <- model.avg(selectionmodeles_inf2, cumsum(weight) <= 0.95)
# confint(avgmod.95)
# sum_avgmod <- summary(avgmod.95)
# sum_avgmod
# var_selec <- row.names(sum_avgmod$coefmat.full %>% as.data.frame() %>% filter(`Pr(>|z|)`<=0.05))[-1]
# 
# fmla_selec <- as.formula(paste("salinization_pc~ ", paste(var_selec, collapse= "+")))
# selec_mod <- glm(fmla_selec, family=gaussian, data = salinization_avg08)
# selec_mod_null <- glm(salinization_pc ~ 1, data = salinization_avg08, family=gaussian)
# summary(selec_mod)
# # pseudoRsquared = 1-(logLik(selec_mod)/logLik(selec_mod_null))
# # pseudoRsquared
# 
# library(pscl)
# pR2(selec_mod) #r2ML = Maximum likelihood pseudo r-squared; r2CU =Cragg and Uhler's pseudo r-squared
# 
# #### checking residuals#### 
# library(DHARMa)
# simulationOutput <- simulateResiduals(fittedModel = full_lm_mean_mod_pc)
# # plot(simulationOutput)
# ks=testUniformity(simulationOutput, plot = F)
# quant=testQuantiles(simulationOutput, plot = F)
# disp=testDispersion(simulationOutput,alternative="two.sided",plot=F) # type="PearsonChisq" ou type="DHARMa" # overdispersion si ratio >1 significatif
# zero=testZeroInflation(simulationOutput,plot=F)
# print(c( ks$p.value,quant$p.value, disp$p.value,zero$p.value))

# # ------------------ using values at river reach - NOT KEPT ####### 
# ##### lmm on raw salinization #####
# lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
# Sys.setenv(OPENBLAS_NUM_THREADS=2)
# modele1_complet <-lmer(salinization~.+(1|HYBAS_ID),
#                        # data=preds %>% select(alkali_basalt:predprob30, salinization, HYBAS_ID),
#                        data=preds %>% select(to_the_sea:predprob30, salinization, HYBAS_ID),
#                        REML=FALSE)
# res_full_model <- summary(modele1_complet)
# View(res_full_model$coefficients)
# # modele1_nul <-lmer(salinization~1+(1|HYBAS_ID),data=preds%>% select(alkali_basalt:predprob30, salinization, HYBAS_ID), REML=FALSE)
# # anova(modele1_complet,modele1_nul)
# 
# ##### lmm on percentage of salinization #####
# modele1_complet <-lmer(salinization_pc~.+(1|HYBAS_ID),
#                        data=preds %>% select(alkali_basalt:predprob30, salinization_pc, HYBAS_ID),
#                        REML=FALSE)
# res_full_model <- summary(modele1_complet)
# View(res_full_model$coefficients)
# modele1_nul <-lmer(salinization~1+(1|HYBAS_ID),data=preds%>% select(alkali_basalt:predprob30, salinization_pc, HYBAS_ID), REML=FALSE)
# anova(modele1_complet,modele1_nul)
# 
# ##### lm on raw salinization #####
# full_lm_mod <-lm(salinization~.,
#                  # data=preds %>% select(alkali_basalt:predprob30, salinization))
#                  data=preds %>% select(to_the_sea:predprob30, salinization))
# 
# summary(full_lm_mod)
# 
# ##### lm on percentage salinization #####
# full_lm_mod_pc <-lm(salinization_pc ~.,
#                     # data=preds %>% select(alkali_basalt:predprob30, salinization_pc))
#                     data=preds %>% select(to_the_sea:predprob30, salinization_pc))
# 
# summary(full_lm_mod_pc)
# 
# # ##### glm on raw salinization #####
# # full_glm_mod <-glm(salinization~.,
# #                    # data=preds %>% select(alkali_basalt:predprob30, salinization), family = gaussian)
# #                    data=preds %>% select(to_the_sea:predprob30, salinization), family = gaussian)
# # 
# # summary(full_glm_mod)
# # 
# # # 1 - (full_glm_mod$null.deviance / full_glm_mod$deviance)
# # 
# # ##### glm on percentage salinization #####
# # full_glm_mod_pc <-glm(salinization_pc ~.,
# #                    # data=preds %>% select(alkali_basalt:predprob30, salinization_pc), family = gaussian)
# #                    data=preds %>% select(to_the_sea:predprob30, salinization_pc), family = gaussian)
# # 
# # summary(full_glm_mod_pc)

# ----------------- Spatial structure  on salinization ####
bbox_terra <- c(-10.5,32.702, 34.856, 71.31)
basins <- st_read("C:/Users/user2380/Documents/Julie/macroclim/data/hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev08_v1c_clean.shp")

basins$HYBAS_ID <- as.character(basins$HYBAS_ID)
salinization_avg08$HYBAS_ID <- as.character(salinization_avg08$HYBAS_ID)

basins_sal <- basins %>%
  inner_join(salinization_avg08, by = "HYBAS_ID")%>%
  select(to_the_sea:predprob30, salinization)

# plot(basins_sal["salinization"])

# Convertir en objet Spatial
basins_sal_sp <- as(basins_sal, "Spatial")

# Créer une matrice de voisinage (queen par exemple)
# neighbors <- poly2nb(basins_sal_sp)
# listw <- nb2listw(neighbors, style = "W") #pb with some polygons
# sum(card(neighbors) == 0)
# no_neighbors <- which(card(neighbors) == 0)
# basins_sal_clean <- basins_sal_sp[-no_neighbors, ]
# # plot(basins_sal_clean)
# neighbors <- poly2nb(basins_sal_clean)
# listw <- nb2listw(neighbors, style = "W")

# voisinage par distance
coords <- coordinates(basins_sal_sp)
# dist_neighbors <- dnearneigh(coords, 0, 10000)  # 0 to 10 km
# listw <- nb2listw(dist_neighbors, style = "W")
# nearest neighbors
knn <- knearneigh(coords, k = 5)
nb <- knn2nb(knn)
listw <- nb2listw(nb, style = "W")


# Si significatif → il y a une autocorrélation spatiale
moran_mc_result <- moran.mc(basins_sal$salinization, listw, nsim = 999)
print(moran_mc_result)

#Modèle spatial autorégressif (SAR)
#Inclut l'effet spatial dans la variable dépendante.
library(spatialreg)
data_model <- basins_sal_sp@data %>%
  select(to_the_sea:predprob30, salinization)
modele_sar <- lagsarlm(salinization~.,
                       data=data_model,
                       listw = listw)
summary(modele_sar)

# Modèle spatial à erreurs (SEM)
# Capture l’autocorrélation dans les résidus du modèle.
modele_sem <- errorsarlm(salinization~.,
                         data=data_model,
                         listw = listw)
summary(modele_sem)

#SAR : si la structure spatiale influence directement la variable réponse
#SEM : si l’influence spatiale vient de facteurs non observés

# Comparer les modèles
AIC(full_lm_mean_mod, modele_sar, modele_sem)

#  Régression géographiquement pondérée (GWR)
# Si tu soupçonnes que les relations changent dans l’espace.

library(GWmodel)
bw <- bw.gwr(salinization~., data = basins_sal_sp, kernel = "exponential", adaptive = T, approach ="AIC")
gwr_model <-  gwr.basic(salinization ~ .,
  data = basins_sal_sp,
  bw = bw,
  kernel = "exponential",
  adaptive = TRUE,
  longlat = TRUE,
  cv = T
)

gwr_model
gwr_model$GW.diagnostic
spplot(gwr_model$SDF, "hft_ix_c09")
spplot(gwr_model$SDF, "q_si_mean")
spplot(gwr_model$SDF, "seawater_area")
spplot(gwr_model$SDF, "slt_pc_cav")
# spplot(gwr_model$SDF, "y")
spplot(gwr_model$SDF, "hft_ix_c09_SE")
spplot(gwr_model$SDF, "seawater_area_SE")

dMat <- gw.dist(dp.locat = coords, p = 2, longlat = TRUE)
# ----------------- Spatial structure  on salinization percentage ####
bbox_terra <- c(-10.5,32.702, 34.856, 71.31)
basins <- st_read("C:/Users/user2380/Documents/Julie/macroclim/data/hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev08_v1c_clean.shp")
# basins <- st_read("/home/julie/macroclim/data/hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev08_v1c_clean.shp")

basins$HYBAS_ID <- as.character(basins$HYBAS_ID)
salinization_avg08$HYBAS_ID <- as.character(salinization_avg08$HYBAS_ID)

basins_sal <- basins %>%
  inner_join(salinization_avg08, by = "HYBAS_ID")%>%
  select(to_the_sea:predprob30, salinization_pc)

# plot(basins_sal["salinization_pc"])

# Convertir en objet Spatial
basins_sal_sp <- as(basins_sal, "Spatial")

# Créer une matrice de voisinage (queen par exemple)
# neighbors <- poly2nb(basins_sal_sp)
# listw <- nb2listw(neighbors, style = "W") #pb with some polygons
# sum(card(neighbors) == 0)
# no_neighbors <- which(card(neighbors) == 0)
# basins_sal_clean <- basins_sal_sp[-no_neighbors, ]
# # plot(basins_sal_clean)
# neighbors <- poly2nb(basins_sal_clean)
# listw <- nb2listw(neighbors, style = "W")

## voisinage par distance
coords <- coordinates(basins_sal_sp)
# dist_neighbors <- dnearneigh(coords, 0, 100000)  # 0 to 100 km
# listw <- nb2listw(dist_neighbors, style = "W")
## nearest neighbors
knn <- knearneigh(coords, k = 5)  
nb <- knn2nb(knn)
listw <- nb2listw(nb, style = "W")

# # Si significatif → il y a une autocorrélation spatiale
# moran_mc_result <- moran.mc(basins_sal$salinization_pc, listw, nsim = 999)
# print(moran_mc_result)

#Modèle spatial autorégressif (SAR)
#Inclut l'effet spatial dans la variable dépendante.
data_model <- basins_sal_sp@data %>%
  select(to_the_sea:predprob30, salinization_pc)
modele_sar <- lagsarlm(salinization_pc~.,
                       data=data_model,
                       listw = listw)
summary(modele_sar)
AIC(modele_sar)

impacts_result <- impacts(modele_sar, listw = listw, R = 1000)  # R for bootstrapping
summary(impacts_result, zstats = TRUE)

# fitted_vals <- fitted(modele_sar)
# observed_vals <- data_model$salinization_pc
# pseudo_r2 <- cor(fitted_vals, observed_vals)^2
# print(pseudo_r2)
# 
# # McFadden's pseudo R²
# null_model <- lagsarlm(salinization_pc ~ 1, data = data_model, listw = listw)
# logLik_full <- logLik(sar_model)
# logLik_null <- logLik(null_model)
# pseudo_r2_mcfadden <- 1 - (as.numeric(logLik_full) / as.numeric(logLik_null))
# print(pseudo_r2_mcfadden)

# # Modèle spatial à erreurs (SEM)
# # Capture l’autocorrélation dans les résidus du modèle.
# modele_sem <- errorsarlm(salinization_pc~.,
#                          data=data_model,
#                          listw = listw)
# summary(modele_sem)
# 
# #SAR : si la structure spatiale influence directement la variable réponse
# #SEM : si l’influence spatiale vient de facteurs non observés
# 
# # Comparer les modèles
# AIC(full_lm_mean_mod_pc, modele_sar, modele_sem)

#  Régression géographiquement pondérée (GWR)
# Si tu soupçonnes que les relations changent dans l’espace.

bw <- bw.gwr(salinization_pc~., data = basins_sal_sp, kernel = "exponential", adaptive = T, approach = "AIC")
bw
gwr_model <-  gwr.basic(salinization_pc ~ .,
                        data = basins_sal_sp,
                        bw = bw,
                        kernel = "exponential",
                        adaptive = T,
                        longlat = TRUE,
                        cv = T
)

gwr_model
gwr_model$GW.diagnostic
# spplot(gwr_model$SDF, "hft_ix_c09")
# spplot(gwr_model$SDF, "impacted_lc")
# spplot(gwr_model$SDF, "bio4_mean")
# spplot(gwr_model$SDF, "seawater_area")
# spplot(gwr_model$SDF, "slt_pc_cav")
# # spplot(gwr_model$SDF, "y")
# spplot(gwr_model$SDF, "hft_ix_c09_SE")
# spplot(gwr_model$SDF, "seawater_area_SE")

gwr_boot <- gwr.bootstrap(
  salinization_pc ~ .,
  data = basins_sal_sp,
  # R=999, longlat=T,
  # dMat = dMat,
  # kernel = "gaussian",
  # k.nearneigh = 1000,
  # adaptive = T,
  # longlat = TRUE,
  # R = 999
)
# 
# dMat <- gw.dist(dp.locat = coords, p = 2, longlat = TRUE)
# 
# gwr_pred <- gwr.predict(
#   salinization_pc ~ .,
#   data = basins_sal_sp,
#   bw = bw,
#   kernel = "exponential",
#   adaptive = TRUE,
#   longlat = TRUE,
#   predictdata = basins_sal_sp
# )
# observed <- basins_sal_sp$salinization_pc
# predicted <- gwr_pred$SDF$prediction
# basins_sal_sp$error <- observed - predicted
# # summary(basins_sal_sp$error)
# quantile(basins_sal_sp$error, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)
# quantile(abs(basins_sal_sp$error), probs = c(0.1, 0.3, 0.7, 0.9), na.rm = TRUE)
# sum(abs(basins_sal_sp@data$error) < 50)/nrow(basins_sal_sp@data)* 100
# boxplot(abs(basins_sal_sp@data$error), outline = F) #outline removes the outlier, caution!!
# spplot(basins_sal_sp, "error",
#        main = "Prediction Error (Observed - Predicted)",
#        at = seq(-85, 98, length.out = 100),  # Set color scale limits here
#        col.regions = colorRampPalette(c("blue", "white", "red"))(99))
# hist((basins_sal_sp@data %>% filter (error <500))$error, breaks=100)

data_model <- basins_sal_sp@data %>%
  select(to_the_sea:predprob30, salinization_pc)

# Sys.time()
# bw0 <- vector()
# for(i in 1:21){
#   formula <- as.formula(paste0("salinization_pc ~ ", names(data_model[i])))
#   bw <- bw.gwr(formula, data = basins_sal_sp,  kernel = "exponential", adaptive = TRUE, approach = "AIC")
#   bw0 <- c(bw0,
#          bw)
# }
# Sys.time()

gw.ms <- gwr.multiscale(salinization_pc ~ .,
                        data = basins_sal_sp,
                        adaptive = T,
                        nlower = 4,
                        max.iterations = 100,
                        # criterion="CVR",
                        approach = "AIC",
                        bws0= rep(16, 21),
                        kernel = "exponential",
                        hatmatrix = F
                        )

gw.ms
gw.ms$GW.diagnostic
# band  width by variable
data.frame(VarName = names(gw.ms$SDF)[1:length(gw.ms[[2]]$bws)], 
           MSGWR_bw = gw.ms[[2]]$bws)

cl <- makeCluster(30)
clusterEvalQ(cl, library(GWmodel)) 

gw.res <- gwr.multiscale(salinization_pc ~ .,
                        data = basins_sal_sp,
                        adaptive = T,
                        nlower = 4,
                        max.iterations = 1000,
                        # criterion="CVR",
                        bws0=gw.ms[[2]]$bws, bw.seled=rep(T, 21),# new line
                        approach = "AIC",
                        kernel = "exponential",
                        hatmatrix = F ,
                        parallel.method = "cluster", parallel.arg = cl
)
stopCluster(cl)

rm(list=setdiff(ls(), c("gw.ms", "gw.res")))

# table of MGWR coefficients
t1 = apply(gw.res$SDF |> 
             select(Intercept:PctBlack) |>
             st_drop_geometry(), 
           2, summary)
# MGWR bandwidths
t2 = gw.ms[[2]]$bws
# join together with a row bind
tab <- rbind(t1, t2)
# add name to last row of tab
rownames(tab)[7] <- "BW"
# transpose tab
tab <- t(round(tab, 3))

save.image(file='environment_gwmodel.RData')