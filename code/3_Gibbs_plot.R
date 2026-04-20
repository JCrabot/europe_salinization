# ------------------ Working environment ####
## Loading packages
pacman::p_load(here, # for the working directory
               dplyr, tidyr, # data wrangling
               lme4, lmerTest,# for glm, not used in the end
               sf, spdep, gstat, sp, ggplot2
) 

## Working directory
wd <- here::here()
setwd(paste0(wd, "/data"))

# ------------------ Loading data ####
# Table with ion concentrations
ratio  <- read.csv(paste0(getwd(), "/Ratios_EC_Cl_Na_Ca_time.csv"),
                   header=T, sep=",", stringsAsFactors = FALSE)

# Salinization as previously predicted
saliniz  <- read.csv(paste0(getwd(), "/output/salinization.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE) 

# ------------------ Data wrangling ####
ratio <- ratio %>%
  select(HYBAS_ID, Chloride, Calcium, Sodium) %>%
  setNames(c("HYBAS_ID", "Cl", "Ca", "Na")) %>%
  # averaging across time
  group_by(HYBAS_ID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  # calculating ratio
  mutate(Na_Ca_ratio = Na/(Ca+Na))

saliniz <- saliniz %>%
  select(HYBAS_ID, EC, pred_EC) %>%
  # averaging at catchment scale (across river reaches)
  group_by(HYBAS_ID) %>%
  summarize(pred_EC = mean(pred_EC, na.rm = T),
            EC = mean(EC, na.rm = T),
            salinization = EC - pred_EC) %>%
  ungroup() %>%
  mutate(salinization_pc = salinization/ (pred_EC+0.1)*100) # percentage of salinization

# joining information on ions ratio and salinization
ratio_salin <- saliniz %>%
  right_join(ratio, by = "HYBAS_ID") %>%
  filter(!is.na(salinization))

# Write table
write.csv(ratio_salin,
          paste0(getwd(),"/output/salinization_ions_ratio.csv"), row.names = TRUE)


# ------------------ Plot Gibbs Diagram ####

gibbs_plot <- ratio_salin %>%
  filter(salinization >0) %>% # Focusing on sites where baseline conductivity is not underestimated
  ggplot(aes(x = Na_Ca_ratio, y = EC, color = salinization)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "#FFF999", high = "#660000") +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "red") +
  xlim(0,1) + 
  scale_y_log10(limits = c(1,1e4)) +
  labs(x = "Na+ / (Na+ + Ca2+)",
       y = "Observed conductivity (µS/cm)") +
  theme_minimal() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))

ggsave(paste0(getwd(), "/output/gibbs_plot.png"), plot = gibbs_plot)
