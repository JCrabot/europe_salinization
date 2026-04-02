#--------Gibbs plot ####
ratio  <- read.csv(paste0(getwd(),"/raw_data/Ratios_EC_Cl_Na_Ca_time.csv"),
                   header=T, sep=",", stringsAsFactors = FALSE) %>%
  select(HYBAS_ID, Chloride, Calcium, Sodium) %>%
  setnames(c("HYBAS_ID", "Cl", "Ca", "Na")) %>%
  # select(HYBAS_ID, Conductivity, Chloride, Calcium, Sodium) %>%
  # setnames(c("HYBAS_ID", "EC", "Cl", "Ca", "Na")) %>%
  group_by(HYBAS_ID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(Na_Ca_ratio = Na/(Ca+Na))

saliniz  <- read.csv(paste0(getwd(), "/output/salinization_june.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE) 

saliniz <- saliniz %>%
  select(HYBAS_ID, EC, pred_EC) %>%
  group_by(HYBAS_ID) %>%
  summarize(pred_EC = mean(pred_EC, na.rm = T),
            EC = mean(EC, na.rm = T),
            salinization = EC - pred_EC) %>%
  ungroup() %>%
  # select(-EC) %>%
  mutate(salinization_pc = salinization/ (pred_EC+0.1)*100)

ratio_salin <- saliniz %>%
  right_join(ratio, by = "HYBAS_ID") %>%
  filter(!is.na(salinization))

##Plot Gibbs Diagram
# (a) TDS vs. Na+/(Na+ + Ca²⁺)

breaks <- 10^(-1:5)

ratio_salin %>%
  filter(salinization >0) %>%
  ggplot(aes(x = Na_Ca_ratio, y = EC, color = salinization)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "#FFF999", high = "#660000") +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "red") +
  xlim(0,1) + 
  # scale_y_log10(breaks = breaks) +
  scale_y_log10(limits = c(1,1e4)) +
  # annotation_logticks() +
  labs(x = "Na+ / (Na+ + Ca2+)",
       y = "Observed conductivity (µS/cm)") +
  theme_minimal() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))

# # (b) TDS vs. Cl⁻/(Cl⁻ + HCO₃⁻)
# ggplot(data, aes(x = TDS, y = Cl_HCO3_Ratio)) +
#   geom_point(size = 3, color = "darkgreen") +
#   geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "red") +
#   labs(title = "Gibbs Diagram (Cl- / (Cl- + HCO3-))",
#        x = "Total Dissolved Solids (mg/L)",
#        y = "Cl- / (Cl- + HCO3-)") +
#   theme_minimal()

write.csv(ratio_salin,
          paste0(getwd(),"/output/salinization_ions_ratio_june.csv"), row.names = TRUE)

