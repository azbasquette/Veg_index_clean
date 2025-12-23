##############################Cleaned Script for VIs##############################

#================= Next troubleshooting: accuracy of VI equations ==========

#ORIGINAL CODE: chris brix v ndvi:

#plt_sp_brix_ndvi <- df_test %>% 
#  filter(param == "brix") %>%
#  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
#  geom_point() +
#  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
#  scale_x_continuous("NDVI", limits = c(0.5,1))+
#  scale_y_continuous("Observed Brix", limits = c(10,30))+
#  scale_color_manual(values = col_all)+
#  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
#           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "brix", "rsq"]$rsq, 2),
#                         cut(df_stats_ndvi[df_stats_ndvi$param == "brix", "p"]$p, 
#                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
#                             labels = c("***", "**", "*", ""))))+
#  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

#:::::::::::Load libraries, data frames, make necessary df alterations::::::::::::::::::::::::::::::::::::

library(plsmod)
library(tidyverse)
library(tidymodels)
library(ggpubr)
library(readxl)
library(dplyr)
library(janitor)
library(lubridate)
library(tidyr)
library(ggplot2)
library(purrr)

berry_spec<-read_excel("~/Desktop/spectra/DATA TABLES_Spectra/berry/rmi_berry_spectra_interpolated_all")%>% slice(-1)

berry_spec <- berry_spec %>%
  mutate(block = str_to_upper(block),  
         date = as.Date(date, format="%m%d%y"),
         variety = case_when(  
           str_to_lower(variety) == "pn" ~ "Pinot Noir",
           str_to_lower(variety) == "cs" ~ "Cabernet Sauvignon",
           str_to_lower(variety) == "tc" ~ "Tinta Cao",
           str_to_lower(variety) == "temp" ~ "Tempranillo",
           str_to_lower(variety) == "syrah" ~ "Syrah",
           str_to_lower(variety) == "sangio" ~ "Sangiovese",
           str_to_lower(variety) == "nebb" ~ "Nebbiolo",
           str_to_lower(variety) == "mourv" ~ "Mourvedre",
           str_to_lower(variety) == "monte" ~ "Montepulciano",
           str_to_lower(variety) == "merlot" ~ "Merlot",
           str_to_lower(variety) == "grenache" ~ "Grenache Noir",
           str_to_lower(variety) == "carig" ~ "Carignan",
           str_to_lower(variety) == "agli" ~ "Aglianico",
           str_to_lower(variety) == "barb" ~ "Barbera",
           str_to_lower(variety) == "zin" ~ "Zinfandel",
           str_to_lower(variety) == "sagran" ~ "Sagrantino",
           str_to_lower(variety) == "tn" ~ "Touriga Nacional",
           TRUE ~ variety  
         ))

berry_spec
berry_chem

#Merging dfs ----------------------------------------
berry_spec <- berry_spec %>% janitor::clean_names()
berry_chem <- berry_chem %>% janitor::clean_names()

berry_spec <- berry_spec %>%
  janitor::clean_names() %>%
  mutate(
    date    = as.Date(as.character(date)),
    block   = as.character(block),
    variety = as.factor(variety),
    rg    = (x675) / (x545),
    ari     = ((1/(x550))-(1/(x700))),
    mari     = (((1/(x550))-(1/(x700)))*(x780)),
    maci = (x940)/(x530),
    nari = ((1/(x550))-(1/(x700)))/((1/(x550))+(1/(x700))),
    #ndai = ((1/(x664))-(1/(x516)))/((1/(x664))+(1/(x516)))
    ndai = ((x664)-(x516))/((x664)+(x516))
  )
berry_spec

berry_chem <- berry_chem %>% 
  mutate(date = as.Date(as.character(date)),
         row = as.character(row),
         vine = as.character(vine),
         block = as.character(block),
         variety = as.factor(variety),
         brix = as.numeric(brix),
         p_h = as.numeric (p_h),
         ta_g_l_tartaric_acid= as.numeric(ta_g_l_tartaric_acid))
berry_chem

#average reps
berry_chem_avg <- berry_chem %>%
  group_by(block, date, variety) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

berry_chem_avg
summary(berry_chem_avg$brix)

#get rid of na's 
berry_spec <- berry_spec %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

berry_chem_avg <- berry_chem_avg %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

berry_spec %>% filter(is.na(date) | is.na(variety))

#check df's
berry_spec
berry_chem_avg

#fix dates
berry_chem_avg <- berry_chem_avg %>%
  dplyr::mutate(date = dplyr::case_when(
    date == as.Date("2023-07-26") ~ as.Date("2023-07-22"),
    date == as.Date("2024-08-06") ~ as.Date("2024-08-05"), 
    date == as.Date("2024-08-20") ~ as.Date("2024-08-19"),  
    TRUE ~ date
  ))
berry_chem_avg
berry_chem_avg%>% count(date)

table(berry_chem_avg$date)
table(berry_spec$date)

berry_chem_avg <- berry_chem_avg %>%
  mutate(variety = case_when(
    variety == "Pinot noir" ~ "Pinot Noir",
    variety == "Cab. Sauv." ~ "Cabernet Sauvignon",
    TRUE ~ variety
  ))

berry_spec <- berry_spec %>%
  mutate(variety = case_when(
    variety == "Carignan" ~ "Carignane",
    TRUE ~ variety
  ))

berry_spec <- berry_spec %>%
  rename_with(~ sub("^x", "", .x), starts_with("x"))  # x400 -> 400

berry_chem_avg %>% count(date)

chem_keep <- berry_chem_avg
chem_keep %>% count(date)

df_data <- berry_spec %>%
  left_join(chem_keep, by = c("date", "block", "variety"))
df_data%>% count(date)

df_reordered <- df_data %>%
  relocate(last_col(offset = 26):last_col(), .after = 4)
df_reordered

colnames(df_reordered)
tail(names(df_reordered),10)
df_reordered <- df_reordered %>%
  mutate(date = case_when(
    date == as.Date("2024-07-26") ~ as.Date("2024-07-22"),
    TRUE ~ date))

#================================ Veg Indicies ==============================

#remove infinites?
df_cor <- df_reordered %>%
  filter(is.finite(mg_l), is.finite(m_in_sample_mg_l),
         is.finite(rg), is.finite(ari),
         is.finite(nari), is.finite(mari),
         is.finite(maci), is.finite(ndai))

#fit_ndai <- lm(mg_l ~ ndai, data = df_cor)
#r2_ndai  <- summary(fit_ndai)$r.squared

# RG vs MALV
fit <- lm(m_in_sample_mg_l ~ rg, data = df_cor)
r2   <- summary(fit)$r.squared
rmse <- sqrt(mean(residuals(fit)^2))
stats_label <- paste0(
  "R² = ", round(r2, 2), "\n",
  "RMSE = ", round(rmse, 2)
)
plt_vi_rg_malv <- df_cor %>% 
  ggplot(
    aes(
      x     = rg,
      y     = m_in_sample_mg_l,
      color = variety
    )
  ) +
  geom_point() +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label,
    size = 4
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Red:Green Index") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  labs(
    title = "Relationship between Malvadin and Red:Green Index",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_rg_malv

#RG VS TOTAL
fit_rg <- lm(mg_l ~ rg, data = df_cor)
r2_rg   <- summary(fit_rg)$r.squared
rmse_rg <- sqrt(mean(residuals(fit_rg)^2))
stats_label <- paste0(
  "R² = ", round(r2_rg, 2), "\n",
  "RMSE = ", round(rmse_rg, 2)
)
plt_vi_rg <- df_cor %>% 
  ggplot(
    aes(
      x     = rg,
      y     = mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Red:Green Index") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Red:Green Index",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_rg

######################## 
#ari VS TOTAL
fit_ari <- lm(mg_l ~ ari, data = df_cor)
r2_ari   <- summary(fit_ari)$r.squared
rmse_ari <- sqrt(mean(residuals(fit_ari)^2))
stats_label_ari <- paste0(
  "R² = ", round(r2_ari, 2), "\n",
  "RMSE = ", round(rmse_ari, 2)
)
plt_vi_ari <- df_cor %>%
  ggplot(
    aes(
      x     = ari,
      y     = mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Anthocyanin Reflectance Index") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_ari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Anthocyanin Reflectance Index",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_ari

#ari VS MALV
fit_ari_malv <- lm(m_in_sample_mg_l ~ ari, data = df_cor)
r2_ari_malv   <- summary(fit_ari_malv)$r.squared
rmse_ari_malv <- sqrt(mean(residuals(fit_ari_malv)^2))
stats_label_malv_ari <- paste0(
  "R² = ", round(r2_ari_malv, 2), "\n",
  "RMSE = ", round(rmse_ari_malv, 2)
)

plt_vi_ari_malv <- df_cor %>%
  ggplot(
    aes(
      x     =ari,
      y     =m_in_sample_mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Anthocyanin Reflectance Index") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_ari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Anthocyanin Reflectance Index",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_ari_malv
########################
#mari VS TOTAL
fit_mari <- lm(mg_l ~ mari, data = df_cor)
r2_mari   <- summary(fit_mari)$r.squared
rmse_mari <- sqrt(mean(residuals(fit_mari)^2))
stats_label_mari <- paste0(
  "R² = ", round(r2_mari, 2), "\n",
  "RMSE = ", round(rmse_mari, 2)
)
plt_vi_mari <- df_cor %>%
  ggplot(
    aes(
      x     = mari,
      y     = mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Modified Anthocyanin Reflectance Index (mARI)") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_mari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Modified Anthocyanin Reflectance Index (mARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_mari

#ari VS MALV
fit_mari_malv <- lm(m_in_sample_mg_l ~ mari, data = df_cor)
r2_mari_malv   <- summary(fit_mari_malv)$r.squared
rmse_mari_malv <- sqrt(mean(residuals(fit_mari_malv)^2))
stats_label_malv_mari <- paste0(
  "R² = ", round(r2_mari_malv, 2), "\n",
  "RMSE = ", round(rmse_mari_malv, 2)
)
plt_vi_mari_malv <- df_cor %>%
  ggplot(
    aes(
      x     = mari,
      y     =m_in_sample_mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Modified Anthocyanin Reflectance Index (mARI)") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_mari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Modified Anthocyanin Reflectance Index (mARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_mari_malv
########################
#nari VS TOTAL
fit_nari <- lm(mg_l ~ nari, data = df_cor)
r2_nari   <- summary(fit_nari)$r.squared
rmse_nari <- sqrt(mean(residuals(fit_nari)^2))
stats_label_nari <- paste0(
  "R² = ", round(r2_nari, 2), "\n",
  "RMSE = ", round(rmse_nari, 2)
)
plt_vi_nari <- df_cor %>%
  ggplot(
    aes(
      x     = nari,
      y     = mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Normalized Anthocyanin Reflectance Index (nARI)") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_nari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Normalized Anthocyanin Reflectance Index (nARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_nari

#nari VS MALV
fit_nari_malv <- lm(m_in_sample_mg_l ~ nari, data = df_cor)
r2_nari_malv   <- summary(fit_nari_malv)$r.squared
rmse_nari_malv <- sqrt(mean(residuals(fit_nari_malv)^2))
stats_label_malv_nari <- paste0(
  "R² = ", round(r2_nari_malv, 2), "\n",
  "RMSE = ", round(rmse_nari_malv, 2)
)
plt_vi_nari_malv <- df_cor %>%
  ggplot(
    aes(
      x     = nari,
      y     =m_in_sample_mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Normalized Anthocyanin Reflectance Index (nARI)") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_nari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Normalized Anthocyanin Reflectance Index (nARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_nari_malv
########################
#maci VS TOTAL
fit_maci <- lm(mg_l ~ maci, data = df_cor)
r2_maci  <- summary(fit_maci)$r.squared
rmse_maci <- sqrt(mean(residuals(fit_maci)^2))
stats_label_maci <- paste0(
  "R² = ", round(r2_maci, 2), "\n",
  "RMSE = ", round(rmse_maci, 2)
)
plt_vi_maci <- df_cor %>%
  ggplot(
    aes(
      x     = maci,
      y     = mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Modified Anthocyanin Content Index (mACI)") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_maci,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Modified Anthocyanin Content Index (mACI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_maci

#maci VS MALV
fit_maci_malv <- lm(m_in_sample_mg_l ~ maci, data = df_cor)
r2_maci_malv   <- summary(fit_maci_malv)$r.squared
rmse_maci_malv <- sqrt(mean(residuals(fit_maci_malv)^2))
stats_label_malv_maci <- paste0(
  "R² = ", round(r2_maci_malv, 2), "\n",
  "RMSE = ", round(rmse_maci_malv, 2)
)
plt_vi_maci_malv <- df_cor %>%
  ggplot(
    aes(
      x     =maci,
      y     =m_in_sample_mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Modified Anthocyanin Content Index (mACI)") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_maci,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Modified Anthocyanin Content Index (mACI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_maci_malv
########################
#ndai VS TOTAL
fit_ndai <- lm(mg_l ~ ndai, data = df_cor)
r2_ndai  <- summary(fit_ndai)$r.squared
rmse_ndai <- sqrt(mean(residuals(fit_ndai)^2))
stats_label_ndai <- paste0(
  "R² = ", round(r2_ndai, 2), "\n",
  "RMSE = ", round(rmse_ndai, 2)
)
plt_vi_ndai <- df_cor %>%
  ggplot(
    aes(
      x     = ndai,
      y     = mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Normalized Difference Anthocyanin Index (NDAI)") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_ndai,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Normalized Difference Anthocyanin Index (NDAI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_ndai

#ndai VS MALV
fit_ndai_malv <- lm(m_in_sample_mg_l ~ ndai, data = df_cor)
r2_ndai_malv   <- summary(fit_ndai_malv)$r.squared
rmse_ndai_malv <- sqrt(mean(residuals(fit_ndai_malv)^2))
stats_label_malv_ndai <- paste0(
  "R² = ", round(r2_ndai_malv, 2), "\n",
  "RMSE = ", round(rmse_ndai_malv, 2)
)
plt_vi_ndai_malv <- df_cor %>%
  ggplot(
    aes(
      x     =ndai,
      y     =m_in_sample_mg_l,
      color = variety
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Normalized Difference Anthocyanin Index (NDAI)") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_ndai,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Normalized Difference Anthocyanin Index (NDAI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_ndai_malv

##############################################################################

####============= This section displays all plots on one figure ==============

###### Not sure about how well the plots transfer over here..
###### Looks like the Maci and Mari are swapped in terms of stats, revisit!!!

# list your veg index columns here
vi_cols <- c("ari", "maci", "ndai", "mari", "rg", "nari")  # adjust names
df_vi_long <- df_reordered %>%
  select(date, all_of(vi_cols), mg_l, m_in_sample_mg_l) %>%
  pivot_longer(
    cols      = all_of(vi_cols),
    names_to  = "vi_name",
    values_to = "vi_value"
  ) %>%
  filter(is.finite(vi_value))
df_stats_total <- df_vi_long %>%
  filter(is.finite(mg_l)) %>%
  group_by(vi_name) %>%
  summarise(
    rsq = summary(lm(mg_l ~ vi_value))$r.squared,
    p   = summary(lm(mg_l ~ vi_value))$coefficients[2, "Pr(>|t|)"],
    .groups = "drop"
  ) %>%
  mutate(
    sig = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    ),
    label = paste("R² =", round(rsq, 2), sig)
  )
plt_vi_total <- df_vi_long %>%
  filter(is.finite(mg_l)) %>%
  ggplot(aes(x = vi_value, y = mg_l, color = as.factor(date))) +
  geom_point() +
  geom_smooth(
    aes(color = NULL),
    color  = "black",
    method = "lm",
    formula = y ~ x,
    se     = FALSE
  ) +
  facet_wrap(~ vi_name, scales = "free_x") +
  scale_x_continuous("Vegetation Index") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  scale_color_manual(values = col_all) +   # drop this line if `col_all` doesn't exist
  
  # annotate per-VI R² and stars in each facet
  geom_text(
    data = df_stats_total,
    aes(
      x     = -Inf,
      y     = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.5,
    size  = 4
  ) +
  labs(
    title = "Relationships between Vegetation Indices and Total Anthocyanins",
    color = "Date"
  ) +
  theme_bw() +
  theme(
    aspect.ratio  = 1,
    legend.title  = element_blank()
  )

plt_vi_total
df_stats_malv <- df_vi_long %>%
  filter(is.finite(m_in_sample_mg_l)) %>%
  group_by(vi_name) %>%
  summarise(
    rsq = summary(lm(m_in_sample_mg_l ~ vi_value))$r.squared,
    p   = summary(lm(m_in_sample_mg_l ~ vi_value))$coefficients[2, "Pr(>|t|)"],
    .groups = "drop"
  ) %>%
  mutate(
    sig = cut(
      p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    ),
    label = paste("R² =", round(rsq, 2), sig)
  )
plt_vi_malv <- df_vi_long %>%
  filter(is.finite(m_in_sample_mg_l)) %>%
  ggplot(aes(x = vi_value, y = m_in_sample_mg_l, color = as.factor(date))) +
  geom_point() +
  geom_smooth(
    aes(color = NULL),
    color  = "black",
    method = "lm",
    formula = y ~ x,
    se     = FALSE
  ) +
  facet_wrap(~ vi_name, scales = "free_x") +
  scale_x_continuous("Vegetation Index") +
  scale_y_continuous("Malvidin (mg/L)") +
  scale_color_manual(values = col_all) +   # optional
  
  geom_text(
    data = df_stats_malv,
    aes(
      x     = -Inf,
      y     = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.5,
    size  = 4
  ) +
  labs(
    title = "Relationships between Vegetation Indices and Malvidin",
    color = "Date"
  ) +
  theme_bw() +
  theme(
    aspect.ratio  = 1,
    legend.title  = element_blank()
  )

plt_vi_malv

#############################################################################
####=======================VI Plots segmented by year

#noticed chris' plots of ndvi show only one year, likely ours will improve if not going across years. 
#can maybe do a single year and a double year one and compare to his saying this is effictive within a year but not effective across years
#do we do the same thing with the hsi model as well? would be cool to be able to say hsi holds upcross years whereas vi's most effictive within one season, which doesnt allow for realtime projection capability.

###########trying to segment by year
vi_cols <- c("ari", "maci", "ndai", "mari", "rg", "nari")

df_vi_long <- df_reordered %>%
  select(date, all_of(vi_cols), mg_l, m_in_sample_mg_l) %>%
  mutate(year = year(date)) %>%    # <-- add year here
  pivot_longer(
    cols      = all_of(vi_cols),
    names_to  = "vi_name",
    values_to = "vi_value"
  ) %>%
  filter(is.finite(vi_value))
make_vi_plot <- function(df_long, response_var, y_label, year_val) {
  # df_long: df_vi_long
  # response_var: "mg_l" or "m_in_sample_mg_l"
  # y_label: pretty label for y-axis
  # year_val: 2023 or 2024
  
  df_year <- df_long %>%
    filter(year == year_val)
  
  # compute stats per VI
  df_stats <- df_year %>%
    filter(is.finite(.data[[response_var]])) %>%
    group_by(vi_name) %>%
    summarise(
      rsq = summary(lm(.data[[response_var]] ~ vi_value))$r.squared,
      p   = summary(lm(.data[[response_var]] ~ vi_value))$coefficients[2, "Pr(>|t|)"],
      .groups = "drop"
    ) %>%
    mutate(
      sig = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "")
      ),
      label = paste("R² =", round(rsq, 2), sig)
    )
  
  # build the plot
  p <- df_year %>%
    filter(is.finite(.data[[response_var]])) %>%
    ggplot(aes(x = vi_value, y = .data[[response_var]], color = as.factor(date))) +
    geom_point() +
    geom_smooth(
      aes(color = NULL),
      color  = "black",
      method = "lm",
      formula = y ~ x,
      se     = FALSE
    ) +
    facet_wrap(~ vi_name, scales = "free_x") +
    scale_x_continuous("Vegetation Index") +
    scale_y_continuous(y_label) +
    scale_color_manual(values = col_all) +   # remove if col_all not defined
    
    geom_text(
      data = df_stats,
      aes(
        x     = -Inf,
        y     = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.5,
      size  = 4
    ) +
    labs(
      title = paste0(
        "Relationships between Vegetation Indices and ",
        y_label,
        " (", year_val, ")"
      ),
      color = "Date"
    ) +
    theme_bw() +
    theme(
      aspect.ratio  = 1,
      legend.title  = element_blank()
    )
  
  p
}
plt_vi_total_2023 <- make_vi_plot(
  df_long   = df_vi_long,
  response_var = "mg_l",
  y_label  = "Total Anthocyanins (mg/L)",
  year_val = 2023
)

plt_vi_total_2024 <- make_vi_plot(
  df_long   = df_vi_long,
  response_var = "mg_l",
  y_label  = "Total Anthocyanins (mg/L)",
  year_val = 2024
)
plt_vi_malv_2023 <- make_vi_plot(
  df_long   = df_vi_long,
  response_var = "m_in_sample_mg_l",
  y_label  = "Malvidin (mg/L)",
  year_val = 2023
)

plt_vi_malv_2024 <- make_vi_plot(
  df_long   = df_vi_long,
  response_var = "m_in_sample_mg_l",
  y_label  = "Malvidin (mg/L)",
  year_val = 2024
)
plt_vi_total_2023
plt_vi_total_2024
plt_vi_malv_2023
plt_vi_malv_2024

############################################################
###================= Segment by Chemistry ================================
# ---- MALVADIN v ALL VIs
chem_var <- "m_in_sample_mg_l"

vi_vars <- c(
  "rg",
  "ari",
  "ndai",
  "maci",
  "nari",  
  "mari"
)

# Optional: nice labels for facets
vi_labels <- c(
  rg   = "Red:Green Index (RG)",
  ari  = "Anthocyanin Reflectance Index (ARI)",
  ndai = "Normalized Difference Anthocyanin Index (NDAI)",
  maci = "Modified Anthocyanin Content Index (mACI)",
  nari  = "Normalized Anthocyanin Reflectance Index (nARI)",
  mari  = "Modified Anthocyanin Reflectance Index (mARI)"
)
df_long <- df_reordered %>%
  select(variety, date, all_of(chem_var), all_of(vi_vars)) %>%
  pivot_longer(
    cols      = all_of(vi_vars),
    names_to  = "vi_name",
    values_to = "vi_value"
  ) %>%
  filter(is.finite(.data[[chem_var]]),
         is.finite(vi_value))

# ---- compute R² and RMSE per VI ----
stats_df <- df_long %>%
  group_by(vi_name) %>%
  group_modify(~ {
    fit  <- lm(m_in_sample_mg_l ~ vi_value, data = .x)
    r2   <- summary(fit)$r.squared
    rmse <- sqrt(mean(residuals(fit)^2))
    
    tibble(
      label = sprintf("R² = %.2f\nRMSE = %.2f", r2, rmse)
    )
  })

# ---- combined figure ----
plt_vi_all <- ggplot(
  df_long,
  aes(
    x     = vi_value,
    y     = m_in_sample_mg_l,
    color = variety
  )
) +
  geom_point() +
  # annotation per panel
  geom_text(
    data = stats_df,
    aes(
      x     = -Inf,
      y     =  Inf,
      label = label
    ),
    hjust = -0.1,
    vjust =  1.1,
    inherit.aes = FALSE,
    size = 3.5
  ) +
  facet_wrap(
    ~ vi_name,
    scales = "free_x",
    labeller = as_labeller(vi_labels)
  ) +
  scale_color_discrete("Variety") +
  labs(
    x = "Vegetation Index value",
    y = "Malvadin (mg/L)",
    title = "Relationship between Malvadin (mg/L) and Vegetation Indices"
  ) +
  theme_bw() +
  theme(
    legend.position  = "right",
    legend.direction = "vertical"
  )

plt_vi_all
#############################
# ---- TOTAL ANTH v ALL VIs

chem_var <- "mg_l"  

vi_vars <- c(
  "rg",
  "ari",
  "ndai",
  "maci",
  "nari",   # replace with your actual VI columns
  "mari"
)

# Optional: nice labels for facets
vi_labels <- c(
  rg   = "Red:Green Index (RG)",
  ari  = "Anthocyanin Reflectance Index (ARI)",
  ndai = "Normalized Difference Anthocyanin Index (NDAI)",
  maci = "Modified Anthocyanin Content Index (mACI)",
  nari  = "Normalized Anthocyanin Reflectance Index (nARI)",
  mari  = "Modified Anthocyanin Reflectance Index (mARI)"
)
df_long <- df_reordered %>%
  select(variety, date, all_of(chem_var), all_of(vi_vars)) %>%
  pivot_longer(
    cols      = all_of(vi_vars),
    names_to  = "vi_name",
    values_to = "vi_value"
  ) %>%
  filter(is.finite(.data[[chem_var]]),
         is.finite(vi_value))

# ---- compute R² and RMSE per VI ----
stats_df <- df_long %>%
  group_by(vi_name) %>%
  group_modify(~ {
    fit  <- lm(mg_l ~ vi_value, data = .x)
    r2   <- summary(fit)$r.squared
    rmse <- sqrt(mean(residuals(fit)^2))
    
    tibble(
      label = sprintf("R² = %.2f\nRMSE = %.2f", r2, rmse)
    )
  })

# ---- combined figure ----
plt_vi_all_total <- ggplot(
  df_long,
  aes(
    x     = vi_value,
    y     = mg_l,
    color = variety
  )
) +
  geom_point() +
  # annotation per panel
  geom_text(
    data = stats_df,
    aes(
      x     = -Inf,
      y     =  Inf,
      label = label
    ),
    hjust = -0.1,
    vjust =  1.1,
    inherit.aes = FALSE,
    size = 3.5
  ) +
  facet_wrap(
    ~ vi_name,
    scales = "free_x",
    labeller = as_labeller(vi_labels)
  ) +
  scale_color_discrete("Variety") +
  labs(
    x = "Vegetation Index value",
    y = "Total Anthocyanins (mg/L)",
    title = "Relationship between Total Anthocyanins (mg/L) and Vegetation Indices"
  ) +
  theme_bw() +
  theme(
    legend.position  = "right",
    legend.direction = "vertical"
  )

plt_vi_all_total

###########################################################################################

## #################################### NOTES #############
#malvadin: 2024 all higher rsquared values except nari ndai and mari
#total anthocyanins higher 2024 r squared in ari maci and rg, higher in 2023 for nari ndai and mari == same
#similarities to be explored here!
#now i want to make sure the number of sampling points in 24 isnt causing anything weird

#####============ Two timepoints per year ========================================

df_reordered <- df_reordered %>%
  group_by(block, date, variety) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")


df_vi_long <- df_reordered %>%
  select(date, all_of(vi_cols), mg_l, m_in_sample_mg_l) %>%
  mutate(
    date = as.Date(date),        # ensure date is Date class
    year = year(date)            # extract year
  ) %>%
  pivot_longer(
    cols      = all_of(vi_cols),
    names_to  = "vi_name",
    values_to = "vi_value"
  ) %>%
  filter(is.finite(vi_value))
str(df_vi_long$year)
table(df_vi_long$year)


#year == 2024 & date %in% as.Date(c("2024-08-19", "2024-08-05"))
make_vi_plot_2024_two_dates <- function(df_reordered, response_var, y_label) {
  
  df_year <- df_reordered %>%
    mutate(
      date = as.Date(date),
      year = year(date)
    ) %>%
    filter(
      year == 2024,
      date %in% as.Date(c("2024-08-19", "2024-08-05"))
    )
  
  
  # compute stats per VI
  df_stats <- df_year %>%
    filter(is.finite(.data[[response_var]])) %>%
    group_by(vi_name) %>%
    summarise(
      rsq = summary(lm(.data[[response_var]] ~ vi_value))$r.squared,
      p   = summary(lm(.data[[response_var]] ~ vi_value))$coefficients[2, "Pr(>|t|)"],
      .groups = "drop"
    ) %>%
    mutate(
      sig = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "")
      ),
      label = paste("R² =", round(rsq, 2), sig)
    )
  
  # build the plot
  p <- df_year %>%
    filter(is.finite(.data[[response_var]])) %>%
    ggplot(aes(x = vi_value, y = .data[[response_var]], color = as.factor(date))) +
    geom_point() +
    geom_smooth(
      aes(color = NULL),
      color  = "black",
      method = "lm",
      formula = y ~ x,
      se     = FALSE
    ) +
    facet_wrap(~ vi_name, scales = "free_x") +
    scale_x_continuous("Vegetation Index") +
    scale_y_continuous(y_label) +
    scale_color_manual(values = col_all) +
    geom_text(
      data = df_stats,
      aes(
        x     = -Inf,
        y     = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.5,
      size  = 4
    ) +
    labs(
      title = paste0(
        "Vegetation Index Relationship with ",
        y_label,
        " (2024 only: Aug 5 & Aug 19)"
      ),
      color = "Date"
    ) +
    theme_bw() +
    theme(
      aspect.ratio  = 1,
      legend.title  = element_blank()
    )
  
  p
}

df_reordered %>% count(date)


plt_vi_total_2024_two_dates <- make_vi_plot_2024_two_dates(
  df_reordered      = df_vi_long,
  response_var = "mg_l",
  y_label      = "Total Anthocyanins (mg/L)"
)

plt_vi_total_2024_two_dates
plt_vi_malv_2024_two_dates <- make_vi_plot_2024_two_dates(
  df_reordered      = df_vi_long,
  response_var = "m_in_sample_mg_l",
  y_label      = "Malvidin (mg/L)"
)

plt_vi_malv_2024_two_dates

########################################################################
#=========== remove july dates from dataset ===========================
#make only post veraison
df_reordered_no_july <- df_reordered %>%
  filter(!date %in% as.Date(c("2024-07-22", "2024-07-26")))

#nari VS TOTAL
fit_nari <- lm(total_anthocyanins ~ nari, data = df_reordered_no_july)
r2_nari   <- summary(fit_nari)$r.squared
rmse_nari <- sqrt(mean(residuals(fit_nari)^2))
stats_label_nari <- paste0(
  "R² = ", round(r2_nari, 2), "\n",
  "RMSE = ", round(rmse_nari, 2)
)
df_cor <- df_reordered_no_july %>%
  filter(is.finite(total_anthocyanins), is.finite(nari))

plt_vi_nari <- df_cor %>%
  ggplot(
    aes(
      x     = nari,
      y     = total_anthocyanins,
      color = variety)
    )+
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Normalized Anthocyanin Reflectance Index (nARI)") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_nari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Normalized Anthocyanin Reflectance Index (nARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_nari

#ari VS MALV
fit_nari_malv <- lm(m_in_sample_mg_l ~ nari, data = df_reordered_no_july)
r2_nari_malv   <- summary(fit_nari_malv)$r.squared
rmse_nari_malv <- sqrt(mean(residuals(fit_nari_malv)^2))
stats_label_malv_nari <- paste0(
  "R² = ", round(r2_nari_malv, 2), "\n",
  "RMSE = ", round(rmse_nari_malv, 2)
)
df_cor <- df_reordered_no_july %>%
  filter(is.finite(m_in_sample_mg_l), is.finite(nari))
plt_vi_nari_malv <- df_cor %>%
  ggplot(
    aes(
      x     = m_in_sample_mg_l,
      y     =nari,
      color = variety,
      group = interaction(variety, date)
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Normalized Anthocyanin Reflectance Index (nARI)") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_nari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Normalized Anthocyanin Reflectance Index (nARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_nari_malv
#ari VS TOTAL

fit_ari <- lm(total_anthocyanins ~ ari, data = df_reordered_no_july)
r2_ari   <- summary(fit_ari)$r.squared
rmse_ari <- sqrt(mean(residuals(fit_ari)^2))
stats_label_ari <- paste0(
  "R² = ", round(r2_ari, 2), "\n",
  "RMSE = ", round(rmse_ari, 2)
)
df_cor <- df_reordered_no_july %>%
  filter(is.finite(total_anthocyanins), is.finite(ari))
plt_vi_ari <- df_cor %>%
  ggplot(
    aes(
      x     = ari,
      y     = total_anthocyanins,
      color = variety,
      group = interaction(variety, date)
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Anthocyanin Reflectance Index") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_ari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Anthocyanin Reflectance Index",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_ari

#ari VS MALV

fit_ari_malv <- lm(m_in_sample_mg_l ~ ari, data = df_reordered_no_july)

r2_ari_malv   <- summary(fit_ari_malv)$r.squared
rmse_ari_malv <- sqrt(mean(residuals(fit_ari_malv)^2))

stats_label_malv_ari <- paste0(
  "R² = ", round(r2_ari_malv, 2), "\n",
  "RMSE = ", round(rmse_ari_malv, 2)
)
df_cor <- df_reordered_no_july %>%
  filter(is.finite(m_in_sample_mg_l), is.finite(ari))

plt_vi_ari_malv <- df_cor %>%
  ggplot(
    aes(
      x     = m_in_sample_mg_l,
      y     =ari,
      color = variety,
      group = interaction(variety, date)
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Anthocyanin Reflectance Index") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_ari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Anthocyanin Reflectance Index",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_ari_malv
#mari VS TOTAL

fit_mari <- lm(total_anthocyanins ~ mari, data = df_reordered_no_july)

r2_mari   <- summary(fit_mari)$r.squared
rmse_mari <- sqrt(mean(residuals(fit_mari)^2))

stats_label_mari <- paste0(
  "R² = ", round(r2_mari, 2), "\n",
  "RMSE = ", round(rmse_mari, 2)
)
df_cor <- df_reordered_no_july %>%
  filter(is.finite(total_anthocyanins), is.finite(mari))

plt_vi_mari <- df_cor %>%
  ggplot(
    aes(
      x     = mari,
      y     = total_anthocyanins,
      color = variety,
      group = interaction(variety, date)
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Modified Anthocyanin Reflectance Index (mARI)") +
  scale_y_continuous("Total Anthocyanins (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_mari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Total Anthocyanins and Anthocyanin Reflectance Index (mARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_mari

#ari VS MALV

fit_mari_malv <- lm(m_in_sample_mg_l ~ mari, data = df_reordered_no_july)

r2_mari_malv   <- summary(fit_mari_malv)$r.squared
rmse_mari_malv <- sqrt(mean(residuals(fit_mari_malv)^2))

stats_label_malv_mari <- paste0(
  "R² = ", round(r2_mari_malv, 2), "\n",
  "RMSE = ", round(rmse_mari_malv, 2)
)
df_cor <- df_reordered_no_july %>%
  filter(is.finite(m_in_sample_mg_l), is.finite(mari))

plt_vi_mari_malv <- df_cor %>%
  ggplot(
    aes(
      x     = m_in_sample_mg_l,
      y     =mari,
      color = variety,
      group = interaction(variety, date)
    )
  ) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  scale_x_continuous("Modified Anthocyanin Reflectance Index (mARI)") +
  scale_y_continuous("Malvadin Concentration (mg/L)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    label = stats_label_malv_mari,
    size  = 4
  ) +
  labs(
    title = "Relationship between Malvadin and Modified Anthocyanin Reflectance Index (mARI)",
    color = "Variety"
  ) +
  theme_bw()

plt_vi_mari_malv

####varying effects on prediction, some better some worse. going to instead rely on the year component

#################################################################################
################=============== END ========================#####################
#################################################################################

