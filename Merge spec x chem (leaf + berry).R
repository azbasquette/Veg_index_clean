# Load libraries
library(plsmod)
library(tidyverse)
library(tidymodels)
library(ggpubr)
library(readxl)
library(writexl)
library(purrr)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(janitor)
library(doFuture)
library(future)
library(pls)
library(tibble)
library(tidyr)

#############################################
##################### BERRY DF #####################
#############################################

#read in spectral data
berry_spec<-read_excel("~/Desktop/spectra/DATA TABLES_Spectra/berry/spectra_berry_23_24_25")%>% slice(-1)

#read in berry composition data
berry_chem <- read_excel("~/Desktop/spectra/DATA TABLES_Chemistry/RMI chem 23-24-25 final.xlsx") 

normalize_brix_to_long <- function(df) {
  cols <- names(df)
  if (all(c("Variety","Date","Brix", "Block") %in% cols)) {
    return(df)  # already long
  }
  date_cols <- grep("^Date\\d+$", cols, value = TRUE)
  brix_cols <- gsub("^Date","Brix", date_cols)
  if (length(date_cols) > 0 && all(brix_cols %in% cols) && "Variety" %in% cols) {
    return(
      df %>%
        pivot_longer(
          cols = all_of(c(date_cols, brix_cols)),
          names_to = c(".value","idx"),
          names_pattern = "([A-Za-z]+)(\\d+)"
        )
    )
  }
  df
}

berry_chem    <- normalize_brix_to_long(berry_chem)

# clean berry df

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

# Clean and standardize column names in both data frames
berry_spec <- berry_spec %>% janitor::clean_names()
berry_chem <- berry_chem %>% janitor::clean_names()

berry_spec <- berry_spec %>%
  mutate(date = as.Date(as.character(date)),
         block = as.character(block),
         variety = as.character(variety))
berry_spec

berry_chem <- berry_chem %>%
  mutate(date = as.Date(as.character(date)),
         row = as.character(row),
         vine = as.character(vine),
         block = as.character(block),
         variety = as.character(variety),
         brix = as.numeric(brix),
         p_h = as.numeric (p_h),
         ta_g_l_tartaric_acid= as.numeric(ta_g_l_tartaric_acid))
berry_chem

#berry_chem <- berry_chem %>%
#  filter(brix >= 12 & brix <= 25)

#average replicates in berry_chem

berry_chem_avg <- berry_chem %>%
  group_by(block, date, variety) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

#write_xlsx(berry_chem_avg, "~/Desktop/Berry_chem_average.xlsx")

berry_chem_avg
summary(berry_chem_avg$brix)
summary(berry_chem$brix)
#rm(berry_chem_avg)
summary(berry_chem_avg$p_h)
summary(berry_chem_avg$ta_g_l_tartaric_acid)

#================================ NOT REMOVING 9/04 HERE
library(dplyr)

# --- Remove rows missing key identifiers ---
berry_chem_avg <- berry_chem_avg %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

berry_spec %>%
  filter(is.na(date) | is.na(block) | is.na(variety))

berry_chem_avg <- berry_chem_avg %>%
  mutate(
    date = as.Date(date),
    block = str_trim(as.character(block)),
    variety = str_trim(as.character(variety))
  )

berry_spec <- berry_spec %>%
  mutate(
    date = as.Date(date),
    block = str_trim(as.character(block)),
    variety = str_trim(as.character(variety))
  )

# --- Fix chemistry dates ---
berry_chem_avg <- berry_chem_avg %>%
  mutate(
    date = case_when(
      date == as.Date("2023-07-26") ~ as.Date("2023-07-22"),
      date == as.Date("2024-08-06") ~ as.Date("2024-08-05"),
      date == as.Date("2024-08-20") ~ as.Date("2024-08-19"),
      TRUE ~ date
    )
  )

# --- Check dates before merge ---
table(berry_chem_avg$date)
table(berry_spec$date)

# --- Fix variety names in chemistry dataframe ---
berry_chem_avg <- berry_chem_avg %>%
  mutate(
    variety = case_when(
      variety == "Pinot noir" ~ "Pinot Noir",
      variety == "Cab. Sauv." ~ "Cabernet Sauvignon",
      TRUE ~ variety
    )
  )

# --- Fix variety names in spectra dataframe ---
berry_spec <- berry_spec %>%
  mutate(
    variety = case_when(
      variety == "Carignan" ~ "Carignane",
      TRUE ~ variety
    )
  )

# --- Rename wavelength columns: x400 -> 400 ---
berry_spec <- berry_spec %>%
  rename_with(~ sub("^x", "", .x), starts_with("x"))

# --- Keep chemistry fields ---
chem_keep <- berry_chem_avg

# --- Merge spectra + chemistry ---
# full_join keeps chemistry-only dates such as 9/04,
# even when no spectra exist for that date
df_data <- berry_spec %>%
  full_join(
    chem_keep,
    by = c("date", "block", "variety")
  )

# --- Move chemistry columns after the 4th column ---
df_reordered <- df_data %>%
  relocate(last_col(offset = 26):last_col(), .after = 4)

# --- Check that 9/04 data remains ---
df_reordered_updated <- df_reordered %>%
  mutate(
    date = as.character(date),
    year = as.integer(substr(date, 1, 4)),
    total_flavonol_mgb = rowSums(
      across(c(myricetin_mg_b, quercetin_mg_b, kaempferol_mg_b)),
      na.rm = TRUE
    )
  ) %>%
  relocate(year, .after = mean) %>%
  relocate(total_flavonol_mgb, .after = total_tannin_mgb)

df_reordered_904 <- df_reordered_updated %>%
  arrange(as.Date(date))

write_xlsx(df_reordered_904, "~/Desktop/df_reordered_904.xlsx" )

######### ===================== remove 9/04 ===================================
berry_chem_avg <- berry_chem_avg %>%
  dplyr::filter(!date %in% as.Date(c("2024-09-04")))

#get rid of na's 
berry_spec <- berry_spec %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

berry_chem_avg <- berry_chem_avg %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

berry_spec %>% filter(is.na(date) | is.na(block) | is.na(variety))

#fix dates
berry_chem_avg <- berry_chem_avg %>%
  dplyr::mutate(date = dplyr::case_when(
    date == as.Date("2023-07-26") ~ as.Date("2023-07-22"),
    date == as.Date("2024-08-06") ~ as.Date("2024-08-05"), 
    date == as.Date("2024-08-20") ~ as.Date("2024-08-19"),  
    TRUE ~ date
  ))
berry_chem_avg

table(berry_chem_avg$date)
table(berry_spec$date)

#####merge spectra with berry composition-------------------- 

# Fix names 
berry_chem_avg <- berry_chem_avg %>%
  mutate(variety = case_when(
    variety == "Pinot noir"        ~ "Pinot Noir",
    variety == "Cab. Sauv."        ~ "Cabernet Sauvignon",
    TRUE                           ~ variety
  ))

berry_spec <- berry_spec %>%
  mutate(variety = case_when(
    variety == "Carignan"          ~ "Carignane",
    TRUE                           ~ variety
  ))

berry_spec <- berry_spec %>%
  rename_with(~ sub("^x", "", .x), starts_with("x"))  # x400 -> 400

# bring in chem fields (keep everything)
chem_keep <- berry_chem_avg

# merge spectra + chem
df_data <- berry_spec %>%
  left_join(chem_keep, by = c("date", "block", "variety"))

# move the chem columns after the 4th column
df_reordered <- df_data %>%
  relocate(last_col(offset = 26):last_col(), .after = 4)
df_reordered

df_reordered <- df_reordered %>%
  filter(
    is.finite(brix),
    is.finite(p_h),
    is.finite(ta_g_l_tartaric_acid)
  )
print(df_reordered)

########### print df_reordered with a year and total flavonol column

df_reordered_updated <- df_reordered %>%
  mutate(
    date = as.character(date),
    
    # --- Extract year from date ---
    year = case_when(
      grepl("23", date) ~ 2023,
      grepl("24", date) ~ 2024,
      grepl("25", date) ~ 2025,
      TRUE ~ NA_real_
    ),
    
    # --- Sum flavonols (mg/berry) ---
    total_flavonol_mgb = 
      myricetin_mg_b + quercetin_mg_b + kaempferol_mg_b
  )
print(df_reordered_updated)

# --- Write to file (CSV) ---
write_xlsx(df_reordered_updated, "~/Desktop/df_reordered_berry.xlsx")

########################### ########################### ########################### 
########################### LEAF DF ################################
########################### ########################### ########################### 

#read in spectral data
leaf_spec<-read_excel("~/Desktop/spectra/DATA TABLES_Spectra/leaf/spectra_leaf_23_24_25.xlsx")%>% slice(-1)

#read in leaf composition data
leaf_chem <- read_excel("~/Desktop/spectra/DATA TABLES_Chemistry/RMI chem 23-24-25 final.xlsx") 

#clean leaf df
library(dplyr)
library(stringr)
library(lubridate)

leaf_spec <- leaf_spec %>%
  mutate(block = str_to_upper(block),  
         date = as.Date(date, format="%m/%d/%y"),
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

#Merging dfs ----------------------------------------

# Clean and standardize column names in both data frames
leaf_spec <- leaf_spec %>% janitor::clean_names()
leaf_chem <- leaf_chem %>% janitor::clean_names()

leaf_spec <- leaf_spec %>%
  mutate(date = as.Date(as.character(date)),
         block = as.character(block),
         variety = as.character(variety))
leaf_spec

leaf_chem <- leaf_chem %>%
  mutate(date = as.Date(as.character(date)),
         row = as.character(row),
         vine = as.character(vine),
         block = as.character(block),
         variety = as.character(variety),
         brix = as.numeric(brix),
         p_h = as.numeric (p_h),
         ta_g_l_tartaric_acid= as.numeric(ta_g_l_tartaric_acid))
leaf_chem

#average replicates in leaf_chem
leaf_chem_avg <- leaf_chem %>%
  group_by(block, date, variety) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

leaf_chem_avg
summary(leaf_chem_avg$brix)
summary(leaf_chem_avg$p_h)
summary(leaf_chem_avg$ta_g_l_tartaric_acid)
summary(df_reordered$'1700')

#get rid of na's 
leaf_spec <- leaf_spec %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

leaf_chem_avg <- leaf_chem_avg %>%
  filter(!is.na(date) & !is.na(block) & !is.na(variety))

leaf_spec %>% filter(is.na(date) | is.na(block) | is.na(variety))

#fix dates
leaf_chem_avg <- leaf_chem_avg %>%
  dplyr::mutate(date = dplyr::case_when(
    date == as.Date("2023-09-07") ~ as.Date("2023-09-05"),
    date == as.Date("2023-09-22") ~ as.Date("2023-09-18"),
    date == as.Date("2024-08-05") ~ as.Date("2024-08-06"), 
    date == as.Date("2024-08-19") ~ as.Date("2024-08-20"), 
    date == as.Date("2024-07-22") ~ as.Date("2024-07-26"), 
    date == as.Date("2024-07-26") ~ as.Date("2024-07-26"),
    TRUE ~ date
  ))

leaf_spec <- leaf_spec %>%
  dplyr::mutate(date = dplyr::case_when(
    date == as.Date("2025-09-11") ~ as.Date("2025-09-10"),
    date == as.Date("2025-08-22") ~ as.Date("2025-08-20"),
    date == as.Date("2025-07-31") ~ as.Date("2025-07-30"), 
    date == as.Date("2024-08-19") ~ as.Date("2024-08-20"),  
    TRUE ~ date
  ))

#####merge spectra with leaf composition-------------------- 

# Fix names 
leaf_chem_avg <- leaf_chem_avg %>%
  mutate(variety = case_when(
    variety == "Pinot noir"        ~ "Pinot Noir",
    variety == "Cab. Sauv."        ~ "Cabernet Sauvignon",
    TRUE                           ~ variety
  ))

leaf_spec <- leaf_spec %>%
  mutate(variety = case_when(
    variety == "Carignan"          ~ "Carignane",
    TRUE                           ~ variety
  ))

leaf_spec <- leaf_spec %>%
  rename_with(~ sub("^x", "", .x), starts_with("x"))  # x400 -> 400

# bring in chem fields (keep everything)
chem_keep <- leaf_chem_avg


# merge spectra + chem
df_data <- leaf_spec %>%
  left_join(chem_keep, by = c("date", "block", "variety"))
df_data

# move the chem columns after the 4th column
df_reordered_leaf <- df_data %>%
  relocate(last_col(offset = 26):last_col(), .after = 4)
df_reordered_leaf


write_xlsx(df_reordered_leaf, "~/Desktop/df_reordered_leaf.xlsx")

