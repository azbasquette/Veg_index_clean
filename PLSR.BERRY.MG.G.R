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

df_reordered <- read_excel("~/Desktop/spectra/DF/df_reordered_berry_final.xlsx") 

##################################################################
#================= PLSR =====================================
############ ANTHOCYANINS ##################################

df_pls <- df_reordered %>%
  mutate(
    date = as.character(date),
    year = as.factor(year),
    variety = as.factor(variety),
    year_strata = as.factor(year)
  )

# --- Data split into train and test ---
set.seed(234)

perm_split <- initial_split(
  df_pls,
  prop = 0.7,
  strata = year_strata
)

perm_train <- training(perm_split)
perm_test  <- testing(perm_split)

# --- V-fold cross validation ---
set.seed(345)

perm_folds <- vfold_cv(
  perm_train,
  strata = year_strata
)

# --- Identify wavelength columns ---
wl_cols <- names(perm_train) %>%
  keep(~ str_detect(.x, "^x?\\d{3,4}$")) %>%
  keep(~ {
    wl_num <- suppressWarnings(as.numeric(str_remove(.x, "^x")))
    !is.na(wl_num) && wl_num >= 350 && wl_num <= 2500
  })

stopifnot(length(wl_cols) > 0)

# --- Recipe function: COMBINED spectra + metadata ---
metadata_cols <- c("year", "variety")

make_pls_rec <- function(outcome_col,
                         data = perm_train,
                         metadata_cols = c("year", "variety")) {
  
  stopifnot(outcome_col %in% names(data))
  
  # predictors INCLUDED in model: wavelengths + metadata
  
  use_as_predictors <- c(wl_cols, metadata_cols)
  
  # everything else becomes ID
  id_cols <- setdiff(
    names(data),
    c(use_as_predictors, outcome_col)
  )
  
  recipe(
    as.formula(paste(outcome_col, "~ .")),
    data = data
  ) %>%
    update_role(
      all_of(id_cols),
      new_role = "id"
    ) %>%
    update_role(
      all_of(outcome_col),
      new_role = "outcome"
    ) %>%
    step_mutate(
      year = as.factor(year),
      variety = as.factor(variety)
    ) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_naomit(all_outcomes())
}

# --- Anthocyanin response variables: combined spectra + metadata ---
pls_rec_totalanth_mgg <- make_pls_rec("mg_g")
pls_rec_m_mgg         <- make_pls_rec("m_in_sample_mg_g")
pls_rec_de_mgg        <- make_pls_rec("de_in_sample_mg_g")
pls_rec_cy_mgg        <- make_pls_rec("cy_in_sample_mg_g")
pls_rec_pet_mgg       <- make_pls_rec("pet_in_sample_mg_g")
pls_rec_peo_mgg       <- make_pls_rec("peo_in_sample_mg_g")

rec_list_anth <- list(
  pls_rec_totalanth_mgg,
  pls_rec_m_mgg,
  pls_rec_de_mgg,
  pls_rec_cy_mgg,
  pls_rec_pet_mgg,
  pls_rec_peo_mgg
)

pred_list_anth <- c(
  "mg_g",
  "m_in_sample_mg_g",
  "de_in_sample_mg_g",
  "cy_in_sample_mg_g",
  "pet_in_sample_mg_g",
  "peo_in_sample_mg_g"
)

# --- Empty output objects ---
df_tune  <- data.frame()
df_vip   <- data.frame()
df_train <- data.frame()
df_test  <- data.frame()

# --- Parallel processing ---
plan(multisession, workers = parallel::detectCores())

# --- Sanity checks ---
class(perm_train$date)

class(rsample::analysis(
  perm_folds$splits[[1]]
)$date)

# Confirm predictors include wavelengths + date/variety
rec_list_anth[[1]]$term_info %>%
  dplyr::filter(role == "predictor") %>%
  dplyr::pull(variable)

# --- PLSR loop ---
for (i in seq_along(rec_list_anth)) {
  cat("In progress:", i, "out of", length(rec_list_anth), "\n")
  
  pls_rec <- rec_list_anth[[i]]
  target  <- pred_list_anth[[i]]
  
  if (is.null(target) || is.na(target) || !nzchar(target)) {
    warning(paste("⚠️ Skipping index", i, "- invalid or missing target name"))
    next
  }
  
  # Model spec
  pls_tuning <- parsnip::pls(num_comp = tune()) %>%
    set_mode("regression") %>%
    set_engine("mixOmics")
  
  # Workflow
  pls_workflow <- workflow() %>%
    add_recipe(pls_rec) %>%
    add_model(pls_tuning)
  
  # Tune model
  comp_grid <- tibble(num_comp = seq(1, 20, 1))
  
  metrics <- yardstick::metric_set(
    yardstick::mae,
    yardstick::rmse,
    yardstick::rsq
  )
  
  tune_results <- tune::tune_grid(
    pls_workflow,
    resamples = perm_folds,
    grid = comp_grid,
    metrics = metrics
  )
  
  # Choose best model
  tuned_best <- tune_results %>%
    select_by_pct_loss(metric = "rmse", limit = 3, 1)
  
  # Store tuning info
  tune <- tune_results %>%
    collect_metrics() %>%
    mutate(
      param = target,
      best_comp = tuned_best$num_comp,
      model_type = "combined"
    )
  
  df_tune <- bind_rows(df_tune, tune)
  
  if (is.na(tuned_best$num_comp) || tuned_best$num_comp < 1) {
    warning(paste("⚠️ Skipping", target, "- no valid component found"))
    next
  }
  
  # Final model & fit
  updated_pls_model <- parsnip::pls(num_comp = tuned_best$num_comp) %>%
    set_mode("regression") %>%
    set_engine("mixOmics")
  
  updated_workflow <- pls_workflow %>%
    update_model(updated_pls_model)
  
  pls_model <- updated_workflow %>%
    fit(data = perm_train)
  
  # Extract weights
  pls_weights <- pls_model %>%
    extract_fit_parsnip() %>%
    tidy()
  
  weights <- pls_weights %>%
    filter(term != "Y", component == tuned_best$num_comp) %>%
    dplyr::select(value) %>%
    dplyr::rename(weight = value)
  
  # ---- VIP section ----
  pred_vars <- pls_rec$term_info %>%
    dplyr::filter(role == "predictor") %>%
    dplyr::pull(variable) %>%
    unique()
  
  pred_vars <- pred_vars[!is.na(pred_vars) & nzchar(pred_vars)]
  
  keep_cols <- intersect(c(target, pred_vars), names(perm_test))
  
  if (!(target %in% keep_cols) || length(setdiff(keep_cols, target)) == 0) {
    warning(paste("⚠️ Skipping", target, "- missing target or no predictors in perm_test"))
    next
  }
  
  perm_test2 <- perm_test %>%
    dplyr::select(all_of(keep_cols))
  
  xnames <- setdiff(names(perm_test2), target)
  
  perm_test2[xnames] <- lapply(
    perm_test2[xnames],
    function(x) suppressWarnings(as.numeric(x))
  )
  
  keep_x <- xnames[
    colSums(!is.na(as.data.frame(perm_test2[xnames]))) > 0
  ]
  
  perm_test2 <- perm_test2 %>%
    dplyr::select(all_of(c(target, keep_x)))
  
  perm_test2 <- perm_test2 %>%
    dplyr::filter(!is.na(.data[[target]]))
  
  na_prop <- sapply(
    perm_test2[setdiff(names(perm_test2), target)],
    function(x) mean(is.na(x))
  )
  
  good_x <- names(na_prop)[na_prop < 0.5]
  
  perm_test2 <- perm_test2 %>%
    dplyr::select(all_of(c(target, good_x)))
  
  if (ncol(perm_test2) < 2 || nrow(perm_test2) < 3) {
    warning(paste("⚠️ Not enough data after NA filtering for", target, "- skipping VIP calc."))
  } else {
    
    xnames2 <- setdiff(names(perm_test2), target)
    
    perm_test2[xnames2] <- lapply(perm_test2[xnames2], function(x) {
      x <- as.numeric(x)
      mu <- mean(x, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      x[is.na(x)] <- mu
      x
    })
    
    keep_var <- sapply(perm_test2[xnames2], function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 2) return(FALSE)
      isTRUE(stats::var(v) > 0)
    })
    
    nzv_x <- xnames2[keep_var]
    
    if (length(nzv_x) == 0 || nrow(perm_test2) < 3) {
      warning(paste("⚠️ Not enough usable predictors/rows for", target, "- skipping VIP calc."))
    } else {
      
      perm_test2_cc <- perm_test2 %>%
        dplyr::select(all_of(c(target, nzv_x)))
      
      n_pred <- length(nzv_x)
      n_obs  <- nrow(perm_test2_cc)
      
      max_allowed <- min(n_pred, n_obs - 1)
      ncomp_fit <- min(tuned_best$num_comp, max_allowed)
      
      if (!is.finite(ncomp_fit) || is.na(ncomp_fit) || ncomp_fit < 1) {
        warning(paste("⚠️ Invalid ncomp after cleaning for", target, "- skipping VIP calc."))
      } else {
        
        formula_text <- paste0("`", target, "` ~ .")
        
        plsr_mod <- pls::plsr(
          as.formula(formula_text),
          ncomp = ncomp_fit,
          data  = perm_test2_cc,
          validation = "LOO"
        )
        
        SS <- c(plsr_mod$Yloadings)^2 * colSums(plsr_mod$scores^2)
        Wnorm2 <- colSums(plsr_mod$loading.weights^2)
        SSW <- sweep(plsr_mod$loading.weights^2, 2, SS / Wnorm2, "*")
        vip_all <- sqrt(nrow(SSW) * rowSums(SSW) / sum(SS))
        
        wvec <- weights$weight
        
        if (length(wvec) != length(vip_all)) {
          wvec <- rep(NA_real_, length(vip_all))
        }
        
        vip <- tibble::tibble(
          term    = rownames(plsr_mod$loading.weights),
          wl      = suppressWarnings(readr::parse_number(term)),
          vip     = vip_all,
          regcoef = as.numeric(coef(plsr_mod, intercept = FALSE, ncomp = ncomp_fit)),
          weight  = wvec,
          ncomp   = ncomp_fit,
          param   = target,
          model_type = "combined"
        )
        
        df_vip <- dplyr::bind_rows(df_vip, vip)
      }
    }
  }
  
  # ---- Training predictions ----
  pls_train <- pls_model %>%
    predict(new_data = perm_train) %>%
    mutate(
      truth = perm_train[[target]],
      block = perm_train$block,
      date = perm_train$date,
      treatment = perm_train$treatment,
      param = target,
      model = "train",
      model_type = "combined"
    )
  
  df_train <- bind_rows(df_train, pls_train)
  
  # ---- Testing predictions ----
  pls_test <- pls_model %>%
    predict(new_data = perm_test) %>%
    mutate(
      truth = perm_test[[target]],
      block = perm_test$block,
      treatment = perm_test$treatment,
      date = perm_test$date,
      param = target,
      model = "test",
      model_type = "combined"
    )
  
  df_test <- bind_rows(df_test, pls_test)
}

# --- Plot testing predictions: anthocyanins combined ---

plot_df <- df_test %>% 
  mutate(
    param = recode(
      param,
      "mg_g" = "Total anthocyanins",
      "m_in_sample_mg_g" = "Malvidin",
      "de_in_sample_mg_g" = "Delphinidin",
      "cy_in_sample_mg_g" = "Cyanidin",
      "pet_in_sample_mg_g" = "Petunidin",
      "peo_in_sample_mg_g" = "Peonidin"
    ),
    param = factor(
      param,
      levels = c(
        "Total anthocyanins",
        "Malvidin",
        "Delphinidin",
        "Cyanidin",
        "Petunidin",
        "Peonidin"
      )
    )
  )

stats_df <- plot_df %>%
  group_by(param) %>%
  summarise(
    rsq = yardstick::rsq_vec(truth = truth, estimate = .pred),
    rmse = yardstick::rmse_vec(truth = truth, estimate = .pred),
    label = paste0(
      "R² = ", round(rsq, 2),
      "\nRMSE = ", round(rmse, 2)
    ),
    .groups = "drop"
  )

label_df <- tibble::tibble(
  param = factor(
    c(
      "Cyanidin",
      "Delphinidin",
      "Malvidin",
      "Peonidin",
      "Petunidin",
      "Total anthocyanins"
    ),
    levels = c(
      "Cyanidin",
      "Delphinidin",
      "Malvidin",
      "Peonidin",
      "Petunidin",
      "Total anthocyanins"
    )
  ),
  panel_label = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)")
)

avptest_anth_combined <- ggplot(plot_df, aes(truth, .pred, color = date)) +
  geom_point() +
  geom_smooth(
    aes(color = NULL),
    method = "lm",
    formula = y ~ x
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_text(
    data = label_df,
    aes(x = -Inf, y = Inf, label = panel_label),
    inherit.aes = FALSE,
    hjust = -0.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    data = stats_df,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.5
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Berry Spectra: Actual vs Predicted Anthocyanins (mg/g) for Testing Dataset",
    x = "Observed value",
    y = "Predicted value",
    color = "Date"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  facet_wrap(~param, scales = "free")

avptest_anth_combined

##################################################################
#================= PLSR =====================================
################ FLAVONOLS ##################################

#df_reordered <- read_excel("~/Desktop/df_final_w904_flav_updated") %>%
#  filter(as.Date(date) != as.Date("2024-09-04"))

df_pls <- df_reordered %>%
  mutate(
    date = as.character(date),
    year = as.factor(year),
    variety = as.factor(variety),
    year_strata = as.factor(year)
  )

# --- Data split into train and test ---
set.seed(234)

perm_split <- initial_split(
  df_pls,
  prop = 0.7,
  strata = year_strata
)

perm_train <- training(perm_split)
perm_test  <- testing(perm_split)

# --- V-fold cross validation ---
set.seed(345)

perm_folds <- vfold_cv(
  perm_train,
  strata = year_strata
)

# --- Identify wavelength columns ---
wl_cols <- names(perm_train) %>%
  keep(~ str_detect(.x, "^x?\\d{3,4}$")) %>%
  keep(~ {
    wl_num <- suppressWarnings(as.numeric(str_remove(.x, "^x")))
    !is.na(wl_num) && wl_num >= 350 && wl_num <= 2500
  })

stopifnot(length(wl_cols) > 0)
# --- Recipe function ---
metadata_cols <- c("year", "variety")

make_pls_rec <- function(outcome_col,
                         data = perm_train,
                         metadata_cols = c("year", "variety")) {
  
  stopifnot(outcome_col %in% names(data))
  
  # predictors INCLUDED in model: wavelengths + metadata
  
  use_as_predictors <- c(wl_cols, metadata_cols)
  
  # everything else becomes ID
  id_cols <- setdiff(
    names(data),
    c(use_as_predictors, outcome_col)
  )
  
  recipe(
    as.formula(paste(outcome_col, "~ .")),
    data = data
  ) %>%
    update_role(
      all_of(id_cols),
      new_role = "id"
    ) %>%
    update_role(
      all_of(outcome_col),
      new_role = "outcome"
    ) %>%
    step_mutate(
      year = as.factor(year),
      variety = as.factor(variety)
    ) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_naomit(all_outcomes())
}

# --- Response variables ---
pls_rec_myr_mgg <- make_pls_rec("myricetin_mg_g")
pls_rec_que_mgg <- make_pls_rec("quercetin_mg_g")
pls_rec_kae_mgg <- make_pls_rec("kaempferol_mg_g")

rec_list_flav <- list(
  pls_rec_myr_mgg,
  pls_rec_que_mgg,
  pls_rec_kae_mgg
)

pred_list_flav <- c(
  "myricetin_mg_g",
  "quercetin_mg_g",
  "kaempferol_mg_g"
)

df_tune <- data.frame()
df_vip <- data.frame()
df_train <- data.frame()
df_test <- data.frame()

plan(multisession, workers = parallel::detectCores())

#Loop starts here=================================================

for (i in 1:length(rec_list_flav)) {
  cat("In progress:", i, "out of", length(rec_list_flav), "\n")
  
  pls_rec <- rec_list_flav[[i]]
  target <- pred_list_flav[[i]]
  
  # ✅ Skip if target is invalid
  if (is.null(target) || is.na(target) || !nzchar(target)) {
    warning(paste("⚠️ Skipping index", i, "- invalid or missing target name"))
    next
  }
  
  # Model spec
  pls_tuning <- parsnip::pls(num_comp = tune()) %>%
    set_mode("regression") %>%
    set_engine("mixOmics")
  
  # Workflow
  pls_workflow <- workflow() %>%
    add_recipe(pls_rec) %>%
    add_model(pls_tuning)
  
  # Tune model
  comp_grid <- tibble(num_comp = seq(1, 20, 1))
  
  
  tuned_pls_results <- pls_workflow %>%
    tune_grid(resamples = perm_folds,
              grid = comp_grid,
              metrics <- yardstick::metric_set(
                yardstick::mae,
                yardstick::rmse,
                yardstick::rsq
              ))
  
  
  # Choose best model
  tuned_best <- tuned_pls_results %>%
    select_by_pct_loss(metric = "rmse", limit = 3, 1)
  
  # Store tuning info
  tune <- tuned_pls_results %>%
    collect_metrics() %>%
    mutate(param = target, best_comp = tuned_best$num_comp)
  df_tune <- bind_rows(df_tune, tune)
  
  # Skip if no usable model
  if (is.na(tuned_best$num_comp) || tuned_best$num_comp < 1) {
    warning(paste("⚠️ Skipping", target, "- no valid component found"))
    next
  }
  
  # Final model & fit
  updated_pls_model <- parsnip::pls(num_comp = tuned_best$num_comp) %>%
    set_mode("regression") %>%
    set_engine("mixOmics")
  
  updated_workflow <- pls_workflow %>%
    update_model(updated_pls_model)
  
  pls_model <- updated_workflow %>%
    fit(data = perm_train)
  
  # Extract weights
  pls_weights <- pls_model %>% extract_fit_parsnip() %>% tidy()
  weights <- pls_weights %>%
    filter(term != "Y", component == tuned_best$num_comp) %>%
    dplyr::select(value) %>%
    dplyr::rename(weight = value)
  
  # Get variable list (predictors and outcome)
  var_list <- pls_rec$term_info %>%
    filter(role %in% c("predictor", "outcome")) %>%
    pull(variable)
  
  # Subset perm_test
  perm_test2 <- dplyr::select(perm_test, all_of(var_list))
  perm_test2 <- perm_test2[, sapply(perm_test2, is.numeric)]
  
  # ✅ Check target exists in perm_test2
  if (!(target %in% names(perm_test2))) {
    warning(paste("⚠️ Target", target, "not in perm_test2 — skipping."))
    next
  }
  
  max_ncomp <- min(ncol(perm_test2) - 1, tuned_best$num_comp)
  
  # PLSR model for VIP (if enough variables)
  # Build test data with target + wavelength predictors only
  pred_vars <- pls_rec$term_info %>%
    dplyr::filter(role == "predictor") %>%
    dplyr::pull(variable)
  
  perm_test2 <- perm_test %>%
    dplyr::select(all_of(c(target, pred_vars)))
  
  # Coerce predictors to numeric (keep target as-is)
  xnames <- setdiff(names(perm_test2), target)
  perm_test2[xnames] <- lapply(perm_test2[xnames], function(x) suppressWarnings(as.numeric(x)))
  
  # Drop predictor columns that are all NA after coercion
  keep_x <- xnames[colSums(!is.na(as.data.frame(perm_test2[xnames]))) > 0]
  perm_test2 <- perm_test2 %>% dplyr::select(all_of(c(target, keep_x)))
  
  # Compute safe component limit
  n_pred <- length(keep_x)
  n_obs  <- nrow(perm_test2)
  max_allowed <- min(n_pred, n_obs - 1)
  
  # Choose components: can't exceed tuned_best or max_allowed
  ncomp_fit <- min(tuned_best$num_comp, max_allowed)
  
  if (is.finite(ncomp_fit) && ncomp_fit >= 1) {
    formula_text <- paste0("`", target, "` ~ .")
    plsr_mod <- pls::plsr(
      as.formula(formula_text),
      ncomp = ncomp_fit,
      data = perm_test2,
      validation = "LOO"
    )
    
    # VIP calculation
    SS      <- c(plsr_mod$Yloadings)^2 * colSums(plsr_mod$scores^2)
    Wnorm2  <- colSums(plsr_mod$loading.weights^2)
    SSW     <- sweep(plsr_mod$loading.weights^2, 2, SS / Wnorm2, "*")
    vip_all <- sqrt(nrow(SSW) * rowSums(SSW) / sum(SS))
    
    vip_terms <- rownames(plsr_mod$loading.weights)
    
    weights_clean <- pls_weights %>%
      filter(term != "Y", component == tuned_best$num_comp) %>%
      dplyr::select(term, value) %>%
      dplyr::rename(weight = value)
    
    vip <- tibble::tibble(
      term    = vip_terms,
      wl      = suppressWarnings(readr::parse_number(vip_terms)),
      vip     = vip_all,
      regcoef = as.numeric(coef(plsr_mod, intercept = FALSE, ncomp = ncomp_fit)),
      ncomp   = ncomp_fit,
      param   = target,
      model_type = model_label
    ) %>%
      left_join(weights_clean, by = "term")
    
    df_vip <- dplyr::bind_rows(df_vip, vip)
  } else {
    warning(paste("⚠️ Not enough data/predictors for", target, "- skipping VIP calc."))
  }
  # Training predictions
  pls_train <- pls_model %>%
    predict(new_data = perm_train) %>%
    mutate(
      truth = perm_train[[target]],
      block = perm_train$block,
      date = perm_train$date,
      treatment = perm_train$treatment,
      param = target,
      model = "train"
    )
  df_train <- bind_rows(df_train, pls_train)
  
  # Testing predictions
  pls_test <- pls_model %>%
    predict(new_data = perm_test) %>%
    mutate(
      truth = perm_test[[target]],
      block = perm_test$block,
      treatment = perm_test$treatment,
      date = perm_test$date,
      param = target,
      model = "test"
    )
  df_test <- bind_rows(df_test, pls_test)
}

# --- Prepare plotting data ---
df_test_berry<-df_test

flavonol_plot <- df_test %>%
  mutate(
    param = recode(
      param,
      "myricetin_mg_g" = "Myricetin (mg/g)",
      "quercetin_mg_g"  = "Quercetin (mg/g)",
      "kaempferol_mg_g" = "Kaempferol (mg/g)"
    )
  ) %>%
  filter(
    param %in% c(
      "Myricetin (mg/g)",
      "Quercetin (mg/g)",
      "Kaempferol (mg/g)"
    )
  ) %>%
  mutate(
    param = factor(
      param,
      levels = c(
        "Myricetin (mg/g)",
        "Quercetin (mg/g)",
        "Kaempferol (mg/g)"
      )
    )
  )

# --- Calculate R² and RMSE ---
flavonol_stats <- flavonol_plot %>%
  group_by(param) %>%
  summarise(
    r_squared = cor(truth, .pred, use = "complete.obs")^2,
    rmse = sqrt(mean((truth - .pred)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    stat_label = paste0(
      "R² = ", round(r_squared, 2),
      "\nRMSE = ", round(rmse, 3)
    )
  )

# --- Panel labels ---
label_df <- tibble::tibble(
  param = factor(
    c("Myricetin (mg/g)", "Quercetin (mg/g)", "Kaempferol (mg/g)"),
    levels = c("Myricetin (mg/g)", "Quercetin (mg/g)", "Kaempferol (mg/g)")
  ),
  label = c("(A)", "(B)", "(C)")
)

# --- Plot ---
flavonol_berry_plot <- ggplot(
  flavonol_plot,
  aes(truth, .pred, color = date)
) +
  geom_point() +
  geom_smooth(
    aes(color = NULL),
    method = "lm",
    formula = y ~ x
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_text(
    data = flavonol_stats,
    aes(
      x = Inf,
      y = -Inf,
      label = stat_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = -0.5,
    size = 3.5
  ) +
  geom_text(
    data = label_df,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Berry Spectra: Actual vs Predicted Myricetin, Quercetin, and Kaempferol (mg/gerry) for Testing Dataset",
    x = "Observed value",
    y = "Predicted value",
    color = "Date"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  facet_wrap(~param, scales = "free", nrow = 1)

flavonol_berry_plot

###### save to plot with leaf:

flavonol_data_berry <- df_test_berry %>%
  mutate(
    param = recode(
      param,
      "myricetin_mg_g"  = "Myricetin",
      "quercetin_mg_g"  = "Quercetin",
      "kaempferol_mg_g" = "Kaempferol"
    ),
    dataset = "Berry"
  ) %>%
  filter(param %in% c("Myricetin", "Quercetin", "Kaempferol"))


##################################################################
#================= PLSR =====================================
################ TANNIN ##################################

df_pls <- df_reordered %>%
  mutate(
    date = as.character(date),
    year = as.factor(year),
    variety = as.factor(variety),
    year_strata = as.factor(year)
  )

# --- Data split into train and test ---
set.seed(234)

perm_split <- initial_split(
  df_pls,
  prop = 0.7,
  strata = year_strata
)

perm_train <- training(perm_split)
perm_test  <- testing(perm_split)

# --- V-fold cross validation ---
set.seed(345)

perm_folds <- vfold_cv(
  perm_train,
  strata = year_strata
)

# --- Identify wavelength columns ---
wl_cols <- names(perm_train) %>%
  keep(~ str_detect(.x, "^x?\\d{3,4}$")) %>%
  keep(~ {
    v <- suppressWarnings(as.numeric(str_remove(.x, "^x")))
    !is.na(v) && v >= 350 && v <= 2500
  })

stopifnot(length(wl_cols) > 0)

# --- Recipe function: combined spectra + year + variety ---
make_pls_rec <- function(outcome_col,
                         data = perm_train,
                         metadata_cols = c("year", "variety")) {
  
  stopifnot(outcome_col %in% names(data))
  
  use_as_predictors <- intersect(
    c(wl_cols, metadata_cols),
    names(data)
  )
  
  id_cols <- setdiff(
    names(data),
    c(use_as_predictors, outcome_col)
  )
  
  recipe(
    as.formula(paste(outcome_col, "~ .")),
    data = data
  ) %>%
    update_role(all_of(id_cols), new_role = "id") %>%
    update_role(all_of(outcome_col), new_role = "outcome") %>%
    step_mutate(
      year = as.factor(year),
      variety = as.factor(variety)
    ) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_naomit(all_outcomes())
}

# --- Response variable ---
pls_rec_tannin_mgg <- make_pls_rec("total_tannin_mgg")

rec_list_tannin <- list(pls_rec_tannin_mgg)
pred_list_tannin <- c("total_tannin_mgg")

model_label <- "combined"

df_tune  <- data.frame()
df_vip   <- data.frame()
df_train <- data.frame()
df_test  <- data.frame()

# --- Parallel processing ---
plan(multisession, workers = parallel::detectCores())

# --- Sanity check: confirm predictors are wavelengths + year + variety ---
rec_list_tannin[[1]]$term_info %>%
  dplyr::filter(role == "predictor") %>%
  dplyr::pull(variable)

# --- PLSR loop ---
for (i in seq_along(rec_list_tannin)) {
  cat("In progress:", i, "out of", length(rec_list_tannin), "\n")
  
  pls_rec <- rec_list_tannin[[i]]
  target  <- pred_list_tannin[[i]]
  
  if (is.null(target) || is.na(target) || !nzchar(target)) {
    warning(paste("⚠️ Skipping index", i, "- invalid or missing target name"))
    next
  }
  
  # Model spec
  pls_tuning <- parsnip::pls(num_comp = tune()) %>%
    set_mode("regression") %>%
    set_engine("mixOmics")
  
  # Workflow
  pls_workflow <- workflow() %>%
    add_recipe(pls_rec) %>%
    add_model(pls_tuning)
  
  # Tune model
  comp_grid <- tibble(num_comp = seq(1, 20, 1))
  
  metrics <- yardstick::metric_set(
    yardstick::mae,
    yardstick::rmse,
    yardstick::rsq
  )
  
  tune_results <- tune::tune_grid(
    pls_workflow,
    resamples = perm_folds,
    grid = comp_grid,
    metrics = metrics
  )
  
  # Choose best model
  tuned_best <- tune_results %>%
    select_by_pct_loss(metric = "rmse", limit = 3, 1)
  
  # Store tuning info
  tune <- tune_results %>%
    collect_metrics() %>%
    mutate(
      param = target,
      best_comp = tuned_best$num_comp,
      model_type = model_label
    )
  
  df_tune <- bind_rows(df_tune, tune)
  
  if (is.na(tuned_best$num_comp) || tuned_best$num_comp < 1) {
    warning(paste("⚠️ Skipping", target, "- no valid component found"))
    next
  }
  
  # Final model & fit
  updated_pls_model <- parsnip::pls(num_comp = tuned_best$num_comp) %>%
    set_mode("regression") %>%
    set_engine("mixOmics")
  
  updated_workflow <- pls_workflow %>%
    update_model(updated_pls_model)
  
  pls_model <- updated_workflow %>%
    fit(data = perm_train)
  
  # Extract weights
  pls_weights <- pls_model %>%
    extract_fit_parsnip() %>%
    tidy()
  
  weights <- pls_weights %>%
    filter(term != "Y", component == tuned_best$num_comp) %>%
    dplyr::select(term, value) %>%
    dplyr::rename(weight = value)
  
  # ---- VIP section ----
  pred_vars <- pls_rec$term_info %>%
    dplyr::filter(role == "predictor") %>%
    dplyr::pull(variable) %>%
    unique()
  
  pred_vars <- pred_vars[!is.na(pred_vars) & nzchar(pred_vars)]
  
  keep_cols <- intersect(c(target, pred_vars), names(perm_test))
  
  if (!(target %in% keep_cols) || length(setdiff(keep_cols, target)) == 0) {
    warning(paste("⚠️ Skipping", target, "- missing target or no predictors in perm_test"))
    next
  }
  
  perm_test2 <- perm_test %>%
    dplyr::select(all_of(keep_cols))
  
  xnames <- setdiff(names(perm_test2), target)
  
  perm_test2[xnames] <- lapply(
    perm_test2[xnames],
    function(x) suppressWarnings(as.numeric(x))
  )
  
  keep_x <- xnames[
    colSums(!is.na(as.data.frame(perm_test2[xnames]))) > 0
  ]
  
  perm_test2 <- perm_test2 %>%
    dplyr::select(all_of(c(target, keep_x))) %>%
    dplyr::filter(!is.na(.data[[target]]))
  
  na_prop <- sapply(
    perm_test2[setdiff(names(perm_test2), target)],
    function(x) mean(is.na(x))
  )
  
  good_x <- names(na_prop)[na_prop < 0.5]
  
  perm_test2 <- perm_test2 %>%
    dplyr::select(all_of(c(target, good_x)))
  
  if (ncol(perm_test2) < 2 || nrow(perm_test2) < 3) {
    warning(paste("⚠️ Not enough data after NA filtering for", target, "- skipping VIP calc."))
  } else {
    
    xnames2 <- setdiff(names(perm_test2), target)
    
    perm_test2[xnames2] <- lapply(perm_test2[xnames2], function(x) {
      x <- as.numeric(x)
      mu <- mean(x, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      x[is.na(x)] <- mu
      x
    })
    
    keep_var <- sapply(perm_test2[xnames2], function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 2) return(FALSE)
      isTRUE(stats::var(v) > 0)
    })
    
    nzv_x <- xnames2[keep_var]
    
    if (length(nzv_x) == 0 || nrow(perm_test2) < 3) {
      warning(paste("⚠️ Not enough usable predictors/rows for", target, "- skipping VIP calc."))
    } else {
      
      perm_test2_cc <- perm_test2 %>%
        dplyr::select(all_of(c(target, nzv_x)))
      
      n_pred <- length(nzv_x)
      n_obs  <- nrow(perm_test2_cc)
      
      max_allowed <- min(n_pred, n_obs - 1)
      ncomp_fit <- min(tuned_best$num_comp, max_allowed)
      
      if (!is.finite(ncomp_fit) || is.na(ncomp_fit) || ncomp_fit < 1) {
        warning(paste("⚠️ Invalid ncomp after cleaning for", target, "- skipping VIP calc."))
      } else {
        
        formula_text <- paste0("`", target, "` ~ .")
        
        plsr_mod <- pls::plsr(
          as.formula(formula_text),
          ncomp = ncomp_fit,
          data  = perm_test2_cc,
          validation = "LOO"
        )
        
        SS <- c(plsr_mod$Yloadings)^2 * colSums(plsr_mod$scores^2)
        Wnorm2 <- colSums(plsr_mod$loading.weights^2)
        SSW <- sweep(plsr_mod$loading.weights^2, 2, SS / Wnorm2, "*")
        vip_all <- sqrt(nrow(SSW) * rowSums(SSW) / sum(SS))
        
        vip_terms <- rownames(plsr_mod$loading.weights)
        
        vip <- tibble::tibble(
          term    = vip_terms,
          wl      = suppressWarnings(readr::parse_number(vip_terms)),
          vip     = vip_all,
          regcoef = as.numeric(coef(plsr_mod, intercept = FALSE, ncomp = ncomp_fit)),
          ncomp   = ncomp_fit,
          param   = target,
          model_type = model_label
        ) %>%
          left_join(weights, by = "term")
        
        df_vip <- dplyr::bind_rows(df_vip, vip)
      }
    }
  }
  
  # ---- Training predictions ----
  pls_train <- pls_model %>%
    predict(new_data = perm_train) %>%
    mutate(
      truth = perm_train[[target]],
      block = perm_train$block,
      date = perm_train$date,
      year = perm_train$year,
      variety = perm_train$variety,
      treatment = perm_train$treatment,
      param = target,
      model = "train",
      model_type = model_label
    )
  
  df_train <- bind_rows(df_train, pls_train)
  
  # ---- Testing predictions ----
  pls_test <- pls_model %>%
    predict(new_data = perm_test) %>%
    mutate(
      truth = perm_test[[target]],
      block = perm_test$block,
      treatment = perm_test$treatment,
      date = perm_test$date,
      year = perm_test$year,
      variety = perm_test$variety,
      param = target,
      model = "test",
      model_type = model_label
    )
  
  df_test <- bind_rows(df_test, pls_test)
}

####prep plot

df_test_berry<-df_test

tannin_plot <- df_test_berry %>%
  filter(param == "total_tannin_mgg") %>%
  mutate(
    param = recode(
      param,
      "total_tannin_mgg" = "Total Tannin (mg/g)"
    )
  )

# --- Calculate R² and RMSE ---
tannin_stats <- tannin_plot %>%
  summarise(
    r_squared = cor(truth, .pred, use = "complete.obs")^2,
    rmse = sqrt(mean((truth - .pred)^2, na.rm = TRUE))
  ) %>%
  mutate(
    stat_label = paste0(
      "R² = ", round(r_squared, 2),
      "\nRMSE = ", round(rmse, 3)
    )
  )

# --- Panel label ---
label_df <- tibble::tibble(
  label = "(A)"
)

# --- Plot ---
tannin_test_plot <- ggplot(tannin_plot, aes(truth, .pred, color = date)) +
  geom_point() +
  geom_smooth(
    aes(color = NULL),
    method = "lm",
    formula = y ~ x
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_text(
    data = tannin_stats,
    aes(
      x = Inf,
      y = -Inf,
      label = stat_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = -0.5,
    size = 3.5
  ) +
  geom_text(
    data = label_df,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Actual vs Predicted Total Tannin (mg/g) for testing dataset",
    x = "Observed value",
    y = "Predicted value",
    color = "Date"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    plot.margin = margin(10, 20, 10, 20)
  )

tannin_test_plot

### save to plot with leaf

tannin_data_berry <- df_test_berry %>%
  mutate(date = as.character(date),
         param = recode(
           param,
           "total_tannin_mgg"  = "Total tannins (mg/gerry)"
         ),
         dataset = "Berry"
  ) %>%
  filter(param %in% c("Total tannins (mg/gerry)"))