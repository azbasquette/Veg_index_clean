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

df_reordered_leaf<-read_excel("~/Desktop/df_reordered_leaf_final.xlsx")

##################################################################
#================= PLSR =====================================
############ PRIMARY CHEMISTRY ##################################

df_pls <- df_reordered_leaf %>%
  mutate(
    date = as.character(date),
    year = as.factor(year),
    variety = as.factor(variety),
    year_strata = as.factor(year)
  ) %>%
  filter(
    !is.na(year),
    !is.na(variety)
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

# --- Recipe function: spectra + year + variety ---
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

# --- Primary chemistry response variables ---
pls_rec_brix <- make_pls_rec("brix")
pls_rec_ph   <- make_pls_rec("p_h")
pls_rec_ta   <- make_pls_rec("ta_g_l_tartaric_acid")

rec_list_primchem <- list(pls_rec_brix, pls_rec_ph, pls_rec_ta)

pred_list_primchem <- c(
  "brix",
  "p_h",
  "ta_g_l_tartaric_acid"
)

model_label <- "spectra_year_variety"

df_tune  <- data.frame()
df_vip   <- data.frame()
df_train <- data.frame()
df_test  <- data.frame()

# --- Parallel processing ---
plan(multisession, workers = parallel::detectCores())

# --- Sanity check: confirm predictors are wavelengths + year + variety ---
rec_list_primchem[[1]]$term_info %>%
  dplyr::filter(role == "predictor") %>%
  dplyr::pull(variable)

# --- PLSR loop ---
for (i in seq_along(rec_list_primchem)) {
  cat("In progress:", i, "out of", length(rec_list_primchem), "\n")
  
  pls_rec <- rec_list_primchem[[i]]
  target  <- pred_list_primchem[[i]]
  
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

# --- Prepare plotting data ---
primary_chem_plot <- df_test %>%
  mutate(
    param = recode(
      param,
      "brix" = "Brix",
      "p_h" = "pH",
      "ta_g_l_tartaric_acid" = "TA (g/L tartaric acid)"
    )
  ) %>%
  filter(
    param %in% c(
      "Brix",
      "pH",
      "TA (g/L tartaric acid)"
    )
  ) %>%
  mutate(
    param = factor(
      param,
      levels = c(
        "Brix",
        "pH",
        "TA (g/L tartaric acid)"
      )
    )
  )

# --- Calculate R² and RMSE ---
primary_chem_stats <- primary_chem_plot %>%
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
    c("Brix", "pH", "TA (g/L tartaric acid)"),
    levels = c("Brix", "pH", "TA (g/L tartaric acid)")
  ),
  label = c("(A)", "(B)", "(C)")
)

# --- Plot ---
primary_chem_test_plot <- ggplot(
  primary_chem_plot,
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
    data = primary_chem_stats,
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
    title = "Leaf Spectra: Actual vs Predicted Brix, pH, and TA for Testing Dataset",
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

primary_chem_test_plot

##################################################################
#================= PLSR =====================================
############ ANTHOCYANIN ##################################

# --- Recipe function: spectra + year + variety ---
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

# --- Anthocyanin response variables ---
pls_rec_totalanth_mgb <- make_pls_rec("mg_berry")
pls_rec_m_mgb         <- make_pls_rec("m_in_sample_mg_b")
pls_rec_de_mgb        <- make_pls_rec("de_in_sample_mg_b")
pls_rec_cy_mgb        <- make_pls_rec("cy_in_sample_mg_b")
pls_rec_pet_mgb       <- make_pls_rec("pet_in_sample_mg_b")
pls_rec_peo_mgb       <- make_pls_rec("peo_in_sample_mg_b")

rec_list_anth <- list(
  pls_rec_totalanth_mgb,
  pls_rec_m_mgb,
  pls_rec_de_mgb,
  pls_rec_cy_mgb,
  pls_rec_pet_mgb,
  pls_rec_peo_mgb
)

pred_list_anth <- c(
  "mg_berry",
  "m_in_sample_mg_b",
  "de_in_sample_mg_b",
  "cy_in_sample_mg_b",
  "pet_in_sample_mg_b",
  "peo_in_sample_mg_b"
)

model_label <- "spectra_year_variety"

df_tune  <- data.frame()
df_vip   <- data.frame()
df_train <- data.frame()
df_test  <- data.frame()

# --- Parallel processing ---
plan(multisession, workers = parallel::detectCores())

# --- Sanity check: confirm predictors are wavelengths + year + variety ---
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

# --- Prepare plotting data ---
anth_plot <- df_test %>%
  mutate(
    param = recode(
      param,
      "mg_berry"           = "Total anthocyanins",
      "m_in_sample_mg_b"   = "Malvidin",
      "de_in_sample_mg_b"  = "Delphinidin",
      "cy_in_sample_mg_b"  = "Cyanidin",
      "pet_in_sample_mg_b" = "Petunidin",
      "peo_in_sample_mg_b" = "Peonidin"
    )
  ) %>%
  mutate(
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

# --- Calculate R² and RMSE ---
anth_stats <- anth_plot %>%
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
    c(
      "Total anthocyanins",
      "Malvidin",
      "Delphinidin",
      "Cyanidin",
      "Petunidin",
      "Peonidin"
    ),
    levels = c(
      "Total anthocyanins",
      "Malvidin",
      "Delphinidin",
      "Cyanidin",
      "Petunidin",
      "Peonidin"
    )
  ),
  label = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)")
)

# --- Plot ---
anth_test_plot <- ggplot(
  anth_plot,
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
    data = anth_stats,
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
    title = "Leaf Spectra: Actual vs Predicted Anthocyanins for Testing Dataset",
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
  facet_wrap(~param, scales = "free", nrow = 2)

anth_test_plot

##################################################################
#================= PLSR =====================================
############ FLAVONOL ##################################

df_reordered<-read_excel("~/Desktop/df_final_w904_flav_updated") %>% 
  dplyr::filter(!date %in% as.Date(c("2024-09-04")))

df_pls <- df_reordered_leaf %>%
  mutate(
    date = as.character(date),
    year = as.factor(year),
    variety = as.factor(variety),
    year_strata = as.factor(year)
  )

df_pls <- df_pls %>%
  mutate(
    date = as.Date(date),
    year = factor(format(date, "%Y"))
  )

sum(is.na(df_pls$year))

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

# --- Recipe function: spectra + year + variety ---
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

# --- Flavonol response variables ---
pls_rec_myricetin_mgb <- make_pls_rec("myricetin_mg_b")
pls_rec_quercetin_mgb <- make_pls_rec("quercetin_mg_b")
pls_rec_kaempferol_mgb <- make_pls_rec("kaempferol_mg_b")

rec_list_flav <- list(
  pls_rec_myricetin_mgb,
  pls_rec_quercetin_mgb,
  pls_rec_kaempferol_mgb
)

pred_list_flav <- c(
  "myricetin_mg_b",
  "quercetin_mg_b",
  "kaempferol_mg_b"
)

model_label <- "spectra_year_variety"

df_tune  <- data.frame()
df_vip   <- data.frame()
df_train <- data.frame()
df_test  <- data.frame()

# --- Parallel processing ---
plan(multisession, workers = parallel::detectCores())

# --- Sanity check: confirm predictors are wavelengths + year + variety ---
rec_list_flav[[1]]$term_info %>%
  dplyr::filter(role == "predictor") %>%
  dplyr::pull(variable)

#********temporary kaempferol check
#*# Run only kaempferol
rec_list_flav <- list(pls_rec_kaempferol_mgb)
pred_list_flav <- c("kaempferol_mg_b")

df_tune  <- data.frame()
df_vip   <- data.frame()
df_train <- data.frame()
df_test  <- data.frame()
#***************

# --- PLSR loop ---
for (i in seq_along(rec_list_flav)) {
  cat("In progress:", i, "out of", length(rec_list_flav), "\n")
  
  pls_rec <- rec_list_flav[[i]]
  target  <- pred_list_flav[[i]]
  
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

######### PLOTTING ######### 

# --- Prepare plotting data ---
flav_plot <- df_test %>%
  mutate(
    param = recode(
      param,
      "myricetin_mg_b"     = "Myricetin",
      "quercetin_mg_b"     = "Quercetin",
      "kaempferol_mg_b"    = "Kaempferol"
    )
  ) %>%
  mutate(
    param = factor(
      param,
      levels = c(
        "Myricetin",
        "Quercetin",
        "Kaempferol"
      )
    )
  )

# --- Calculate R² and RMSE ---
flav_stats <- flav_plot %>%
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
    c(
      "Myricetin",
      "Quercetin",
      "Kaempferol"
    ),
    levels = c(
      "Myricetin",
      "Quercetin",
      "Kaempferol"
    )
  ),
  label = c("(A)", "(B)", "(C)")
)

# --- Plot ---
flav_plot <- flav_plot %>%
  mutate(
    date = factor(date)
  )

flav_leaf_plot <- ggplot(
  flav_plot,
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
    data = flav_stats,
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
    title = "Leaf Spectra: Actual vs Predicted Flavonol Composition (mg/berry) for Testing Dataset",
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

flav_leaf_plot

###### save to plot
flavonol_data_leaf <- df_test %>%
  mutate(
    param = recode(
      param,
      "myricetin_mg_b"  = "Myricetin",
      "quercetin_mg_b"  = "Quercetin",
      "kaempferol_mg_b" = "Kaempferol"
    ),
    dataset = "Leaf",
    date = as.character(date)
  ) %>%
  filter(param %in% c("Myricetin", "Quercetin", "Kaempferol"))

flavonol_data_berry <- flavonol_data_berry %>%
  mutate(
    date = as.character(date)
  )
##create plot

combined_plot_data <- bind_rows(
  flavonol_data_berry,
  flavonol_data_leaf
) %>%
  mutate(
    dataset = factor(dataset, levels = c("Berry", "Leaf")),
    param = factor(
      param,
      levels = c("Myricetin", "Quercetin", "Kaempferol")
    ),
    date = factor(date, levels = sort(unique(date)))
  )
###### stats

plot_stats <- combined_plot_data %>%
  group_by(dataset, param) %>%
  summarise(
    r_squared = cor(truth, .pred, use = "complete.obs")^2,
    rmse = sqrt(mean((truth - .pred)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    stat_label = paste0(
      "R² = ", round(r_squared, 2),
      "\nRMSE = ", round(rmse, 3)
    ),
    x_pos = Inf,
    y_pos = -Inf,
    hjust = 1.05,
    vjust = -0.5
  )

###### panel labels

label_df <- combined_plot_data %>%
  distinct(dataset, param) %>%
  arrange(dataset, param) %>%
  mutate(
    panel_label = paste0("(", LETTERS[row_number()], ")")
  )

###### combined plot

combined_plot <- ggplot(
  combined_plot_data,
  aes(x = truth, y = .pred, color = date)
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
    data = plot_stats,
    aes(
      x = x_pos,
      y = y_pos,
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
      label = panel_label
    ),
    inherit.aes = FALSE,
    hjust = -0.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(dataset ~ param, scales = "free", nrow = 2) +  labs(
    x = "Observed value",
    y = "Predicted value",
    color = "Date",
    title = "Observed vs Predicted Flavonols from Berry and Leaf Spectra"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  )

combined_plot

##################################################################
#================= PLSR =====================================
############ TANNIN ##################################
df_pls <- df_reordered_leaf %>%
  mutate(
    date = as.character(date),
    year = as.factor(year),
    variety = as.factor(variety),
    year_strata = as.factor(year)
  )

df_pls <- df_pls %>%
  mutate(
    date = as.Date(date),
    year = factor(format(date, "%Y"))
  )

sum(is.na(df_pls$year))

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


# --- Recipe function: spectra + year + variety ---
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

# --- Tannin response variable ---
pls_rec_tannin_mgb <- make_pls_rec("total_tannin_mgb")

rec_list_tannin <- list(
  pls_rec_tannin_mgb
)

pred_list_tannin <- c(
  "total_tannin_mgb"
)

model_label <- "spectra_year_variety"

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

# --- Prepare plotting data ---
df_test_leaf <-df_test

tannin_plot <- df_test_leaf %>%
  mutate(
    param = recode(
      param,
      "total_tannin_mgb" = "Total tannins"
    )
  ) %>%
  mutate(
    param = factor(
      param,
      levels = c("Total tannins")
    )
  )

# --- Calculate R² and RMSE ---
tannin_stats <- tannin_plot %>%
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
    "Total tannins",
    levels = c("Total tannins")
  ),
  label = "(A)"
)

# --- Plot ---
tannin_test_plot <- ggplot(
  tannin_plot,
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
    title = "Leaf Spectra: Actual vs Predicted Total Tannins for Testing Dataset",
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

tannin_test_plot

### save to plot 

tannin_data_leaf <- df_test_leaf %>%
  mutate(
    param = recode(
      param,
      "total_tannin_mgb" = "Total tannins (mg/berry)"
    ),
    dataset = "Leaf",
    date = as.character(date)
  ) %>%
  filter(param == "Total tannins (mg/berry)")

## create plot

combined_plot_data <- bind_rows(
  tannin_data_berry,
  tannin_data_leaf
) %>%
  mutate(
    dataset = factor(dataset, levels = c("Berry", "Leaf")),
    param = factor(param, levels = c("Total tannins (mg/berry)")),
    date = factor(
      date,
      levels = sort(unique(date))
    )
  )
###### stats

plot_stats <- combined_plot_data %>%
  group_by(dataset, param) %>%
  summarise(
    r_squared = cor(truth, .pred, use = "complete.obs")^2,
    rmse = sqrt(mean((truth - .pred)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    stat_label = paste0(
      "R² = ", round(r_squared, 2),
      "\nRMSE = ", round(rmse, 3)
    ),
    x_pos = Inf,
    y_pos = -Inf,
    hjust = 1.05,
    vjust = -0.5
  )

###### panel labels

label_df <- combined_plot_data %>%
  distinct(dataset, param) %>%
  arrange(dataset, param) %>%
  mutate(
    panel_label = paste0("(", LETTERS[row_number()], ")")
  )

###### combined plot

combined_plot <- ggplot(
  combined_plot_data,
  aes(x = truth, y = .pred, color = date)
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
    data = plot_stats,
    aes(
      x = x_pos,
      y = y_pos,
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
      label = panel_label
    ),
    inherit.aes = FALSE,
    hjust = -0.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(dataset ~ param, scales = "free", nrow = 1) +  labs(
    x = "Observed value",
    y = "Predicted value",
    color = "Date",
    title = "Observed vs Predicted Tannins from Berry and Leaf Spectra"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  )

combined_plot
