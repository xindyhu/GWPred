##################################################################################
#File name: MICE.R
#Author: Cindy Hu
#Date: Jan 2025
#Purpose: MICE imputation for censored concentrations
# R v4.4.2
##################################################################################

source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes
MCLs <- c(10, 5, 60, 300, 4000)
names(MCLs) <- metal.codes

impute_missing_values <- function(metal.code) {
  # Load Data
  master <- readRDS(paste0("Data_Files/", metal.code, "_df_PredictorsSelected.rds")) %>%
    select(-c(censored.conc, ros.conc)) %>%
    # make aquifer a factor variable
    mutate(aquifer = as.factor(aquifer))
  if(!"detect.limit"%in%colnames(master)){
    # bring back limit of detection
    sup_master <- readRDS(paste0("R_Output/", metal.code, "_RawPredictors.rds")) %>%
      st_drop_geometry() %>%
      dplyr::select(location.id, detect.limit)
    master <- master %>%
      left_join(sup_master, by = "location.id")
  }
  
  percent_missing <- data.frame(
    column = names(master),
    percent_missing = sapply(master, function(x)
      mean(is.na(x)) * 100)
  ) %>%
    arrange(desc(percent_missing))
  
  ### Step 1, missingness has meaning, set missing drainage to zero
  master$drainage[is.na(master$drainage)] = 0
  
  ### Step 2, factor variables, create a missing category
  factor_column_names <- names(master)[sapply(master, is.factor)]
  character_column_names <- names(master)[sapply(master, is.character)]
  factor_column_to_fill <- percent_missing %>%
    filter(column %in% factor_column_names |
             column %in% character_column_names) %>%
    filter(percent_missing > 0) %>%
    pull(column)
  
  master <- master %>%
    mutate(across(
     all_of(factor_column_to_fill),
      ~ forcats::fct_na_value_to_level(.x, "missing")
    ))
  
  ### Step 3, numeric variables, impute missing values with KNN
  numeric_column_to_fill <- percent_missing %>%
    filter(!(column %in% factor_column_names) &
             !(column %in% character_column_names)) %>%
    filter(percent_missing > 0) %>%
    pull(column) %>%
    # do not impute conc, well.depth
    setdiff(c('conc', 'well.depth'))
  
  if(length(numeric_column_to_fill)>0){
    master_imp <- master %>%
    rename(lon = long) %>%
    VIM::kNN(
      variable = numeric_column_to_fill,
      k = 5,
      dist_var = c("lon", "lat"),
      impNA = FALSE
    ) %>%
    dplyr::select(-ends_with("_imp"))
  
  percent_missing_imp <- data.frame(
    column = names(master_imp),
    percent_missing = sapply(master_imp, function(x)
      mean(is.na(x)) * 100)
  ) %>%
    arrange(desc(percent_missing))
  } else{
    master_imp <- master
  }
  
  ### Step 4, impute censored conc with MICE
  
  # if censored samples don't have detect.limit, drop the observation
  master_imp <- master_imp %>% filter(!(censored & is.na(detect.limit)))
  master_imp$conc[master_imp$censored] <- NA
  
  if(nrow(master_imp)<nrow(master)){
    print("pause to check sample size in imputation input.")
  }
  
  print(paste("percent censored is ",
              round(sum(master_imp$censored) / nrow(master_imp) * 100, 2)))
  # We run the mice code with 0 iterations
  imp <- mice(master_imp, maxit = 0)
  # extractr predictorMatrix and methods of imputation
  #predM <- as_tibble(imp$predictorMatrix)
  predM <- matrix(0, ncol = ncol(master_imp), nrow = ncol(master_imp),
                  dimnames = list(names(master_imp), names(master_imp)))
  # Impute conc using only a curated subset of predictors
  important_predictors <- get_important_predictors(metal.code)
  missing_vars <- setdiff(important_predictors, names(master_imp))
  if (length(missing_vars) > 0) {
    stop(paste("Important predictors missing from master_imp:", paste(missing_vars, collapse = ", ")))
  }
  predM["conc", important_predictors] <- 1
  # specify variables not to impute
  cols_to_zero <- percent_missing_imp %>%
    # columns with too much missingness
    filter(percent_missing > 10) %>%
    pull(column) %>%
    # these need to be imputed, so exclude them from cols_to_zero
    setdiff(c("conc"))
  # Use mutate(across()) to set them to 0
  predM <- as_tibble(predM) %>%
    mutate(across(all_of(cols_to_zero), ~ 0))
  # If you need to convert back to a matrix
  predM <- as.matrix(predM)
  rownames(predM) <- colnames(master_imp)
  colnames(predM) <- colnames(master_imp)
  meth <- imp$method
  meth[cols_to_zero] <- ""
  meth["conc"] <- "conc_below_limit"  # Assign custom method
  
  mice.impute.conc_below_limit <- function(y, ry, x, ...) {
    # Ensure "detect.limit" exists in x
    if (!"detect.limit" %in% colnames(x)) {
      stop("Error: 'detect.limit' column is missing from predictor matrix.")
    }
    detect_limit <- x[!ry, "detect.limit"]  # Extract detection limits for missing rows
    if (any(is.na(detect_limit))) {
      stop("Missing detection limit for censored rows.")
    }
    # Base values using PMM
    pmm_vals <- do.call(mice.impute.pmm, list(y = y, ry = ry, x = x))
    
    # Generate random values below the detection limit
    imputed_values <- EnvStats::rlnormTrunc(
      length(pmm_vals),
      min = pmin(0.001, detect_limit / 2),
      # Lower bound (assumes non-negative conc)
      max = detect_limit,
      # Upper bound (detection limit)
      meanlog = log(pmm_vals + 0.001),
      sdlog = sd(log(y[ry]), na.rm = TRUE)
    )  # Use observed SD
    
    return(imputed_values)
  }
  exists("mice.impute.conc_below_limit")
  # MICE imputation step
  imp2 <- tryCatch({
    message("Attempting MICE imputation (first try)...")
    
    # First attempt
    mice(
      master_imp,
      maxit = 5,
      predictorMatrix = predM,
      method = meth,
      print = TRUE,
      seed = 123
    )
  }, error = function(e) {
    message("Error encountered: ", e$message)
    message("Handling singularity: Removing near-zero variance variables and retrying MICE...")
    
    
    # if failed due to singularity, try again after removing the vars with near zero variance
    master_imp <- master_imp %>%
      dplyr::select(-(caret::nearZeroVar(master_imp, freqCut = 999 / 1)))
    # We run the mice code with 0 iterations
    imp <- mice(master_imp, maxit = 0)
    # extractr predictorMatrix and methods of imputation
    #predM <- as_tibble(imp$predictorMatrix)
    predM <- matrix(0, ncol = ncol(master_imp), nrow = ncol(master_imp),
                    dimnames = list(names(master_imp), names(master_imp)))
    # Impute conc using only a curated subset of predictors
    missing_vars <- setdiff(important_predictors, names(master_imp))
    if (length(missing_vars) > 0) {
      stop(paste("Important predictors missing from master_imp:", paste(missing_vars, collapse = ", ")))
    }
    predM["conc", important_predictors] <- 1
    # specify variables not to impute
    cols_to_zero <- percent_missing_imp %>%
      # columns with too much missingness
      filter(percent_missing > 10) %>%
      pull(column) %>%
      # these need to be imputed, so exclude them from cols_to_zero
      setdiff(c("conc"))
    # Use mutate(across()) to set them to 0
    predM <- as_tibble(predM) %>%
      mutate(across(all_of(cols_to_zero), ~ 0))
    # If you need to convert back to a matrix
    predM <- as.matrix(predM)
    rownames(predM) <- colnames(master_imp)
    colnames(predM) <- colnames(master_imp)
    meth <- imp$method
    meth[cols_to_zero] <- ""
    meth["conc"] <- "conc_below_limit"  # Assign custom method
    
    
    message("Retrying MICE imputation (after adjustments)...")
    
    # second attempt
    mice(
      master_imp,
      maxit = 5,
      predictorMatrix = predM,
      method = meth,
      print =  TRUE,
      seed = 123
    )
  })
  
  # check imputated values are not all zero
  complete(imp2, "long") %>%
    filter(censored == TRUE) %>%
    pull(conc) %>%
    summary()

  
  # inspect quality of imputations
  bind_rows(master_imp, 
            complete(imp2, "long")) %>%
    mutate(.imp = ifelse(is.na(.imp), "original", .imp)) %>%
    #visualize density plot, by .imp and censored
    ggplot(aes(x = .imp, y = conc)) +
    geom_jitter(aes(color = censored), width = 0.25, alpha = 0.5) +
    # log transform y axis
    scale_y_sqrt() +
    # rename x-axis to imputation, rename y-axis to concentration
    labs(x = "Imputation", y = "Concentration") +
    # add the title metal.code +
    ggtitle(paste0("Imputation of ", metals[metal.code])) +
    # add a horizontal line at the detection limit 5
    geom_hline(
      yintercept = quantile(master_imp$detect.limit, 0.95, na.rm = TRUE),
      linetype = "dashed"
    ) +
    theme_minimal(base_size = 9)
  
  ggsave(
    paste0(
      "R_Output/",
      metal.code,
      "_imputation_qualitycheck_plot.png"
    ),
    width = 6,
    height = 4,
    dpi = 300
  )
  
  # check the relationship is unchanged
  # master %>%
  #   filter(conc>0 & censored == FALSE) %>% 
  #   mutate(.imp = 0) %>%
  #   bind_rows(complete(imp2, "long"))%>%
  #   ggplot(aes(x = TRI.water.impact, y = log(conc+0.001))) +
  #   geom_smooth(method = "lm", se = FALSE) +
  #   geom_point(aes(color = as.factor(censored)), alpha = 0.5) +
  #   stat_poly_eq(
  #     aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
  #     formula = y ~ x,
  #     parse = TRUE,
  #     label.x = "right",
  #     label.y = "top"
  #   )+
  #   facet_wrap(~.imp) +
  #   labs(x = "TRI", y = "Concentration") +
  #   theme_minimal(base_size = 9)
  # 
  # save imputed data
  saveRDS(complete(imp2, "long"),
          paste0("R_Output/", metal.code, "_imputed_data.rds"))
}

# vectorize across all five metals
purrr::map(metal.codes, impute_missing_values)

# Create Table 1
# read in imputated data, compute summary statistcs, such as number of wells, percent not censorsored, concentration 
# min 25th percentile, median, 75th percentile, and max
# % above MCL
imputed_data <- purrr::map_dfr(metal.codes, ~ readRDS(paste0("R_Output/", .x, "_imputed_data.rds")), .id = "metal") %>%
  mutate(metal = metal.codes[as.integer(metal)]) %>%
  group_by(metal) %>%
  summarise(
    n_wells = n_distinct(location.id),
    n_samples = n(),
    percent_not_censored = sum(!censored) / n() * 100,
    conc_min = min(conc, na.rm = TRUE),
    conc_25th = quantile(conc, 0.25, na.rm = TRUE),
    conc_median = median(conc, na.rm = TRUE),
    conc_75th = quantile(conc, 0.75, na.rm = TRUE),
    conc_max = max(conc, na.rm = TRUE),
    percent_above_MCL = sum(conc > MCLs[metal]) / n() * 100
  )
# save table 1
write.csv(imputed_data, "R_Output/Table1_imputed_data_summary.csv", row.names = FALSE)


##################################################################################
#### Table S5. Reporting (censoring) limits for censored observations ############
##################################################################################
# Responds to Reviewer 1, comment 1a: "Are there multiple censoring limits for
# each of the 5 trace elements? For each constituent, what were the censoring
# limit(s) and number of samples for each limit?"
#
# The tables describe the analytic dataset that enters the MICE imputation above.
# Two sources of the censoring limit are distinguished, because they are not the
# same thing:
#   (i)  "reported"    - the detection/quantitation limit reported with the
#                        record in the Water Quality Portal (detect.limit,
#                        carried through scripts 1-3);
#   (ii) "kNN-filled"  - records with no reported limit, whose detect.limit is
#                        filled in Step 3 of impute_missing_values() together
#                        with the other numeric covariates (VIM::kNN, k = 5,
#                        nearest neighbours in lon/lat). These are NOT reported
#                        limits and are flagged separately below.
# The limit actually used to bound the imputed concentration is the analytic
# limit = reported where available, kNN-filled otherwise. One record = one well
# (script 1 keeps the most recent sample per location.id).
#
# The analytic limits are read back from the saved imputation output
# (R_Output/<metal>_imputed_data.rds), so the counts are exactly those used in
# the imputation. If that file is absent the code falls back to reported limits
# only and warns.

# Assemble reported and analytic limits for one element ------------------------
get_censoring_input <- function(metal.code) {
  # same loading steps as impute_missing_values(), before any imputation
  master <- readRDS(paste0("Data_Files/", metal.code, "_df_PredictorsSelected.rds")) %>%
    dplyr::select(-c(censored.conc, ros.conc))
  if (!"detect.limit" %in% colnames(master)) {
    sup_master <- readRDS(paste0("R_Output/", metal.code, "_RawPredictors.rds")) %>%
      st_drop_geometry() %>%
      dplyr::select(location.id, detect.limit)
    master <- master %>%
      left_join(sup_master, by = "location.id")
  }
  master <- master %>%
    dplyr::select(location.id, conc, censored, detect.limit) %>%
    dplyr::rename(limit_reported = detect.limit)

  imp_file <- paste0("R_Output/", metal.code, "_imputed_data.rds")
  if (file.exists(imp_file)) {
    # first completed dataset; detect.limit is identical across the 5 imputations
    analytic_limits <- readRDS(imp_file) %>%
      filter(.imp == 1) %>%
      dplyr::select(location.id, limit_analytic = detect.limit)
    master <- master %>% left_join(analytic_limits, by = "location.id")
  } else {
    warning(paste0(imp_file, " not found; Table S5 falls back to reported ",
                   "limits and drops censored records without one."))
    master <- master %>%
      mutate(limit_analytic = limit_reported) %>%
      filter(!(censored & is.na(limit_analytic)))
  }

  master %>%
    mutate(
      # guard against floating-point representations of the reported limits
      limit_reported = signif(limit_reported, 6),
      limit_analytic = signif(limit_analytic, 6),
      limit_source   = ifelse(is.na(limit_reported), "kNN-filled", "reported")
    )
}

# One row per element x censoring limit ----------------------------------------
censoring_limit_table <- function(metal.code) {
  master <- get_censoring_input(metal.code)
  n_total <- nrow(master)
  cens    <- master %>% filter(censored)

  cens %>%
    group_by(limit_analytic) %>%
    summarise(
      n_censored     = n(),
      n_limit_reported   = sum(limit_source == "reported"),
      n_limit_knn_filled = sum(limit_source == "kNN-filled"),
      .groups = "drop"
    ) %>%
    arrange(limit_analytic) %>%
    mutate(
      element              = metal.code,
      analyte              = metals[metal.code],
      pct_of_censored      = 100 * n_censored / nrow(cens),
      pct_of_all_samples   = 100 * n_censored / n_total,
      n_samples_analytic   = n_total,
      n_censored_total     = nrow(cens),
      n_distinct_limits    = dplyr::n_distinct(limit_analytic)
    ) %>%
    dplyr::select(
      element, analyte,
      censoring_limit_ugL = limit_analytic,
      n_censored, n_limit_reported, n_limit_knn_filled,
      pct_of_censored, pct_of_all_samples,
      n_samples_analytic, n_censored_total, n_distinct_limits
    )
}

# (a) Element-level summary: one row per element --------------------------------
tableS5_summary <- purrr::map_dfr(metal.codes, function(metal.code) {
  master <- get_censoring_input(metal.code)
  cens   <- master %>% filter(censored)
  dl     <- cens$limit_analytic
  tab    <- sort(table(dl), decreasing = TRUE)
  tibble(
    element                 = metal.code,
    analyte                 = metals[metal.code],
    MCL_or_HBSL_ugL         = MCLs[metal.code],
    n_samples_analytic      = nrow(master),
    n_censored              = nrow(cens),
    pct_censored            = round(100 * nrow(cens) / nrow(master), 1),
    n_distinct_limits       = length(unique(dl)),
    limit_min_ugL           = min(dl),
    limit_median_ugL        = median(dl),
    limit_max_ugL           = max(dl),
    most_common_limit_ugL   = as.numeric(names(tab)[1]),
    pct_at_most_common      = round(100 * as.numeric(tab[1]) / length(dl), 1),
    n_limits_covering_90pct = which(cumsum(as.numeric(tab)) / length(dl) >= 0.90)[1],
    n_censored_limit_knn_filled = sum(cens$limit_source == "kNN-filled")
  )
})

# (b) Full table: every distinct censoring limit, for every element -------------
tableS5_full <- purrr::map_dfr(metal.codes, censoring_limit_table) %>%
  mutate(across(c(pct_of_censored, pct_of_all_samples), ~ round(.x, 2)))

# (c) Condensed table for the SI: limits carrying at least `pool_threshold`
#     percent of an element's censored records are listed individually; the
#     remaining (rare) limits are pooled into a single "Other" row.
pool_threshold <- 1  # percent of censored records
tableS5_condensed <- tableS5_full %>%
  group_by(element, analyte) %>%
  arrange(desc(n_censored), censoring_limit_ugL, .by_group = TRUE) %>%
  group_modify(function(d, key) {
    keep  <- d %>% filter(pct_of_censored >= pool_threshold) %>%
      mutate(censoring_limit_ugL = as.character(censoring_limit_ugL))
    other <- d %>% filter(pct_of_censored <  pool_threshold)
    if (nrow(other) > 0) {
      keep <- bind_rows(keep, tibble(
        censoring_limit_ugL = sprintf(
          "Other (%d limits, %s-%s)", nrow(other),
          format(min(other$censoring_limit_ugL), scientific = FALSE),
          format(max(other$censoring_limit_ugL), scientific = FALSE)
        ),
        n_censored         = sum(other$n_censored),
        n_limit_reported   = sum(other$n_limit_reported),
        n_limit_knn_filled = sum(other$n_limit_knn_filled),
        pct_of_censored    = sum(other$pct_of_censored),
        pct_of_all_samples = sum(other$pct_of_all_samples),
        n_samples_analytic = d$n_samples_analytic[1],
        n_censored_total   = d$n_censored_total[1],
        n_distinct_limits  = d$n_distinct_limits[1]
      ))
    }
    keep %>% mutate(cum_pct_of_censored = round(cumsum(pct_of_censored), 1))
  }) %>%
  ungroup() %>%
  mutate(across(c(pct_of_censored, pct_of_all_samples), ~ round(.x, 1)))

write.csv(tableS5_summary,   "R_Output/TableS5a_censoring_limits_summary.csv", row.names = FALSE)
write.csv(tableS5_condensed, "R_Output/TableS5b_censoring_limits_condensed.csv", row.names = FALSE)
write.csv(tableS5_full,      "R_Output/TableS5c_censoring_limits_full.csv", row.names = FALSE)
