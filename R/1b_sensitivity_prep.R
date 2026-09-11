suppressMessages({library(dplyr); library(readr); library(jsonlite); library(data.table)})
UP <- "/mnt/user-data/uploads/GWPred/data/R_Output"
OUT <- "/tmp/gw/out"
metals <- c("As","Cd","Li","Mn","Sr")
diss_set <- c("Dissolved","Filtered, lab","Filterable","Pot. Dissolved")
remove_list <- c("date","conc","DL.missing","detect.limit","censored","is.imputed",
                 "data.source","lon","lat","location.id","well.depth",".imp",".id")
factor_cols_all <- c("KB","surfgeo","lith","landcover_500m","aquifer","aq_rocktype",
                     "drainage_class_int","hydgrp_int","soilorder_int","str_int","weg_int")

get_best_model <- function(csv){
  read_csv(csv, show_col_types=FALSE) %>%
    filter(mean_val < min_plus_sd[1]) %>%
    slice_max(mean_val, n=1) %>%
    dplyr::select(trees, tree_depth, learn_rate, sample_size, loss_reduction) %>%
    slice(1)
}

for (m in metals){
  df <- readRDS(file.path(UP, paste0(m,"_imputed_data.rds")))
  df$logconc <- log10(df$conc)
  df$is.imputed <- as.integer(df$censored)
  # fraction join
  fm <- fread(file.path(UP,"frac_map", paste0(m,"_fraction_map.csv")))
  fm <- as.data.frame(fm)[,c("location.id","sample.fraction")]
  n_wells <- length(unique(df$location.id))
  df <- df %>% left_join(fm, by="location.id")
  matched <- df %>% distinct(location.id, sample.fraction)
  n_unmatched <- sum(is.na(matched$sample.fraction))
  df$frac_group <- ifelse(is.na(df$sample.fraction), NA,
                     ifelse(df$sample.fraction %in% diss_set, "dissolved","total"))
  # factor levels present in this element
  fl <- list()
  for (fc in factor_cols_all){
    if (fc %in% names(df)) fl[[fc]] <- levels(df[[fc]])
  }
  write(toJSON(fl, auto_unbox=TRUE), file.path(OUT, paste0(m,"_factlevels.json")))
  hp <- get_best_model(file.path(UP, paste0(m,"_df_hyperparameter.csv")))
  write(toJSON(as.list(hp), auto_unbox=TRUE), file.path(OUT, paste0(m,"_hp.json")))
  # convert date to char, factors to character for csv
  df$date <- as.character(df$date)
  for (fc in factor_cols_all) if (fc %in% names(df)) df[[fc]] <- as.character(df[[fc]])
  fwrite(df, file.path(OUT, paste0(m,"_analytic.csv.gz")))
  # coverage per well (imp1)
  w1 <- df[df$.imp==1,]
  fg <- table(w1$frac_group, useNA="ifany")
  cat(sprintf("%s: wells=%d unmatched_frac=%d | imp1 frac_group: %s | hp: trees=%s depth=%s lr=%s subsamp=%s gamma=%s\n",
      m, n_wells, n_unmatched, paste(names(fg),fg,sep="=",collapse=" "),
      hp$trees, hp$tree_depth, hp$learn_rate, hp$sample_size, hp$loss_reduction))
}
