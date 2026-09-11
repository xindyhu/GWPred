############################################################################################
#File name: Predictor_Selection.R
#Author: Jennifer Sun, Cindy Hu
#Date: Jan 2025
#Purpose: Final dataframe cleaning steps, including manual and automatic predictor selection
############################################################################################

source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

select_predictors <- function(metal.code) {
  print(paste0('Begin variable selection for ', metal.code))
  # Read the refined data set
  df <- readRDS(paste0("R_Output/", metal.code, "_RawPredictors.rds"))
  percent_missing <- data.frame(column = names(df),
                                percent_missing = sapply(df, function(x)
                                  mean(is.na(x)) * 100)) %>%
    arrange(desc(percent_missing))
  
  ## 1. Remove unwanted predictor variables -----------------------------------------------------------------------------------------------------------------
  variables_missing_gt10pct = percent_missing %>%
    filter(percent_missing > 10) %>%
    pull(column) %>%
    setdiff("well.depth") # keep well depth
  
  df <- df %>%
    dplyr::select(-all_of(variables_missing_gt10pct)) # variables with > 10% missingness
  
  ## 2. Remove variables that are replicated within raw datasets -------------------------------------------------------------------------------------
  # duplicates from statsgo soil properties; use calculated mean instead
  prefixes <- str_replace(colnames(df), "(L|H|mean)$", "") %>%
    # find the ones that appear three times
    tibble(value = .) %>%
    count(value, name = 'count') %>%
    filter(count == 3) %>%
    pull(value)
  df <- df %>%
    dplyr::select(-c(
      paste0(prefixes, 'L'),
      paste0(prefixes, 'H'),
      'PFLATLOW',
      'PFLATUP'
    ))
  
  # df = subset(df, select=-c(cec_025, cec_05, cec_050, clay_025, clay_05, clay_2550, clay_3060,
  #                           ec_025, ec_05, ksat_05, ph_025, ph_05, ph_2550, ph_3060, sand_025, sand_05,
  #                           sand_2550, sand_3060, silt_025, silt_05, silt_2550, silt_3060, paws_025, paws_050,
  #                           min_ksat, max_ksat)) # duplicates from #16. soil property rasters
  
  
  ## 3. Remove highly correlated variables -----------------------------------------------------------------------------------------------------------------------
  # remove irrelevant variables before running correlation analysis
  df_corData <- df %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c(
      'lon', 'lat', 'conc', 'well.depth', 'DL.missing'
    ))) %>% # keep key vars
    dplyr::select(-any_of(c(
      'TRI.total.impact', 'TRI.water.impact', 'SEMS'
    ))) # anthropogenic input variables that may be correlated but we will keep
  
  # identify highly correlated variables
  df_corMat = cor(df_corData, method = 'pearson') # manually inspect
  # automatic removal of correlated variables - did not use (instead, manually select which of the correlated variables to remove)
  corVars <- caret::findCorrelation(df_corMat, cutoff = 0.9)
  corVars <- colnames(df_corData)[corVars]
  
  # manually selected parameters to remove based on correlation coefficient (view correlation matrix)
  if (metal.code == 'Sr') {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt',
          'no3_pub',
          'AVG_KV'
        )
      ))
  } else if (metal.code == 'Li') {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt',
          'AVG_NO200',
          'no3_pub'
        )
      ))
  } else if (metal.code == 'As') {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt',
          'aq_rocktype'
        )
      ))
  } else if (metal.code %in% c('Cd', 'Mn')) {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt'
        )
      ))
  }
  
  # Save final dataframe
  saveRDS(df,
          paste0("Data_Files/", metal.code, "_df_PredictorsSelected.rds"))
  print(paste0('Variable selection for ', metal.code, ' complete'))
}


# vectorize across all five metals
purrr::map(metal.codes, select_predictors)


############################################################################################
# Depth representativeness of the training wells
# Added September 2026 in response to Reviewer 1: what are the depth ranges of the wells
# in the dataset, what groundwater do the prediction surfaces represent, and how closely
# does the well dataset represent the depths of domestic drinking-water usage?
#
# Outputs, written to R_Output/:
#   TableS9a_depth_distribution.csv        n, min, p25, median, p75, max of well depth,
#                                          overall and by data provider
#   TableS9b_depth_vs_usgs_domestic.csv    training-well depth against the USGS
#                                          domestic-supply depth surfaces, nationally
#   TableS9c_depth_vs_usgs_by_aquifer.csv  the same comparison by principal aquifer
#   FigureS21_depth_vs_domestic.png        depth distribution against the modelled
#                                          domestic-supply interval
#
# UNITS: well depths in the Water Quality Portal are reported in feet (WellDepthMeasure
# MeasureUnitCode is "ft" for every retained record; see 1_All_Data_Prep.R). The USGS
# depth-of-drinking-water grids are also in feet. No unit conversion is applied.
#
# WELL TYPE: the reviewer asked for this broken out by well type. WQP station files carry
# no well-use field -- MonitoringLocationTypeName is the plain value "Well" for >99% of
# the retained sites -- and 1_All_Data_Prep.R does not retain it in any case. Recovering
# domestic / public / monitoring status would require the NWIS site service
# (dataRetrieval::readNWISsite, field well_use_cd) and would cover only the NWIS sites.
# The breakout below is therefore by data provider, which is what actually determines
# whether a depth is reported at all.
############################################################################################

library(terra)
library(sf)
library(ggplot2)

## ---------------------------------------------------------------------------------------
## 1. Assemble the analytic wells (the datasets that 4_MICE.R consumes)
## ---------------------------------------------------------------------------------------
wells <- purrr::map_dfr(metal.codes, function(m) {
  d <- readRDS(paste0("Data_Files/", m, "_df_PredictorsSelected.rds"))
  data.frame(
    element     = m,
    location.id = as.character(d$location.id),
    long        = d$long,
    lat         = d$lat,
    well.depth  = d$well.depth,       # feet
    data.source = as.character(d$data.source),
    aq_code     = suppressWarnings(as.integer(as.character(d$aquifer))),
    stringsAsFactors = FALSE
  )
})

## ---------------------------------------------------------------------------------------
## 2. Table S9a -- depth distribution
## ---------------------------------------------------------------------------------------
depth_summary <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      n_wells    = n(),
      n_depth    = sum(!is.na(well.depth)),
      pct_depth  = 100 * mean(!is.na(well.depth)),
      min_ft     = suppressWarnings(min(well.depth, na.rm = TRUE)),
      p25_ft     = quantile(well.depth, 0.25, na.rm = TRUE, names = FALSE),
      median_ft  = median(well.depth, na.rm = TRUE),
      p75_ft     = quantile(well.depth, 0.75, na.rm = TRUE, names = FALSE),
      max_ft     = suppressWarnings(max(well.depth, na.rm = TRUE)),
      .groups    = "drop"
    ) %>%
    mutate(across(c(min_ft, max_ft), ~ ifelse(is.finite(.x), .x, NA_real_)))
}

tabS9a <- bind_rows(
  depth_summary(wells, element) %>% mutate(stratum = "All wells", .after = element),
  depth_summary(wells, element, data.source) %>% rename(stratum = data.source)
) %>% arrange(element, stratum)

write.csv(tabS9a, "R_Output/TableS9a_depth_distribution.csv", row.names = FALSE)
print(as.data.frame(tabS9a), digits = 4)

## ---------------------------------------------------------------------------------------
## 3. USGS depth of groundwater used for drinking-water supplies (ver. 2.0, January 2026)
##    Kauffman, Degnan, Belitz, Stackelberg and Erickson, https://doi.org/10.5066/P94640EM
##    Downloaded once and cached under CoVar/.
## ---------------------------------------------------------------------------------------
usgs_dir <- "CoVar/USGS_DepthDrinkingWater"
dir.create(usgs_dir, showWarnings = FALSE, recursive = TRUE)
usgs_zip <- file.path(usgs_dir, "domestic_grids.zip")
usgs_url <- paste0("https://www.sciencebase.gov/catalog/file/get/5e43efc3e4b0edb47be84c3d",
                   "?f=__disk__4f%2F4e%2F22%2F4f4e22143258d3c373dd1f9e3841afae8b88db91")

if (!file.exists(file.path(usgs_dir, "domestic_top_open.asc"))) {
  if (!file.exists(usgs_zip)) {
    options(timeout = max(1800, getOption("timeout")))
    download.file(usgs_url, usgs_zip, mode = "wb")   # ~82 MB
  }
  unzip(usgs_zip, files = c("domestic_top_open.asc", "domestic_bottom_open.asc"),
        exdir = usgs_dir)
}

## The .prj files are old-style ESRI ALBERS descriptors that terra does not parse; the
## parameters they contain (Albers, NAD83/GRS80, standard parallels 29 30 and 45 30,
## central meridian -96, latitude of origin 23, metres) are USGS CONUS Albers = EPSG:5070.
crs_5070 <- "EPSG:5070"
top_open <- terra::rast(file.path(usgs_dir, "domestic_top_open.asc"))
bot_open <- terra::rast(file.path(usgs_dir, "domestic_bottom_open.asc"))
terra::crs(top_open) <- crs_5070
terra::crs(bot_open) <- crs_5070

pts <- sf::st_as_sf(wells, coords = c("long", "lat"), crs = 4269, remove = FALSE) %>%
  sf::st_transform(crs_5070) %>%
  terra::vect()

wells$usgs_top_ft <- terra::extract(top_open, pts)[, 2]
wells$usgs_bot_ft <- terra::extract(bot_open, pts)[, 2]

## ---------------------------------------------------------------------------------------
## 4. Table S9b -- national comparison
##
## usgs_bot_ft is a moving MEDIAN of the depth to the bottom of the open interval of
## domestic-supply wells. If the training wells sampled the same depth population, about
## half of them would be deeper than it at their own location; that is the interpretable
## statistic. The fraction falling strictly inside [top, bottom] is also reported, but the
## modelled open interval is thin (national median thickness ~45 ft) relative to the spread
## of well depths, so a low value there is not by itself evidence of bias.
## ---------------------------------------------------------------------------------------
cmp <- wells %>%
  filter(!is.na(well.depth), !is.na(usgs_top_ft), !is.na(usgs_bot_ft),
         usgs_bot_ft > usgs_top_ft)

tabS9b <- cmp %>%
  group_by(element) %>%
  summarise(
    n                 = n(),
    p25_ft     = quantile(well.depth, 0.25, na.rm = TRUE, names = FALSE),
    median_ft  = median(well.depth, na.rm = TRUE),
    p75_ft     = quantile(well.depth, 0.75, na.rm = TRUE, names = FALSE),
    median_usgs_top   = median(usgs_top_ft),
    median_usgs_bot   = median(usgs_bot_ft),
    pct_within        = 100 * mean(well.depth >= usgs_top_ft & well.depth <= usgs_bot_ft),
    pct_shallower     = 100 * mean(well.depth <  usgs_top_ft),
    pct_deeper        = 100 * mean(well.depth >  usgs_bot_ft),
    median_diff_ft    = median(well.depth - usgs_bot_ft),
    median_ratio      = median(well.depth / usgs_bot_ft),
    pct_gt_2x_bottom  = 100 * mean(well.depth > 2 * usgs_bot_ft),
    pct_gt_1000ft     = 100 * mean(well.depth > 1000),
    .groups = "drop"
  )
write.csv(tabS9b, "R_Output/TableS9b_depth_vs_usgs_domestic.csv", row.names = FALSE)
print(as.data.frame(tabS9b), digits = 4)

## ---------------------------------------------------------------------------------------
## 5. Table S9c -- by principal aquifer
##    Pooled over unique wells; a well contributing to more than one element is counted once.
## ---------------------------------------------------------------------------------------
aq_lookup <- sf::read_sf(dsn = "CoVar/aquifrp025_nt00003", layer = "aquifrp025") %>%
  sf::st_drop_geometry() %>%
  transmute(aq_code = as.integer(AQ_CODE), aq_name = as.character(AQ_NAME)) %>%
  distinct(aq_code, .keep_all = TRUE)

tabS9c <- cmp %>%
  distinct(location.id, .keep_all = TRUE) %>%
  left_join(aq_lookup, by = "aq_code") %>%
  group_by(aq_code, aq_name) %>%
  summarise(
    n               = n(),
    median_well_ft  = median(well.depth),
    median_usgs_top = median(usgs_top_ft),
    median_usgs_bot = median(usgs_bot_ft),
    pct_within      = 100 * mean(well.depth >= usgs_top_ft & well.depth <= usgs_bot_ft),
    pct_deeper      = 100 * mean(well.depth >  usgs_bot_ft),
    median_diff_ft  = median(well.depth - usgs_bot_ft),
    .groups = "drop"
  ) %>%
  filter(n >= 100) %>%
  arrange(desc(n))
write.csv(tabS9c, "R_Output/TableS9c_depth_vs_usgs_by_aquifer.csv", row.names = FALSE)
print(as.data.frame(tabS9c), digits = 4)

## ---------------------------------------------------------------------------------------
## 6. Figure S21
## ---------------------------------------------------------------------------------------
fig_df <- cmp %>% filter(well.depth > 0)
band <- cmp %>%
  summarise(top = median(usgs_top_ft), bot = median(usgs_bot_ft))

p <- ggplot(fig_df, aes(x = well.depth)) +
  annotate("rect", xmin = band$top, xmax = band$bot, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.45) +
  geom_density(aes(colour = element), linewidth = 0.6) +
  scale_x_log10(breaks = c(10, 30, 100, 300, 1000, 3000, 10000),
                labels = scales::comma) +
  labs(x = "Well depth (ft, log scale)", y = "Density", colour = "Element",
       caption = paste0("Shaded band: median modelled domestic-supply open interval ",
                        "(", round(band$top), "-", round(band$bot), " ft) at the ",
                        "training-well locations.")) +
  theme_bw(base_size = 11)
ggsave("R_Output/FigureS21_depth_vs_domestic.png", p, width = 7, height = 4.5, dpi = 300)

saveRDS(wells, "R_Output/wells_depth_representativeness.rds")
cat("Depth representativeness analysis complete.\n")
