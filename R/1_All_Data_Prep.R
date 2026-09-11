###########################################################################
#File name: All_Data_Prep.R
#Author: Jennifer Sun, Jonas LaPier, Cindy Hu
#Date: Jan 2025
#Purpose: All steps for WQP well data preparation
###########################################################################

#### Attach Libraries and Set Working Directory ####
source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))
mapUSm <- st_read("CoVar/US48/US_48states.shp") # load projected map of the US (in m)

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes
MCLs <- c(10, 5, 60, 300, 4000)
names(MCLs) <- metal.codes

## MCL Values ----
# Arsenic -> 10 ug/L
# Cadmium -> 5 ug/L
# Lead -> 15 ug/L (action level for drinking water)
# Manganese -> 50 ug/L (secondary standard)
# Uranium -> 30 ug/L (MCL)
# Strontium -> 4000 ug/L (HBSL)
# Lithium -> 60 ug/L (drinking water HBSL)

data_prep_function <- function(metal.code) {
  ## 1. Load WQP Files -----------------------------------------------------------------------------------------------------------------
  #  Read in the station and sample result csv files obtained from the WQP. Locations outside the contiguous US are excluded.
  
  ## Load Data
  metal.stn = read.csv(paste0("WQP/", metal.code, "-station.csv"))
  metal.smpl = read.csv(paste0("WQP/", metal.code, "-result.csv"))
  # Merge files based on location identifier
  metal.large = metal.smpl |>
    left_join(
      metal.stn,
      by = c(
        'OrganizationIdentifier',
        'OrganizationFormalName',
        'MonitoringLocationIdentifier',
        'ProviderName'
      )
    )
  
  ## 2. Format WQP data ------------------------------------------------------------------------------------------------------------------
  # Rename columns and match up data types and units
  
  metal.df = metal.large |>
    # Keep only important columns (all well depth units are in feet)
    dplyr::select(
      "MonitoringLocationIdentifier",
      "MonitoringLocationTypeName",
      "LatitudeMeasure",
      "LongitudeMeasure",
      "AquiferName",
      "AquiferTypeName",
      "WellDepthMeasure.MeasureValue",
      "ActivityIdentifier",
      "ActivityStartDate",
      "ActivityTypeCode",
      "ResultMeasureValue",
      "ResultMeasure.MeasureUnitCode",
      "ResultDetectionConditionText",
      "DetectionQuantitationLimitMeasure.MeasureValue",
      "DetectionQuantitationLimitMeasure.MeasureUnitCode",
      "ActivityMediaSubdivisionName",
      'ResultSampleFractionText',
      'ResultStatusIdentifier',
      'ResultValueTypeName',
      'ResultLaboratoryCommentText',
      'DetectionQuantitationLimitTypeName',
      'ProviderName'
    ) |>
    # Rename columns
    dplyr::rename(
      "location.id" = "MonitoringLocationIdentifier",
      "location.type" = "MonitoringLocationTypeName",
      "latitude" = "LatitudeMeasure",
      "longitude" = "LongitudeMeasure",
      "aquifer.name" = "AquiferName",
      "aquifer.type" = "AquiferTypeName",
      "well.depth" = "WellDepthMeasure.MeasureValue",
      "activity" = "ActivityIdentifier",
      "date" = "ActivityStartDate",
      "activity.type" = "ActivityTypeCode",
      "conc" = "ResultMeasureValue",
      "conc.unit" = "ResultMeasure.MeasureUnitCode",
      "result.condition" = "ResultDetectionConditionText",
      "detect.limit" = "DetectionQuantitationLimitMeasure.MeasureValue",
      "detect.unit" = "DetectionQuantitationLimitMeasure.MeasureUnitCode",
      "sample.media" = "ActivityMediaSubdivisionName",
      "sample.fraction" = "ResultSampleFractionText",
      "result.status" = "ResultStatusIdentifier",
      "result.type" = "ResultValueTypeName",
      "result.comments" = "ResultLaboratoryCommentText",
      "DL.type" = "DetectionQuantitationLimitTypeName",
      "data.source" = "ProviderName"
    ) |>
    # convert conc, well.depth, detect.limit from character to numeric
    mutate(
      conc = as.numeric(as.character(conc)),
      well.depth = as.numeric(as.character(well.depth)),
      detect.limit = as.numeric(as.character(detect.limit))
    ) |>
    # add a year column, and convert date to date type
    mutate(date = lubridate::ymd(date), year = lubridate::year(date)) |>
    # if conc.unit is mg/l or ppm, convert conc to 1000*conc, then reset conc.unit to ug/L
    mutate(
      conc = ifelse(conc.unit %in% c("mg/l", "ppm"), 1000 * conc, conc),
      conc.unit = ifelse(conc.unit %in% c("mg/l", "ppm"), "ug/l", conc.unit)
    ) |>
    # if detect.unit is mg/l or ppm, convert detect.limit to 1000*detect.limit, then reset detect.unit to ug/L
    mutate(
      detect.limit = ifelse(
        detect.unit %in% c("mg/l", "ppm"),
        1000 * detect.limit,
        detect.limit
      ),
      detect.unit = ifelse(detect.unit %in% c("mg/l", "ppm"), "ug/l", detect.unit)
    )
  
  #### 3. Data Cleaning ---------------------------------------------------------------------------------------------------------------------
  metal.df1 <- metal.df |>
    # subset to samples collected in 1990 and later
    filter(year >= 1990) |>
    # remove concentrations that are negative
    filter(conc >= 0 | is.na(conc)) |>
    # Remove rows with other units (i.e. mg/kg, pCi/L)
    filter(conc.unit == "ug/l" |
             (conc.unit == "" & detect.unit == "ug/l"))  |>
    # keep rows that match string "Well" in location.type
    filter(grepl("Well", location.type)) |>
    # keep rows where activity.type is Sample or Sample-Routine
    filter(activity.type %in% c("Sample", "Sample-Routine")) |>
    # Remove records that are not labeled as groundwater in sample.media column
    filter(sample.media %in% c("Ground Water", "Groundwater")) |>
    # Remove records that are the wrong sample fraction type
    filter(!sample.fraction %in% c("Suspended", "", "Bed Sediment")) |>
    # remove records with systematic contamination
    filter(!(result.condition %in% c('Systematic Contamination')))
  # Cut to values within CONUS
  metal.df1 <- st_as_sf(metal.df1,
                        coords = c("longitude", "latitude"),
                        crs = st_crs(mapUSm))# setting initial crs
  # Determine if each point is within the polygon
  within_polygon <- st_contains(st_union(mapUSm), metal.df1, sparse = FALSE)
  metal.df1$insideUS <- within_polygon[1, ]
  metal.df2 <- metal.df1[metal.df1$insideUS == TRUE, ]
  #drop geometry
  metal.df2 <- st_drop_geometry(metal.df2) |>
    as.data.frame()
  ## 4. Detection Limit Handling -------------------------------------------------------------------------------------------------------------------
  metal.df2 <- metal.df2 %>%
    # remove records where detection limit is below 0 or above the MCL
    filter(is.na(detect.limit) |
             (detect.limit > 0 & detect.limit < MCLs[metal.code])) |>
    mutate(DL.missing = as.numeric(detect.limit == 0 |
                                     is.na(detect.limit)))
  
  print(metal.code)
  print(paste0('Unique wells: ', length(unique(
    metal.df2$location.id
  ))))
  print(paste0('Total records: ', length(metal.df2$conc)))
  
  ## 5. Select well value (merge multiple measurements at the same well location) -----------------------------------------------------------------------------
  # select most recent value (and use median value when there are multiple measurements on the same day)
  metal.df3 <- metal.df2 |>
    group_by(location.id) |>
    arrange(desc(date)) |> # select data from the most recent date for each well
    slice(1) |>
    ungroup() |>
    # keep relevant columns
    dplyr::select(c(
      location.id,
      date,
      conc,
      well.depth,
      location.type,
      sample.fraction,
      data.source,
      DL.missing,
      detect.limit
    )) |>
    # join back lat/long
    left_join(metal.df %>% dplyr::select(location.id, latitude, longitude) %>% distinct(),
              by = 'location.id')
  
  
  ## 6. Save cleaned data ---------------------------------------------------------------------------------------------------------------------------
  saveRDS(metal.df3,
          paste0("Data_Files/", metal.code, "_Combined_WQP.rds"))
  print(paste("WQP Data Saved for", metal.code))
}

# vectorize across all five metals
purrr::map(metal.codes, data_prep_function)
