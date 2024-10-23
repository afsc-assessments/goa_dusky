## Data base queries on the side ##
## July 16, 2024

## load ----
library(tidyverse)
library(here)
library(keyring)
## data and switches ----

yr <- 2024
species_code_race <- c(30152, 30150) # dusky (and dark_rockfish) & dusky and dark rockfish
species_name <- "Dusky"
sci_name <- "Sebastes_variabilis"
# design-based (GAP) index from bottom-trawl survey
## Need to run the GAP index to get the subareas

# Old way to get GAP index
keyring::key_list()
#db_afsc = afscdata::connect("afsc")
db_akfin = afscdata::connect()

db_race_goa_biom_area <- dplyr::tbl(db_akfin, dplyr::sql("afsc.race_biomassareaaigoa")) %>%
  dplyr::rename_with(tolower) %>%
  dplyr::filter(year <= yr,
                species_code %in% c(species_code_race),
                survey == "GOA" ) %>%
  dplyr::mutate(area = case_when(regulatory_area_name == "WESTERN GOA" ~ "Western",
                                 regulatory_area_name == "CENTRAL GOA" ~ "Central",
                                 regulatory_area_name == "EASTERN GOA" ~ "Eastern")) %>%
  dplyr::group_by(year, area, haul_count) %>%
  dplyr::summarise(n_catch = sum(catch_count, na.rm=T),
                   biomass = sum(area_biomass, na.rm=T),
                   var = sum(biomass_var, na.rm=T ) ) %>%
  dplyr::collect()

db_race_goa_biom_total <- db_race_goa_biom_area %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(haul_count = sum(haul_count, na.rm=T),
                   n_catch = sum(n_catch, na.rm=T),
                   biomass = sum(biomass, na.rm=T),
                   var = sum(var, na.rm=T) ) %>%
  dplyr::mutate(area = "Total")

db_race_goa_biom <- full_join(db_race_goa_biom_area, db_race_goa_biom_total) %>%
  dplyr::arrange(year, area)

#write.csv(db_race_goa_biom, here::here(year, 'data','raw','goa_area_bts_biomass_data.csv'), row.names= F )

# VAST
libs_vast <- c("INLA", "TMB", "DHARMa", "FishStatsUtils", "VAST")
lapply(libs_vast, library, character.only = TRUE)
