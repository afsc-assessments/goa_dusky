# load ----
# devtools::unload("afscassess")
library(afscdata)
library(afscassess)

# globals ----
year = 2023
species = "DUSK"
area = "goa"
TAC <- c(3676, 5389, 5372)


# setup
setup_folders(year)
accepted_model(2022, "m22.3a", year)

# data ----
goa_dusk(year, off_yr = TRUE)
clean_catch(year, species=species, TAC=TAC)
bts_biomass(year=year)
