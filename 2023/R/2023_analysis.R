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
bts_biomass(year=year, rmv_yrs = c(1984,1987))

# run projection
proj_ak(year=year, last_full_assess=2022, folder="harvest_proj", species="dusky", region="goa", off_yr=TRUE)



# setup for next year
setup_folders(year+1)
accepted_model(2022, "m22.3a", year +1)
