## SARA FILE##
library(tidyverse)
library(here)

AYR = 2024
model = "m22.5a"

SARAdir <- here::here(AYR, 'data', 'sara')

if(file.exists(SARAdir)==F) {
  dir.create(paste0(SARAdir), showWarnings = T)
}

## Files and data sets used in SARA file

exec_sum <- read.csv(here::here(AYR, model, 'processed', 'exec_summ.csv'))
bio_rec_f <- read.csv(here::here(AYR, model, 'processed', 'bio_rec_f.csv'))
mcmc_dat <- read.csv(here::here(year, model, "tables", "t_recr_biom_mcmc.csv"))
survey_b <- read.csv(here::here(AYR, model, 'processed', 'survey.csv'))
catch <- read.csv(here::here(AYR, 'data', 'output', 'fish_catch.csv'))
surv_yrs_bind <- data.frame(year = 1977:2023)

best_f <- 0.13678918568669

#### Static SARA info ----
stock <- "DUSKY"
stock_name <- "Dusky rockfish - Gulf of Alaska"
fmp <- "GOA"
AYR <- AYR #assessment year
SYR <- 2023
asmt_type <- "Operational full"
asmt_mo <- "Dec"
tier <- "3a"
num_sexes <- 1
num_fish <- 1
rec_mult <- 1000000
recage <- 4
complex <- "NA"
asmt_mod_cat <- "6 - Statistical Catch-at-Age"
asmt_mod <- "Custom SCAA: Custom Statistical Catch-at-Age Model"
mod_version <- "22.5a"
lead_lab <- "AFSC"
email <- "kristen.omori@noaa.gov"
review_result <- "Full acceptance" #
catch_input_dat <- 4 #options: 0 - None, 1 - Major gaps preclude use, 2 - Major gaps in some sectors(s), 3 - Minor gaps across sectors, 4 - Minor gaps in some sector(s), 5 - Near complete knowledge)
abund_input_dat <- 3 # options: 0 - None, 1 - Uncertain or expert opinion, 2 - Standardized fishery-dependent, 3 - Limited fishery-independent, 4 - Comprehensive fishery-independent, 5 - Absolute abundance)
biol_input_dat <- 2 # options: 0 - None, 1 - Proxy-based, 2 - Empirical and proxy-based, 3 - Mostly empirical estimates, 4 - Track changes over time, 5 - Comprehensive over time and space)
sizeage_comp_input_dat <- 4 # options: 0 - None, 1 - Major gaps preclude use, 2 - Support data-limited only, 3 - Gaps, bus supports age-structured assessment, 4 - Support fishery composition, 5 - Very complete)
ecosys_link <- 2

# F section
Fbasis <- "F for Fully-Selected Fish"
Funit <- "Fully-selected F"
Flim_basis <- "F35% as proxy"  # trying to figure out if this is FOFL or FABC
Fmsy_basis <- "F35% as proxy" # trying to figure out if this is FOFL or FABC

# Biomass section
b_basis <- "Mature Female Biomass"
b_unit <- "Metric Tons"
est_method <- "Credible" # (options Tiers 1-3 only: Asymptotic, Credible, Bootstrapped, "NA" for Tiers 4-6)
interval_size <- "95" # Specify size of confidence interval (options Tiers 1-3 only: 50 to 99, "NA" for Tiers 4-6)
b_msy_basis <- "B35%" # (options Tiers 1-3 only: Direct estimate, S_MSY escapement, Average Survey CPUE, B40%, B35%, B30%, "NA" for Tiers 4-6)

# other wording
surv_desc <- "Gulf of Alaska Shelf and Slope Groundfish Bottom Trawl" # ("NA" for Tier 6)

stock_status <- "SAFE report indicates that this stock was not subjected to overfishing in 2023 or is being overfished in 2024."

# numbers
Fest <- bio_rec_f %>% filter(year %in% (AYR-1)) %>% mutate(F = round(F, 3)) %>% pull(F)  # last full F rate
Flim <- exec_sum %>% filter(item == "FOFL") %>% mutate(y1 = round(as.numeric(y1),3) ) %>% pull(y1) # OFL spec year (not new rec spec)
Fmsy <- Flim # OFL last year

b_est <- bio_rec_f %>% filter(year %in% (AYR)) %>% mutate(sp_biom = round(sp_biom, 0)) %>% pull(sp_biom)
b_lci <- mcmc_dat %>% filter(Year %in% (AYR)) %>% mutate(sp_lci = round(sp_lci,0)) %>% pull(sp_lci)
b_uci <-  mcmc_dat %>% filter(Year %in% (AYR)) %>% mutate(sp_uci = round(sp_uci,0)) %>% pull(sp_uci)
b_msy <- exec_sum %>% filter(item == "B35%") %>% pull(y3) # (Tiers 1-3 only, "NA" for Tiers 4-6)

# time series
fish_yrs <- bio_rec_f %>% pull(year)
F_yrs <- bio_rec_f %>% pull(year)
ages <- recage:33

recruit <- bio_rec_f %>% mutate(recruits = round(recruits, 4)) %>% pull(recruits)
spawn_biom <- bio_rec_f %>% mutate(sp_biom = round(sp_biom, 0)) %>% pull(sp_biom)
tot_biom <- bio_rec_f %>% mutate(tot_biom = round(tot_biom, 0)) %>% pull(tot_biom)
F_mort <- bio_rec_f %>% mutate(F = round(F, 4)) %>% pull(F)
tot_catch <- catch %>% mutate(catch = round(catch, 0)) %>% pull(catch)
surv_biom <- survey_b %>% mutate(biomass = round(biomass)) %>% full_join(surv_yrs_bind) %>% arrange(year) %>% pull(biomass)

# Compile the .dat file ----
cat(
  "#Stock Assessment Results Archive (SARA) file for stocks managed under the North Pacific Fisheries Management Council", "\n",
  "#This form is only required for Tier 1-6 stocks in a full year, or Tier 1-3 stocks in a partial year", "\n",
  "#There are four required sections to update: Assessment Summary, Fishing Mortality Estimates, Biomass Estimates, and Time Series Estimates with an optional fifth Survey Estimates section", "\n",
  "#Please fill in the text (surrounded by", '" "', ') or data as values in the line after each field marked with a # and capitalized name (e.g., #STOCK, the next line should be', paste0("\"", "AFT","\""), "\n",
  "#Note that several of the fields below may be static or only change occasionally, please see the SARA_HQ_Lookup Google sheet for more details", "\n",
  "#ASSESSMENT_SUMMARY -------------------------------------------------------------------------------------------------", "\n",
  "#STOCK","\n",
  paste0("\"", stock,"\""), "\n",
  "#STOCK_NAME","\n",
  paste0("\"", stock_name,"\""), "\n",
  "#REGION","\n",
  paste0("\"", fmp,"\""), "\n",
  "#ASMT_TYPE", "\n",
  paste0("\"", asmt_type,"\""), "\n",
  "#ASMT_YEAR", "\n",
  AYR, "\n",
  "#ASMT_MONTH", "\n",
  paste0("\"", asmt_mo,"\""), "\n",
  #'"', asmt_mo, '"', "\n",
  "#TIER", "\n",
  paste0("\"", tier,"\""), "\n",
  #'"', tier, '"', "\n",
  "#NUM_SEXES", "\n",
  num_sexes, "\n",
  "#NUM_FISHERIES", "\n",
  num_fish, "\n",
  "#REC_MULT", "\n",
  format(rec_mult, scientific = F), "\n",
  "#RECAGE", "\n",
  recage, "\n",
  "#COMPLEX", "\n",
  paste0("\"", complex,"\""), "\n",
  "#LAST_DATA_YEAR", "\n",
  AYR, "\n",
  "#ASMT_MODEL_CATEGORY", "\n",
  paste0("\"", asmt_mod_cat,"\""), "\n",
  "#ASMT_MODEL", "\n",
  paste0("\"", asmt_mod,"\""), "\n",
  "#MODEL_VERSION", "\n",
  paste0("\"", mod_version,"\""), "\n",
  "#ENSEMBLE", "\n",
  paste0("\"", "NA","\""), "\n",
  "#LEAD_LAB", '"', "\n",
  paste0("\"", lead_lab,"\""), "\n",
  "#POC_EMAIL", "\n",
  paste0("\"", email,"\""), "\n",
  "#REVIEW_RESULT", "\n",
  paste0("\"", review_result,"\""), "\n",
  "#CATCH_INPUT_DATA", "\n",
  catch_input_dat, "\n",
  "#ABUNDANCE_INPUT_DATA", "\n",
  abund_input_dat, "\n",
  "#BIOLOGICAL_INPUT_DATA", "\n",
  biol_input_dat, "\n",
  "#SIZEAGE_COMP_INPUT_DATA", "\n",
  sizeage_comp_input_dat, "\n",
  "#ECOSYSTEM_LINKAGE", "\n",
  ecosys_link, "\n",
  "#FISHING_MORTALITY_ESTIMATES ----------------------------------------------------------------------------------------", "\n",
  "#F_YEAR", "\n",
  AYR-1, "\n",
  "#F_BASIS", "\n",
  paste0("\"", Fbasis,"\""), "\n",
  "#F_UNIT", "\n",
  paste0("\"", Funit,"\""), "\n",
  "#BEST_F_ESTIMATE", "\n",
  Fest, "\n",
  "#F_LIMIT", "\n",
  Flim, "\n",
  "#F_LIMIT_BASIS", "\n",
  paste0("\"", Flim_basis,"\""), "\n",
  "#F_MSY", "\n",
  Fmsy, "\n",
  "#F_MSY_BASIS", "\n",
  paste0("\"", Fmsy_basis,"\""), "\n",
  "#BIOMASS_ESTIMATES --------------------------------------------------------------------------------------------------", "\n",
  "#B_YEAR", "\n",
  AYR, "\n",
  "#B_BASIS", "\n",
  paste0("\"", b_basis,"\""), "\n",
  "#B_UNIT", "\n",
  paste0("\"", b_unit,"\""), "\n",
  "#BEST_B_ESTIMATE", "\n",
  b_est, "\n",
  "#LOWER_B_ESTIMATE", "\n",
  b_lci, "\n",
  "#UPPER_B_ESTIMATE", "\n",
  b_uci, "\n",
  "#ESTIMATE_METHOD", "\n",
  paste0("\"", est_method,"\""), "\n",
  "#INTERVAL_SIZE", "\n",
  interval_size, "\n",
  "#B_MSY", "\n",
  b_msy, "\n",
  "#B_MSY_BASIS", "\n",
  paste0("\"", b_msy_basis,"\""), "\n",
  "#TIME_SERIES_ESTIMATES  ----------------------------------------------------------------------------------------------", "\n",
  "#FISHERYYEAR", "\n",
  fish_yrs, "\n",
  "#AGE", "\n",
  ages, "\n",
  "#RECRUITMENT", "\n",
  recruit, "\n",
  "#SPAWNBIOMASS", "\n",
  spawn_biom, "\n",
  "#TOTALBIOMASS", "\n",
  tot_biom, "\n",
  "#TOTFSHRYMORT", "\n",
  F_mort, "\n",
  "#TOTALCATCH", "\n",
  tot_catch, "\n",
  "#SURVEYDESC", "\n",
  paste0("\"", surv_desc,"\""), "\n",
  "#STOCKNOTES", "\n",
  stock_status, "\n",
  "#SURVEY_ESTIMATES [OPTIONAL] ------------------------------------------------------------------------------------------", "\n",
  surv_biom,
  file=paste0(SARAdir,"/DUSKGOA",AYR,"_HQ.dat"))

# End
