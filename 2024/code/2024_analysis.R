# 2024 dusky rockfish assessment
# kristen.omori@noaa.gov
# based on ben.williams@noaa.gov

### load ----
#devtools::unload("afscdata")
#devtools::unload("afscassess")
#devtools::install_github("afsc-assessments/afscdata", force = TRUE)
#devtools::install_github("BenWilliams-NOAA/afscassess", force = TRUE)
#https://afsc-assessments.github.io/sop/data.html#passwords
#devtools::install_github(repo = "afsc-assessments/spmR", force = TRUE)

library(afscdata)
library(afscassess)
library(tidyverse)
library(ggplot2)
library(here)
library(keyring)
library(rema)
library(R2admb)
#library(spmR)
theme_set(afscassess::theme_report())

### globals ----

# static
species <- "DUSK"
area = "goa"
afsc_species1 =  30150
afsc_species2 = 30152
norpac_species = 330
cas_code = c(154, 172)
species_group_code = c("PEL7", "PELS", "DUSK")
fishery = "fish"
dat_name = "goa_dusk"

'%nin%' <- Negate('%in%')

keyring::key_list()
akfin_user = keyring::key_list("akfin")$username
akfin_pwd = keyring::key_get("akfin", keyring::key_list("akfin")$username)
afsc_user = keyring::key_list("afsc")$username
afsc_pwd = keyring::key_get("afsc", keyring::key_list("afsc")$username)

# adjust/ update
year = 2024
TAC_yrs <- c(year-3, year-2, year-1)
TAC <- c(5389, 5372, 7917) # previous 3 years
OFL_2023 <-9638
# Can check TAC values with AKFIN: AKR.V_CAS_TAC
prev_model_name = "m22.3a"

rec_age = 4
plus_age = 30
len_bins = read.csv(here::here(year, "data", "user_input", "lbins.csv"))
len_bins = min(len_bins$len_bins):max(len_bins$len_bins) # use csv 'lbins.csv'

# source temporary functions
# data functions: stores only concat_dat_dusk() functions
source(here::here(year, "code", "data_functions.R"))
source(here::here(year, "code", "functions_2024.R"))
source(here::here(year, "code", "figure_functions.R"))
source(here::here(year, "code", "table_functions.R"))

# setup
afscassess::sp_switch(species)

# Run only once! (folder setup and accepted model)
# setup folder structure - only run this once
afscdata::setup_folders(year)
# setup of base folder (copies over the previously accepted model)
afscdata::accepted_model(2022, prev_model_name, year) # don't run if already setup
# setup .tpl files for data files
afscassess::setup_tpl(year)
# manually create or copy over the following folders:
# code
# .exe, .cpp, .tpl, and .ctl (if applicable) files


### data ----
# data query (folder 'raw')
afscdata::goa_dusk(year, off_yr= FALSE)
#goa_dusk(year, akfin_user, akfin_pwd, afsc_user, afsc_pwd, off_yr= FALSE)
# next year fix goa_dusk to include 'akfin_load_date'
# AKFIN_LOAD_DATE" AS "akfin_load_date"

# Need area specific GAP index:
akfin = connect()
afsc_species = c(30150, 30152)
table_temp = dplyr::tbl(akfin, dplyr::sql("afsc.race_biomassareaaigoa")) %>%
  dplyr::rename_with(tolower) %>%
  dplyr::filter(year <= 2023,
                species_code %in% c(afsc_species),
                survey == 'GOA') %>%
  dplyr::collect()
vroom::vroom_write(table_temp, here::here(year, "data", "raw", "goa_area_bts_biomass_data.csv"),
                   delim = ",")

# query TAC values from AKFIN (historical harvest specs)
# doesn't work
afscdata::q_specs(year = year, species = species, area = "GOA", db = akfin, save = TRUE, print_sql=TRUE)
# download area specific specs from akfin
harv_specs <- read.csv(here::here(year, "data", "user_input", "goa_dusk_harvest_specs.csv"))

# query fmp species incidental catch from rockfish targeted fisheries
q_fmp_target_dusk(year, area = "GOA", target = "k", db= akfin, save = TRUE)

# query non-target catch (non-FMP species catch estimate)
afscdata::q_nontarget(year= year, target="k", area="GOA", db=akfin, save=TRUE)

# query psc
# dusky specific because had to arrange years in order
#targeted species: 'p' = pollock-mid, 'b' = pollock-bottom, x' = rex, 'h' = shallow flats,
#   'k' = rockfish, 'w' = arrowtooth, 'c' = pcod, 'i' = halibut
q_psc_dusk(year= year, target = "k", area = "GOA", db = akfin, save = TRUE)

# query non-commercial catch
q_non_comm_dusk(year = year, area= "GOA", species = species, db= akfin, save = TRUE)

# split fractions (AKFIN answers)

# user_input files to either obtain new each full assessment or copy over
# Copy from previous assessment: reader_tester.csv, lbins.csv, fixed_catch.csv, goa_dusk_catch_1977_1990.csv
# User input new: VAST index, ABC and TAC (or append)

# assemble and format data files
# located in data > output

# fishery catch (note: this automates the in-year estimation)
## output/yld_ratio.csv has the expansion factor ($ratio) and the 3-year catch/TAC ratio ($yld)
## which are used for in-year and next-two-year catches, respectively

# output: fish_catch.csv and yld_rat.csv
suppressWarnings(afscassess::clean_catch(year = year,
                                         species = species,
                                         TAC = TAC))
# fishery discard rate
fish_discard_dusk(year = year)

# fishery age comp
# base case (currently used)
# output: fish_age_comp.csv
afscassess::fish_age_comp(year = year,
                          rec_age = rec_age,
                          plus_age = plus_age,
                          rmv_yrs = c(1987, 1989))

# Why are the numbers changing from this year compared to 2022?

# fishery length comp
# output: fish_length_comp.csv
afscassess::fish_length_comp(year = year,
                             rec_age = rec_age,
                             lenbins = len_bins, #lengths
                             rmv_yrs = c(1987, 1989, 1990, year)) # makes sure removes current assessment year

# trawl survey biomass (not used only for reference)
# output: goa_total_bts_biomass.csv
# Note: run for VAST comparison, but use VAST bts index for assessment model
vast_bts <- read.csv(here::here(year, "data", "user_input", "vast_lognormal.csv")) %>%
  dplyr::filter(Stratum == "Stratum_1") %>%
  dplyr::select(year, biomass = total_biomass, se) %>%
  dplyr::mutate(biomass = biomass/1000, # changing units kg --> mt
                se = se/1000,
                lci = biomass - 1.96 * se,
                uci = biomass + 1.96 * se) %>%
  dplyr::select(year, biomass, se, lci, uci) %>%
  vroom::vroom_write(here::here(year, "data", "output", "goa_total_bts_biomass_vast.csv"), ",")

#afscassess::bts_biomass(year = year,
#                        rmv_yrs = c(1984, 1987),
#                        file = "vast_lognormal_goa.csv",
#                        id = "vast")

# get from gap products
#rema_bts_biomass
afscassess::bts_biomass(year = year,
                        rmv_yrs = c(1984, 1987),
                        id= "rema")

gap_dat <- read.csv(here::here(year, "data", "output", "goa_total_bts_biomass_rema.csv")) %>%
  dplyr::mutate(strata = "goa",
                cv = se/biomass ) %>%
  dplyr::select(strata, year, biomass, cv)


input_db_rema <- prepare_rema_input( model_name = 'db_rema',
                                     biomass_dat = gap_dat,
                                     end_year = max(gap_dat$year),
                                     zeros = list(assumption = 'NA' ) )
mod_db_rema <- fit_rema(input_db_rema)
out_db_rema <- tidy_rema(mod_db_rema)

write.csv(out_db_rema$total_predicted_biomass, here::here(year, "data", "output", "goa_biom_rema.csv"))
# trawl survey age comp
# input: goa_bts_agecomp_data.csv and goa_bts_specimen_data.csv
#file.copy(from=here::here(year, 'data','raw','bts_specimen_data.csv'),
#          to = here::here(year, 'data','raw','goa_bts_age_specimen_data.csv'))
# output: goa_bts_age_comp.csv

afscassess::bts_age_comp(year = year,
                         area = "goa",
                         rec_age = rec_age,
                         plus_age = plus_age,
                         rmv_yrs = c(1984,1987))


# trawl survey length comps
# input: goa_bts_sizecomp_data.csv and goa_bts_length_data.csv
# output: goa_bts_sizecomp.csv

afscassess::bts_length_comp(year = year,
                            area = "goa",
                            bysex= NULL,
                            rmv_yrs = c(1984,1987),
                            #sa_index = 2,
                            lenbins = len_bins)

## manually rename a file; lookup funs want "bts"
#file.rename(from=here::here(year, 'data','output','goa_ts_length_comp.csv'),
#            to = here::here(year, 'data','output','goa_bts_length_comp.csv'))

## Requires running admb; Models are in data > model folder
# Make sure all the .exe, .tpl, .ctl files in all folders

# ageing error matrix
# User input: reader_tester.csv
# runs ageage model
# output: ae_model_csv, ae_mtx_100.csv, ae_SD.csv
afscassess::age_error(year = year,
                      reader_tester = "reader_tester.csv",
                      species = species,
                      rec_age = rec_age,
                      plus_age = plus_age)

# current size-age matrix
# requires running the ageing error matrix first
# input: ae_model.csv, goa_bts_specimen_data.csv, goa_bts_length_data.csv
# runs length_sd and vonb models
# output laa_stats.csv, lbar_params.csv, saa.csv
afscassess::size_at_age(year = year,
                        area = "goa",
                        rec_age = rec_age,
                        lenbins = len_bins)

# weight-at-age
# runs allometric and wvonb models
# copy .ctl file over too for wvonb model
# output: alpha_beta_lw.csv and waa.csv, waa_stats.csv, wal_stats.csv, Wbar_params.csv
afscassess::weight_at_age(year = year,
                          rec_age = rec_age,
                          area = "goa")

# note: must provide file for VAST or DB estimates are output
# requires running the ageing error matrix first
# input: ae_model.csv, goa_bts_specimen_data.csv
# If need a rema comparison model

# concatenate dat file, for now writing it to output folder in data
# model folders:
base_folder <- "m22.3a_base"
alt_model_folder <- "m22.5a"

concat_dat_dusk(year = year,
                species = species,
                area = "goa",
                folder = "m22.3a_base", # name folders based on model name
                dat_name = dat_name,
                rec_age = rec_age,
                plus_age = plus_age,
                spawn_mo = 3,
                bts_biom_id = "vast") # bts_biom_id = NULL (original or vast or rema)

# Need to output the parameter estimates from outside model into one file for SAFE

# returns years in datasets to check
data_years_fn(yr = year, area = "goa")

### Run models ----

# Copy alternative models over to main year's folder
#if(dir.exists(here::here(year, "alt_models")) == FALSE ) {dir.create(here::here(year, "alt_models")) }
if(dir.exists(here::here(year, "base_2024"))) {
  file.copy(from = list.files(here::here(year, "alt_models_september", "m1.base"), full.names= TRUE),
            to = here::here(year, "base_2024"))
}
if(dir.exists(here::here(year, "m22.5a"))) {
  file.copy(from = list.files(here::here(year, "alt_models_september", "m3.srvproj_v3"), full.names= TRUE),
            to = here::here(year, "m22.5a"))
}

# Delete all insides of the folders except:
# 1. folder structure
# 2. base.cpp, base.exe, base.tpl, goa_dr_yyyy.ctl
# 3. .ctl file: Line 2 --> change file name to match new .dat file (with current assessment year)
#               Line 5 --> change end year
#               Line 59 --> change yield_ratio to yld from new yield_ratio
# 4. make sure any new model code updates in admb are compiled again

base_folder = base_mod = "m22.3a_base"
alt_model_folder = alt_folder = alt_mod <- "m22.5a"

## Base (previously accepted model) with new data

# goa_dusk_yyyy.dat for each model
concat_dat_dusk(year = year,
                species = species,
                area = "goa",
                folder = "m22.3a_base", # name folders based on model name
                dat_name = dat_name,
                rec_age = rec_age,
                plus_age = plus_age,
                spawn_mo = 3,
                bts_biom_id = "vast") # bts_biom_id = NULL (original or vast or rema)

year_mod <- 2024

# Examine initial model runs
# Run base and alternative accepted model
# Converge? Check gradient and LL (if not converged use pin files)
# Pin files: save a converged model .par file as .pin file and rerun

# Model 1 base
# Change working directory
base_mod = "m22.3a_base"
setwd(here::here(year, "m22.3a_base")) ; getwd()

R2admb::run_admb('base', verbose= TRUE)

# Model 2 alternative
# Skipping just adding the lognormal survey error... PT and SSC accepted model 3

# Model 3 alternative
#lognormal survey error + corrected projection of the recruitment age starting year
srvproj_mod = "m22.5a" #"m3.srvproj"
setwd(here::here(year, "m22.5a")) ; getwd()

R2admb::run_admb('base', verbose= TRUE)

### Model evaluation and comparison (initial examination) ---
base_mod = "m22.3a_base"
alt_mod = srvproj_mod = "m22.5a"

# saving extra files for model comparison
report_admb_dusk(model_year = year_mod, year = year, model= base_mod)
report_admb_dusk(model_year = year_mod, year = year, model = alt_mod)

# run to process results with and without mcmc
process_results_dusk(year = year, folder = base_mod, admb_name = "base", dat_name = "goa_dusk",
                     mcmc = FALSE, rec_age = rec_age, plus_age = plus_age, len_bins = len_bins)

process_results_dusk(year = year, folder = alt_mod, admb_name = "base", dat_name = "goa_dusk",
                     mcmc = FALSE, rec_age = rec_age, plus_age = plus_age, len_bins = len_bins)

process_results_dusk(year = year, folder = "base", admb_name = "base", dat_name = "goa_dr_2022",
                     mcmc = FALSE, rec_age = rec_age, plus_age = plus_age, len_bins = len_bins)
# note that the nll are not comparable between the two models because the model
#      structure is different

### Figures for basic comparison and data fit (data composition fits) ---
#  (will rerun these again after mcmc to ensure final values are used)
models_comp = c("m22.5a", "m22.3a_base", "base")

plot_compare_survey_dusk(ayr= year, models = models_comp, legend_x = .15, legend_y = .85, save = TRUE, final = FALSE)
# still need to put the error bars from vast model in legend
plot_compare_biom_f_dusk(ayr = year, models = c("m22.5a", "m22.3a_base", "base"), save = TRUE, final= FALSE)

### Figures for Data Composition Fits ---
# Run data composition fit ensure data fits look ok
afscassess::plot_comps(year = year, folder = alt_mod, save = TRUE)

# For other plots see: https://github.com/BenWilliams-NOAA/afscassess/blob/main/R/plots.R

### MCMC ---

# Model 3 alternative
#lognormal survey error + corrected projection of the recruitment age starting year

# setup new folder inside model
alt_mod = "m22.5a" #"m3.srvproj"

if (!dir.exists(here::here(year, alt_mod, "mcmc"))){
  dir.create(here::here(year, alt_mod, "mcmc"), recursive=TRUE)
}

file.copy(c(paste0("base.tpl"),
            paste0("base.exe"),
            paste0(dat_name, "_", year, ".dat"),
            paste0(dat_name, ".ctl"),
            "MAT.dat"),
          here::here(year, alt_mod, "mcmc"),
          overwrite = TRUE)

# run mcmc
setwd(here::here(year, alt_mod, "mcmc")) ; getwd()

R2admb::compile_admb('base', verbose=FALSE ) # also might need to run 1st and rerun with pin file
R2admb::run_admb('base', verbose= TRUE)

# run prelim mcmc (takes only a few min)
mcmcruns <- 1e5
mcmcsave <- mcmcruns/ 50
# check gradients/ convergence

# full mcmc runs
# final model: run -mcmc 10,000,000 -mcsave 2000 (~6 hrs) Ben has a virtual machine --> ask IT for access to a virtual machine
mcmcruns <- 1e7
mcmcsave <- mcmcruns/ 2000 # --> should have been done mcmcruns/5000 to save 2000, now I have saved 5000

# mcmc runs = 1e6, mcsave 2000 (I think saves 500?) --> total time: 20 min (in admb command center)

R2admb::run_admb('base', verbose= TRUE, mcmc = TRUE,
                 mcmc.opts = R2admb::mcmc.control(mcmc = mcmcruns,
                                                  mcsave = mcmcsave))

st_time <- Sys.time()
system(paste0('base', '.exe',' -mcmc ',mcmcruns,' -mcsave ',mcmcsave))
system(paste0('base','.exe',' -mceval')) # very fast
end_time <- Sys.time()
tot_time <- end_time - st_time

# rerun process results
process_results_dusk(year = year, folder = paste0(alt_mod,"/mcmc"), admb_name = "base", dat_name = "goa_dusk",
                     rec_age= rec_age, plus_age= plus_age, mcmc= TRUE, mcsave = mcmcsave, n_mcmc= mcmcruns, len_bins = len_bins)

# process other results
mcmc_eval_summary_dusk(year = year, model= alt_mod)

### Figures for Data Composition Fits and Model Outputs ---
# Run these now to ensure data fits look ok
afscassess::plot_comps(year = year, folder = paste0(alt_mod,"/mcmc"), save = TRUE)
afscassess::plot_selex(year = year, folder = paste0(alt_mod,"/mcmc"), save = TRUE)
plot_params_dusk(year = year, folder = paste0(alt_mod,"/mcmc"), model_name = "base",
                             pars = c("q_srv1", "ABC", "nattymort", "log_mean_rec", "tot_biom", "F40","spawn_biom"), save = TRUE)
## fix params so bars go up higher
afscassess::plot_phase(year = year, folder = paste0(alt_mod,"/mcmc"), model_name = "base", save = TRUE)
afscassess::plot_rec(year = year, folder = paste0(alt_mod,"/mcmc"), rec_age = rec_age, save = TRUE)
afscassess::plot_rec_ssb(year= year, folder = paste0(alt_mod,"/mcmc"), rec_age = rec_age, save = TRUE)
afscassess::plot_swath(year = year, folder = paste0(alt_mod,"/mcmc"), save = TRUE)

# afscassess::base_plots():
# plot_catch plot_survey plot_selex plot_params plot_phase
# plot_rec plot_rec_ssb plot_swath plot_comps

### Retrospective Analysis ---
# for testing
#mcmcruns_ret <- 10000
#mcmcsave_ret <- mcmcruns_ret / 5

# for final run
mcmcruns_ret = 500000  # Could change these, but 500,000 is a manageable number to deal with
mcmcsave_ret = mcmcruns_ret / 250

run_retro_dusk(year = year, model = alt_mod, tpl_name = "base", dat_name = "goa_dusk",
               n_retro = 10, mcmc_on = TRUE, mcmc_num = mcmcruns_ret, mcsave_num = mcmcsave_ret)

# Plot retrospective plots for spawning biomass and total biomass (and save Mohn's rho)
# saves plots in the model folder > figs and mohn's rho csv in model folder > processed
st_time <- Sys.time()

plot_retro_dusk(year = year, folder = alt_mod, n_retro=10, save_indiv=TRUE, rhos = c("ssb", "totbiom"))

end_time <- Sys.time() ; tot_time <- end_time - st_time

### Projections ---

# Run projections in own "proj" folder (in model folder)
afscassess::proj_ak(year=year, species= "dusk", region=area,
                    rec_age=rec_age, folder=alt_mod, off_yr=FALSE)
# output folders: author_f, max_f,
# output in processed folder: mcscen_ssb, mscen_f, mscen_yld, exec_summ (for executive summary)

# see pop line 414 start of projection plotting
# Jim's projection model (not quite workig): https://github.com/afsc-assessments/spmR

### Best F (for SARA file) ---

# Use this year's model to recalculate the F from the previous year that would have
#     produced naa, waa, etc from terminal year of model (from GOA POP & 2022 Dusky)

projdat <- readLines(here::here(year, alt_mod, "proj.dat"))
ages_proj <- as.numeric(unlist(strsplit(projdat[grep("#_Natural_Mortality", projdat)], " ") )[-1] )
naa_proj <- as.numeric(unlist(strsplit(projdat[grep("#_Numbers_at_age_end_year", projdat)+1], " ") )[-1] )
waa_proj <- as.numeric(unlist(strsplit(projdat[grep("#_Wt_at_age_spawners", projdat)+1], " ") )[-1] )
saa_proj <- as.numeric(unlist(strsplit(projdat[grep("#_Selectivity_fishery_scaled_to_max_at_one", projdat)+1], " ") )[-1] )
natmort <- 0.07
rep_dat <- readLines(here::here(year, "base", "base.rep"))
exec_summ_last <- read.csv(here::here(year-1, "harvest_proj", "processed", "exec_summ.csv"))
ofl_last <- OFL_2023 # ASK BEN is this supposed to be the values from 2023 assessment (for 2024) or from 2022 assessment (for 2023)
f40_last <- as.numeric(exec_summ_last %>% filter(item == "FABC") %>% pull(y3))

afscassess::best_f(data = data.frame(age = ages_proj, naa = naa_proj, waa= waa_proj, saa= saa_proj),
                   m = natmort,
                   type = 1, ## one sex, one gear,
                   f_ratio= NULL, ## only one gear
                   last_ofl = ofl_last, ## last model's OFL for last year (2022 for 2022, https://meetings.npfmc.org/CommentReview/DownloadFile?p=b1e05a01-b24d-4eaa-9a99-06b118df5f57.pdf&fileName=GOA%20POP%20PRESENTATION.pdf)
                   last_f = f40_last)

### Run apportionment ---
#setup folder
# download/ query split EGOA fraction from gap_products/ AKFIN answers
egoa_split_frac <- read.csv(here::here(year, 'data', 'user_input', 'goa_dusk_egoa_biomass_split_fractions.csv')) %>%
  dplyr::select(year = Year, eyak_biom = Eastern.Biomass..mt., eyak_frac = Eastern.Fraction,
                wyak_biom = Western.Biomass..mt., wyak_fract = Western.Fraction)

biom_gap_dat <- read.csv(here::here(year,'data','raw','goa_area_bts_biomass_data.csv')) %>%
  dplyr::rename(area= regulatory_area_name ) %>%
  dplyr::group_by(year, area) %>%
  dplyr::summarise(biomass= sum(area_biomass, na.rm=T)/1000,
                   se = sqrt(sum(biomass_var, na.rm=T))/1000) %>%
  dplyr::mutate(area = tolower(stringr::str_remove(area, " GOA")))

stata_names <- data.frame('area' = as.factor(c("Total","Western","Central","Eastern")),
           'Stratum' = c("Stratum_1", "Stratum_2", "Stratum_3", "Stratum_4"))

biom_vast_dat <- read.csv(here::here(year, 'data','user_input','vast_lognormal_index.csv')) %>%
  dplyr::left_join( stata_names ) %>%
  dplyr::select(year = Time, area, biomass = Estimate, se = Std..Error.for.Estimate) %>%
  dplyr::mutate(area = tolower(area),
                biomass = biomass/1000,
                se = se/1000) %>% # units bio0mass = kg (need to convert to t later)
  dplyr::filter(year >= 1990,
                biomass > 0 ,
                !area == "total")

run_apport_dusk(year = year, model = alt_mod, biom_dat = biom_vast_dat,
                biom_name = "vast", frac_dat = egoa_split_frac)
run_apport_dusk(year = year, model = alt_mod, biom_dat = biom_gap_dat,
                biom_name = "gap", frac_dat = egoa_split_frac)
end_year_apport_dusk(year= year, model = alt_mod, biom_name = "vast",abc = NULL)

egoa_frac <- read.csv(here::here(year, alt_mod, "apport", "apport_prop_vast.csv")) %>% mutate(value = as.character(wyak))
egoa_split_plot <- ggplot(egoa_frac %>% filter(!is.na(wyak)), aes(x= year, y= wyak)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label =value), vjust = -0.2) +
  labs(y="Proportion given to WYAK")

ggplot2::ggsave(plot = egoa_split_plot, here::here(year, alt_mod, "figs", "egoa_split_fig.png"),
                width = 4, height = 3.5, units = "in", dpi = 200)

# figure for apportionment
apport_catch_plot(year = year, model = alt_mod, model_name= "vast")

### Tables for SAFE ---
# Other tables not already produced in tables.Rmd or above code

# formatted data comp table
data_comps_table(year = year, folder = alt_mod, comp_type = c("fac", "fsc", "sac", "ssc")[1] ) # fac table
data_comps_table(year = year, folder = alt_mod, comp_type = c("fac", "fsc", "sac", "ssc")[2] ) # fsc table
data_comps_table(year = year, folder = alt_mod, comp_type = c("fac", "fsc", "sac", "ssc")[3] ) # sac table
data_comps_table(year = year, folder = alt_mod, comp_type = c("fac", "fsc", "sac", "ssc")[4] ) # ssc table

fac_resid <- pearson_resid_plot(comp_type = c("fac", "fsc", "sac", "ssc")[1] ,
                                 label = 'Age', outlier=3, title_name ="Fishery")
fsc_resid <- pearson_resid_plot(comp_type = c("fac", "fsc", "sac", "ssc")[2] ,
                                label = 'Length', outlier=3, title_name ="Fishery")
sac_resid <- pearson_resid_plot(comp_type = c("fac", "fsc", "sac", "ssc")[3] ,
                                label = 'Age', outlier=3, title_name ="Survey")

resid_temp <- cowplot::align_plots(fac_resid,sac_resid, fsc_resid)
pearson_plot <- cowplot::ggdraw() +
  cowplot::draw_grob(resid_temp[[1]], x = 0, y = 0.5, width = 0.5, height = 0.5) +
  cowplot::draw_grob(resid_temp[[2]], x = 0, y = 0, width = 0.5, height = 0.5) +
  cowplot::draw_grob(resid_temp[[3]], x = .5, y = .3, width = 0.5, height = 0.7)

ggplot2::ggsave(plot = pearson_plot, here::here(year, alt_mod, "figs", "pearson.png"),
                width = 8, height = 8, units = "in", dpi = 200)

# likelihood, penalty, prior, param est table
model_compare_ll_penalties(year, base_pr_mod = "base", base_mod= "m22.3a_base", alt_mod = alt_mod)
# output each component separately; need to put together in single csv file (with a line separation) for SAFE
# save final combined csv as "likelihood.csv"

# key parameter estimates from mle and mcmc
mcmc_eval_summary_dusk(year= year, model= alt_mod)
# model > tables> t_pars_mle_mcmc.csv

# estimates numbers
numbers_mat_selex_dusk(year = year, model = alt_mod, model.name = "base", rec_age = rec_age, plus_age = plus_age)

# time series comparison of ssb, biom, catch rate, recruits, with past model
mod_compare_dusk(year = year, prev_mod = "base", pref_mod= alt_mod, admb_name = "base")

# times series of recruits and biom from mcmc
recr_B_tbl_dusk(year= year, model= alt_mod, model_name= "base", rec_age= rec_age)
# output: models>tables> t_recr_biom_mcmc.csv

# 50% maturity by age and length:
# look in base.std file: a50 = 10.3
#$L_\infty$= 48.260, $\kappa$ = 0.182, $t_0$ = 0.401
# Linf(1-exp(-(k(t-t0))))
48.260*(1-exp(-(0.182*(10.3-0.040))))
55.004*(1-exp(-(0.124*(10.3- (-0.072)))))

### Figures for SAFE ---
# Other figures for SAFE

# Plot bts survey catch of past most recent years
sp_dist_bts_plot(years = c(2019, 2021, 2023), save=TRUE)

# Plot of catch
plot_catch_dusk(year = year, folder= alt_mod, save=TRUE)

# Final survey plot with rema
models_comp = c("m22.5a", "m22.3a_base", "base")
plot_compare_survey_dusk(ayr= year, models = models_comp, legend_x = .15, legend_y = .85,
                         save = TRUE, final = TRUE, add_rema = "TRUE")


# Biomass, F time series model comparison

plot_compare_biom_F_dusk(year = year, pref_mod = alt_mod,  base_current = "m22.3a_base", base_pr_mod = "base", admb_name = "base")

# recruitment deviations

plot_recr_devs_dusk(year = year, pref_mod = alt_mod, save = TRUE)

# catch rate
plot_catch_rate(year = year, admb_name = "base", model = alt_mod)

# final ssb and tb plot
plot_final_biomass(year = year, model = alt_mod, save = TRUE )

# Apportionment + subarea ABC investigation

catch_wy <- read.csv(here::here(year, "data", "raw", "fish_catch_data.csv")) %>%
  dplyr::filter(year == 2024,
                fmp_subarea %in% c("WY")) %>%
  dplyr::mutate(month = month(week_end_date)) %>%
  dplyr::group_by(fmp_subarea, trip_target_name, fmp_gear, month) %>%
  dplyr::summarise(weight = sum(weight_posted, na.rm=T))

catch_wy_july <- read.csv(here::here(year, "data", "raw", "fish_catch_data.csv")) %>%
  dplyr::filter(year == 2024,
                fmp_subarea %in% c("WY")) %>%
  dplyr::mutate(month = month(week_end_date)) %>%
  dplyr::filter(month == 7) %>%
  dplyr::select(year, species_name, retained_or_discarded, trip_target_name, fmp_gear,
                fmp_subarea, week_end_date, weight = weight_posted, vessel_id)

# now we can directly pull the TAC values: AKR.V_CAS_TAC


####### Stopped here (and above to clean up and put as function)
# q_specs (harvest specs in https://github.com/afsc-assessments/afscdata/blob/main/R/queries.R)

## Params table??


# setup for next year
afscassess::setup_folders(year+1)
afscassess::accepted_model(2024, "m22.5a", year +1)

