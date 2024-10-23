# 2024 dusky rockfish assessment
# kristen.omori@noaa.gov
# based on ben.williams@noaa.gov

# load ----
#library(groundfishr)
#devtools::unload("groundfishr")
#devtools::install_github("afsc-assessments/afscdata", force = TRUE)
#devtools::install_github("BenWilliams-NOAA/afscassess", force = TRUE)

library(afscdata)
library(afscassess)
library(tidyverse)
library(ggplot2)
library(here)
library(keyring)
library(rema)
theme_set(afscassess::theme_report())

library(scico)
#ggplot2::theme_set(theme_report())
#library(funcr)

# globals ----
year = 2024
species <- "DUSK"
area = "goa"
afsc_species1 =  30150
afsc_species2 = 30152
norpac_species = 330
TAC <- c(3700, 3676, 5389)
fishery = "fsh"
#admb_home = "C:/ADMB-13.0" # I use this because I have multiple admb versions
#curr_mdl_fldr = "2020.1-2023"
#prev_mdl_fldr = "2020.1-2021"
prev_model_name = "m22.3a"
dat_name = "goa_pop"

rec_age = 4
plus_age = 25
length

'%nin%' <- Negate('%in%')

keyring::key_list()
afsc_user = keyring::key_list("akfin")$username
afsc_pwd = keyring::key_get("akfin", keyring::key_list("akfin")$username)
akfin_user = keyring::key_list("afsc")$username
akfin_pwd = keyring::key_get("afsc", keyring::key_list("afsc")$username)

# setup(year = year)
# folder setup
# run only once

afscassess::sp_switch(species)
# setup folder structure - only run this once
# afscdata::setup_folders(year)
accepted_model(2022, prev_model_name, year) # don't run if already setup

# data ----
# data query
goa_dusk(year, akfin_user, akfin_pwd, afsc_user, afsc_pwd)

groundfishr::clean_catch(year=year, species=species, TAC=TAC)
groundfishr::age_error(reader_tester = "reader_tester.csv", species, year, rec_age=rec_age, plus_age = plus_age)
groundfishr::fish_age_comp(year=year, rec_age=rec_age, plus_age=plus_age)
groundfishr::ts_age_comp(year=year, rec_age=rec_age, plus_age=plus_age, rmv_yrs = c(1984, 1987))
groundfishr::fish_length_comp(year=year, rec_age=rec_age, lenbins = 'lbins.csv')
groundfishr::ts_length_comp(year=year, lenbins = 'lbins.csv', bysex = FALSE,  rmv_yrs = c(1984, 1987))
groundfishr::size_at_age(year=year, rec_age=rec_age, lenbins="lbins.csv")
groundfishr::weight_at_age(year=year, rec_age=rec_age)

# note: must provide file for VAST or DB estimates are output

# Base model w/updated data & diff surveys
# design-based model
groundfishr::ts_biomass(year=year, rmv_yrs = c(1984, 1987))
# base model with design-based survey
groundfishr::concat_dat(year=year, species=species, model="db", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST default
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_default.csv', alt = "m22")
groundfishr::concat_dat(year=year, species=species, model="m22", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST lognormal
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_lognormal.csv', alt = "m22a")
groundfishr::concat_dat(year=year, species=species, model="m22a", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)


# increase length comp ----
# design-based model
ts_biomass(year=year, rmv_yrs = c(1984, 1987))
groundfishr::fish_length_comp(year=year, rec_age=rec_age, lenbins = 'lbins2.csv')
groundfishr::ts_length_comp(year=year, lenbins = 'lbins2.csv', bysex = FALSE,  rmv_yrs = c(1984, 1987))
groundfishr::size_at_age(year=year, rec_age=rec_age, lenbins="lbins2.csv")
groundfishr::concat_dat(year=year, species=species, model="db.1", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST default
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_default.csv', alt = "m22.1")
groundfishr::fish_length_comp(year=year, rec_age=rec_age, lenbins = 'lbins2.csv')
groundfishr::ts_length_comp(year=year, lenbins = 'lbins2.csv', bysex = FALSE,  rmv_yrs = c(1984, 1987))
groundfishr::size_at_age(year=year, rec_age=rec_age, lenbins="lbins2.csv")
groundfishr::concat_dat(year=year, species=species, model="m22.1", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST lognormal
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_lognormal.csv', alt = "m22.1a")
groundfishr::fish_length_comp(year=year, rec_age=rec_age, lenbins = 'lbins2.csv')
groundfishr::ts_length_comp(year=year, lenbins = 'lbins2.csv', bysex = FALSE,  rmv_yrs = c(1984, 1987))
groundfishr::size_at_age(year=year, rec_age=rec_age, lenbins="lbins2.csv")
groundfishr::concat_dat(year=year, species=species, model="m22.1a", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)


# increase age comp ----
plus_age = 30
# design-based model
ts_biomass(year=year, rmv_yrs = c(1984, 1987))
groundfishr::age_error(reader_tester = "reader_tester.csv", species, year, rec_age=rec_age, plus_age = plus_age)
groundfishr::fish_age_comp(year=year, rec_age=rec_age, plus_age=plus_age)
groundfishr::ts_age_comp(year=year, rec_age=rec_age, plus_age=plus_age, rmv_yrs = c(1984, 1987))
groundfishr::fish_length_comp(year=year, rec_age=rec_age, lenbins = 'lbins.csv')
groundfishr::ts_length_comp(year=year, lenbins = 'lbins.csv', bysex = FALSE,  rmv_yrs = c(1984, 1987))
groundfishr::size_at_age(year=year, rec_age=rec_age, lenbins="lbins.csv")
groundfishr::weight_at_age(year=year, rec_age=rec_age)
groundfishr::concat_dat(year=year, species=species, model="db.2", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST default
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_default.csv', alt = "m22.2")
groundfishr::concat_dat(year=year, species=species, model="m22.2", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST lognormal
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_lognormal.csv', alt = "m22.2a")
groundfishr::concat_dat(year=year, species=species, model="m22.2a", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)


# increase age & length comp ----
plus_age = 30
# design-based model
ts_biomass(year=year, rmv_yrs = c(1984, 1987))
groundfishr::fish_length_comp(year=year, rec_age=rec_age, lenbins = 'lbins2.csv')
groundfishr::ts_length_comp(year=year, lenbins = 'lbins2.csv', bysex = FALSE,  rmv_yrs = c(1984, 1987))
groundfishr::size_at_age(year=year, rec_age=rec_age, lenbins="lbins2.csv")
groundfishr::concat_dat(year=year, species=species, model="db.3", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST default
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_default.csv', alt = "m22.3")
groundfishr::concat_dat(year=year, species=species, model="m22.3", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# VAST lognormal
ts_biomass(year=year, rmv_yrs = c(1984, 1987), file = 'vast_lognormal.csv', alt = "m22.3a")
groundfishr::concat_dat(year=year, species=species, model="m22.3a", dat_name='goa_dr_2022', spawn_mo=3, rec_age=rec_age, plus_age=plus_age)

# results ----

process_results(year, model="db", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=25, mcmc = 100000, mcsave=100)
process_results(year, model="m22", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=25, mcmc = 100000, mcsave=100)
process_results(year, model="m22a", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=25, mcmc = 100000, mcsave=100)

# plus length
process_results(year, model="db.1", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=25, mcmc = 100000, mcsave=100)
process_results(year, model="m22.1", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=25, mcmc = 100000, mcsave=100)
process_results(year, model="m22.1a", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=25, mcmc = 100000, mcsave=100)

# plus age
process_results(year, model="db.2", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=30, mcmc = 100000, mcsave=100)
process_results(year, model="m22.2", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=30, mcmc = 100000, mcsave=100)
process_results(year, model="m22.2a", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=30, mcmc = 100000, mcsave=100)

# plus and length
process_results(year, model="db.3", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=30, mcmc = 100000, mcsave=100)
process_results(year, model="m22.3", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=30, mcmc = 100000, mcsave=100)
process_results(year, model="m22.3a", model_name="base", dat_name="goa_dr_2022",
                rec_age=4, plus_age=30, mcmc = 100000, mcsave=100)

plot_compare_survey(year,
                    models = c('2022, db', '2022, db.1', '2022, m22', '2022, m22.1', '2022, m22.1a', '2022, m22.2'))
plot_compare_biomass(year,
                     models = c('2022, db', '2022, db.1', '2022, m22', '2022, m22.1', '2022, m22.1a', '2022, m22.2'))

plot_compare_survey(year,
                    models = c('2022, db', '2022, m22', '2022, m22a'))
plot_compare_biomass(year,
                     models = c('2022, db', '2022, m22', '2022, m22a'))

plot_compare_survey(year,
                    models = c('2022, db.1', '2022, m22.1', '2022, m22.1a'))
plot_compare_biomass(year,
                     models = c('2022, db.1', '2022, m22.1', '2022, m22.1a'))

plot_compare_survey(year,
                    models = c('2022, db.2', '2022, m22.2', '2022, m22.2a'))
plot_compare_biomass(year,
                     models = c('2022, db.2', '2022, m22.2', '2022, m22.2a'))

plot_compare_survey(year,
                    models = c('2022, db.3', '2022, m22.3', '2022, m22.3a'))
plot_compare_biomass(year,
                     models = c('2022, db.3', '2022, m22.3', '2022, m22.3a'))

plot_compare_biomass(year,
                     models = c('2022, db.3', '2022, m22.1a', '2022, m22.2a', '2022, m22.3a'))

# retrospective

run_retro(year, model = "m22.3a", tpl_name = "base", n_retro = 10, mcmc = 100000, mcsave = 100)

groundfishr::plot_comps(year, model = "m22.3a")
rockfishr::base_plots(year, model, model_name='base', rec_age = 4)
rockfishr::plot_swath(year, model)
afscassess:plot_swath(year, model)
rockfishr::plot_retro_survey(year, model)

rockfishr::param_table(year, model, model_name = 'base')

STD = read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)
mceval = read.csv(here::here(year, model, "processed", "mceval.csv"))

params <- function(year, STD, mceval, param){

  # doh - use the same names!!!
  param2 = ifelse(param == "nattymort", "natmort", param)

  std = STD %>%
    dplyr::filter(name == !!param)

  if(param == "spawn_biom_proj") {
    param2 = paste0("spawn_biom_proj_", year + 1)
    std = STD %>%
      dplyr::filter(name == !! param) %>%
      dplyr::slice_head(1)
  }
  mceval %>%
    dplyr::select(param2) %>%
    dplyr::summarise(mean = mean(!!dplyr::sym(param2)),
                     median = median(!!dplyr::sym(param2)),
                     sd = sd(!!dplyr::sym(param2)),
                     lci = quantile(!!dplyr::sym(param2), 0.025),
                     uci = quantile(!!dplyr::sym(param2), 0.975)) %>%
    dplyr::bind_cols(std) %>%
    dplyr::mutate(name := !! param) %>%
    dplyr::select(name, value, mean, median, std.dev, sd, lci, uci)


}
names(mceval)


model = "m22.3a"
rockfishr::plot_biomass(year, model)


# data files
REP = readLines(here::here(year, "base", paste0("base", ".rep")))
proj = readLines(here::here(year, model, "proj.dat"), warn=FALSE)

# helper function
split_item <- function(name, one = NULL){

  if(is.null(one)){
    x = strsplit(proj[grep(name, proj)], " ")[[1]]
    as.numeric(x[2:length(x)])
  } else {
    x = strsplit(proj[grep(name, proj) + 1], " ")[[1]]
    as.numeric(x[2:length(x)])
  }
}

# pull data
age = split_item("#_Natural_Mortality")
naa = split_item("#_Numbers_at_age_end_year", 1)
waa = split_item("#_Wt_at_age_spawners", 1)
saa = split_item("#_Selectivity_fishery_scaled_to_max_at_one", 1)

m = as.numeric(strsplit(proj[grep("#_Natural_Mortality", proj) + 1], " ")[[1]][[2]])
ofl = as.numeric(strsplit(REP[grep(paste0("OFL for ", year-1), REP) + 1] [[2]]," "))
f40 = as.numeric(strsplit(REP[grep("F_40", REP) + 1] [[1]]," "))

dat = list(best = data.frame(age, naa, waa, saa), m = m, last_ofl = ofl, last_f = f40)

groundfishr::best_f(data=dat$best, m=dat$m, last_ofl=dat$last_ofl, last_f=dat$last_f)


