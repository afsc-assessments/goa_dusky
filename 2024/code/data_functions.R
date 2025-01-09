# data functions (2024)
# concat_dat_dusk()
# data_years_fn()
# fish_discard_dusk()
# queries:
# q_psc_dusk()
# q_fmp_target_dusk()
# q_non_comm_dusk

concat_dat_dusk <- function(year, species, area = "goa", folder, dat_name="goa_dusk", rec_age, plus_age, spawn_mo = 3,
                       maturity = NULL, n_ageage = 1, n_sizeage = 1, n_fleets = 1, n_ts = NULL, n_lls = NULL,
                       retro = NULL,  ryear = NULL, bts_biom_id =NULL){

  # create directory
  if (!dir.exists(here::here(year, folder))){
    dir.create(here::here(year, folder), recursive=TRUE)
  }


  if(length(grep(paste0(area,"_lls"),
                 list.files(here::here(year, "data", "output")), value=TRUE)) > 0){
    llslc = read.csv(here::here(year, "data", "output", paste0(area, "_lls_length_comp.csv")))
    llsb = read.csv(here::here(year, "data", "output", paste0(area, "_lls_biomass.csv")))
  }

  if(!is.null(maturity)){
    mature = as.vector(read.csv(paste0(here::here(year, "data", "user_input", maturity))) %>%
                         dplyr::rename_all(tolower) %>%
                         dplyr::select(-age))
  }

  if(is.null(ryear)) { ryear = year }

  fishery = grep("fish", list.files(here::here(year, 'data', "output")), value=TRUE)
  survey = grep("bts_", list.files(here::here(year, 'data', "output")), value=TRUE)
  ll_survey = grep("lls_", list.files(here::here(year, 'data', "output")), value=TRUE)

  catch = read.csv(here::here(year, "data", "output", grep("catch", fishery, value=TRUE)))
  waa = read.csv(here::here(year, "data", "output", "waa.csv"))
  saa = read.csv(here::here(year, "data", "output", "saa.csv"))
  ae = read.csv(here::here(year, "data", "output", "ae_model.csv"))
  fishac = read.csv(here::here(year, "data", "output", grep("age", fishery, value=TRUE)))
  fishlc = read.csv(here::here(year, "data", "output", grep("length", fishery, value=TRUE)))
  tsac = read.csv(here::here(year, "data", "output", grep("age", survey, value=TRUE)))
  tslc = read.csv(here::here(year, "data", "output", grep("size", survey, value=TRUE)))      ############ CHANGED 'length' to size
  if(is.null(bts_biom_id)) {
    tsb = read.csv(here::here(year, "data", "output",
                              grep("total_bts_biomass", survey, value=TRUE)))
  } else {
    tsb = read.csv(here::here(year, "data", "output",
                              grep(paste0("goa_total_bts_biomass_", bts_biom_id, ".csv"), survey, value=TRUE)))
  }

  if(!is.null(retro)) {

    catch = catch %>%
      dplyr::filter(Year <= ryear)
    fishac = fishac %>%
      dplyr::filter(year < ryear)
    fishlc = fishlc %>%
      dplyr::filter(year < ryear)
    tsac = tsac %>%
      dplyr::filter(year < ryear)
    tslc = tslc %>%
      dplyr::filter(year < ryear)
    tsb = tsb %>%
      dplyr::filter(year < ryear)
  }

  if(length(ll_survey) > 0){
    llsrpw = read.csv(here::here(year, "data", "output", grep("biomass", ll_survey, value=TRUE)))
    llsslc = read.csv(here::here(year, "data", "output", grep("length", ll_survey, value=TRUE)))
    llsrpn = read.csv(here::here(year, "data", "output", grep("numbers", ll_survey, value=TRUE)))
  }

  # tsb %>%
  #   dplyr::select(-X)  %>%
  #   dplyr::mutate(lci = t - sd *1.96,
  #                 uci = t + sd *1.96) -> tsb
  # names(tsb) <- c("year", "biomass", "se", "lci", "uci")
  m_nages = nrow(ae)
  nages = length(rec_age:plus_age)

  # get length bin info
  lbin = as.numeric(gsub("[^0-9.]", "",  colnames(tslc)))
  lbin = lbin[!is.na(lbin)]
  nlenbins = length(lbin)

  if(is.null(n_ageage)){
    n_ageage = 1
  }

  if(is.null(n_sizeage)){
    n_sizeage = 1
  }

  sep = "# ========================================================="
  rwt = "# -"

  # header ----
  header = c(rwt,
             sep,
             paste0("#", area, " ", species, " Rockfish .dat file for ADMB optimization"),
             paste ("# New data provided on:", read.table(file = here::here(year, "data", "raw", "data_called.txt"),
                                                          sep = "\t")[2,1]),
             "# Notes:",
             "#   ~ Weight-at-age and length-age transition matrix automatically updated",
             "#   ~ Formatted to conduct automated retrospective analysis",
             "#   ~ Does not use most recent years fishery size data",
             "#   ~ Does not use fishery size data in years when ages are expected",
             sep,
             "#",
             "#")

  # model inputs ----

  if(is.null(maturity)){
    mipv <- c(sep,
              "# Model input parameters/vectors",
              sep,
              "# Start and end years, recruitment age, number of age and length bins",
              "# Model start year (styr):",
              as.character(min(catch$year)),
              "# Model end year (endyr): #!",
              as.character(ryear),
              "# Age at recruitment (rec_age): #-",
              as.character(rec_age),
              "# Number of ages in data (nages_D):",
              as.character(nages),
              "# Number of ages in model (nages_M):",
              as.character(m_nages),
              "# Number of length bins (nlenbins):",
              as.character(nlenbins),
              "# Number of age-age transition matrices (n_ageage_mat):",
              as.character(n_ageage),
              "# Number of size-age transition matrices (n_sizeage_mat):",
              as.character(n_sizeage),
              "# Length bin labels (len_bin_labels):",
              paste(lbin, collapse=" "),
              "# Spawn month (spawn_fract):",
              as.character(spawn_mo),
              "#",
              "#")

  } else {
    mipv <- c(sep,
              "# Model input parameters/vectors",
              sep,
              "# Start and end years, recruitment age, number of age and length bins",
              "# Model start year (styr):",
              as.character(min(catch$year)),
              "# Model end year (endyr): #!",
              as.character(ryear),
              "# Age at recruitment (rec_age): #-",
              as.character(rec_age),
              "# Number of ages in data (nages_D):",
              as.character(nages),
              "# Number of ages in model (nages_M):",
              as.character(m_nages),
              "# Number of length bins (nlenbins):",
              as.character(nlenbins),
              "# Number of age-age transition matrices (n_ageage_mat):",
              as.character(n_ageage),
              "# Number of size-age transition matrices (n_sizeage_mat):",
              as.character(n_sizeage),
              "# Length bin labels (len_bin_labels):",
              paste(lbin, collapse=" "),
              "# Spawn month (spawn_fract):",
              as.character(spawn_mo),
              "#",
              "#")

    mat = c(sep,
            "# Proportion mature at age (p_mature):",
            paste0("#! ",
                   paste(mature$mature, collapse = " ")),
            "#-",
            "",
            "")
  }

  waa = c(sep,
          "# Weight-at-age (wt):",
          paste(waa$x, collapse=" "),
          "#",
          "#")

  # fishery catch ----
  fishery_catch = c(sep,
                    "# Fishery catch (mt): obs_catch(styr,endyr)",
                    sep,
                    paste0("#! ", paste(min(catch$year):year, collapse=" ")),
                    paste(catch$catch, collapse=" "),
                    "#-",
                    "",
                    "")
  # cpue ----
  # not currently used for northern rockfish
  cpue = c(sep,
           "# CPUE Data",
           sep,
           "# Number of CPUE years",
           "0",
           "# CPUE observations (leave blank if 0)",
           "",
           "")

  # trawl biomass ----
  trawl_biomass = c(sep,
                    "# Trawl Survey Biomass",
                    sep,
                    "#! Number of trawl surveys: nyrs_srv1",
                    as.character(nrow(tsb)),
                    "#- Trawl survey years: yrs_srv1(1,nyrs_srv1) #!",
                    paste(tsb$year, collapse=" "),
                    "#- Observed trawl survey biomass (mt): obs_srv1_biom(1,nyrs_srv1) #!",
                    paste(tsb$biomass, collapse=" "),
                    "#- SE of observed trawl survey biomass: obs_srv1_se(1,nyrs_srv1) #!",
                    paste(tsb$se, collapse=" "),
                    "#- Lower CI, 1.96*SE #!",
                    paste(tsb$lci, collapse=" "),
                    "#- Upper CI, 1.96*SE #!",
                    paste(tsb$uci, collapse=" "),
                    "#-",
                    "",
                    "")
  # long line survey biomass ----

  if(exists("llsrpw")){
    ll_biomass = c(
      sep,
      "# Longline Survey Biomass",
      sep,
      "# Number of longline surveys: nyrs_srv2",
      as.character(nrow(llsb)),
      "# Longline survey years: yrs_srv2(1,nyrs_srv2)",
      paste(llsb$year, collapse=" "),
      "# Observed longline survey biomass (mt): obs_srv2_biom(1,nyrs_srv2)",
      paste(llsb$rpw, collapse=" "),
      "# SE of observed longline survey biomass: obs_srv2_se(1,nyrs_srv2)",
      paste(llsb$sd, collapse=" "),
      "# Lower CI, 1.96*SE",
      paste(llsb$lci, collapse=" "),
      "# Upper CI, 1.96*SE",
      paste(llsb$uci, collapse=" "),
      "",
      "")
  } else {
    ll_biomass = c(
      sep,
      "# Longline Survey Biomass",
      sep,
      "# Number of longline surveys: nyrs_srv2",
      "1",
      "# Longline survey years: yrs_srv2(1,nyrs_srv2)",
      "1999",
      "# Observed longline survey biomass (mt): obs_srv2_biom(1,nyrs_srv2)",
      "1000",
      "# SE of observed longline survey biomass: obs_srv2_se(1,nyrs_srv2)",
      "100",
      "# Lower CI, 1.96*SE",
      "10",
      "# Upper CI, 1.96*SE",
      "10000",
      "",
      "")
  }

  # fishery age comp ----
  fac <- c(
    sep,
    "# Fishery Age Composition",
    sep,
    "#! Number of years: nyrs_fish_age",
    as.character(nrow(fishac)),
    "#- Fishery age comp years: yrs_fish_age #!",
    paste(fishac$year, collapse=" "),
    "#- Number of samples: nsamples_fish_age(1,nyrs_fish_age) #!",
    paste(fishac$n_s, collapse=" "),
    "#- Number of hauls: nhauls_fish_age(1,nyrs_fish_age) #!",
    paste(fishac$n_h, collapse=" "),
    "#- Index for age-age error matrix #!",
    paste(fishac$AA_Index, collapse=" "),
    "#- Observed fishery age compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages) #!",
    collapse_row(dplyr::select(fishac, -year, -n_s, -n_h, -AA_Index)),
    "#-",
    "",
    "")

  # trawl survey age comp ----

  tsac <- c(sep,
            "# Trawl Survey Age Composition",
            sep,
            "#! Number of years: nyrs_srv1_age",
            as.character(nrow(tsac)),
            "#- Trawl Survey age comp years: yrs_srv1_age #!",
            paste(tsac$year, collapse=" "),
            "#- Number of samples: nsamples_srv1_age(1,nyrs_srv1_age) #!",
            paste(tsac$n_s, collapse=" "),
            "#- Number of hauls: nhauls_srv1_age(1,nyrs_srv1_age) #!",
            paste(tsac$n_h, collapse=" "),
            "#- Index for age-age error matrix #!",
            paste(tsac$AA_Index, collapse=" "),
            "#- Observed trawl survey age compositions (proportions at age): oac_srv1(1,nyrs_srv1_age,1,nages) #!",
            collapse_row(dplyr::select(tsac, -year, -n_s, -n_h, -AA_Index)),
            "#-",
            "",
            "")

  # fishery length comp ----
  flc <- c(
    sep,
    "# Fishery Size Composition",
    sep,
    "#! Number of years:",
    as.character(nrow(fishlc)),
    "#- Fishery size comp years: #!",
    paste(fishlc$year, collapse=" "),
    "#- Number of samples:  #!",
    paste(fishlc$n_s, collapse=" "),
    "#- Number of hauls:  #!",
    paste(fishlc$n_h, collapse=" "),
    "#- Index for size-age error matrix #!",
    paste(fishlc$SA_Index, collapse=" "),
    "#- Observed fishery size compositions (proportions at age)#!",
    collapse_row(dplyr::select(fishlc, -year, -n_s, -n_h, -SA_Index)),
    "#-",
    "",
    "")

  # trawl survey size comp ----
  tslc <- c(
    sep,
    "# Trawl Survey Size Composition",
    sep,
    "#! Number of years:",
    as.character(nrow(tslc)),
    "#- Survey Years: #!",
    paste(tslc$year, collapse=" "),
    "#- Number of samples:#!",
    paste(tslc$n_s, collapse=" "),
    "#- Number of hauls: #!",
    paste(tslc$n_h, collapse=" "),
    "#- Index for size-age error matrix #!",
    paste(tslc$SA_Index, collapse=" "),
    "#- Observed survey size compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages) #!",
    collapse_row(dplyr::select(tslc, -year, -n_s, -n_h, -SA_Index)),
    "#-",
    "",
    "")

  # longline survey size comp ----
  if(exists("llslc")){

    llsc <- c(sep,
              "# Longline Survey Size Composition",
              sep,
              "# Number of years: nyrs_srv2_size",
              as.character(nrow(llslc)),
              "# Longline Survey size comp years: yrs_srv1_size",
              paste(llslc$year, collapse=" "),
              "# Number of samples: nsamples_srv2_size(1,nyrs_srv2_size)",
              paste(llslc$n_s, collapse=" "),
              "# Number of hauls: nhauls_srv2_size(1,nyrs_srv2_size)",
              paste(llslc$n_h, collapse=" "),
              "# Index for size-age error matrix",
              paste(llslc$SA_Index, collapse=" "),
              "# Observed longline survey size compositions (proportions at length): osc_srv2(1,nyrs_srv2_size,1,nlenbins)",
              collapse_row(dplyr::select(llslc, -year, -n_s, -n_h, -SA_Index)),
              "",
              "")
  } else {
    llsc <- c(sep,
              "# Longline Survey Size Composition, NOT USED IN MODEL, include one year of fake data",
              sep,
              "# Number of years: nyrs_srv2_size",
              "1",
              "# Longline Survey size comp years: yrs_srv1_size",
              "1999",
              "# Number of samples: nsamples_srv2_size(1,nyrs_srv2_size)",
              "99",
              "# Number of hauls: nhauls_srv2_size(1,nyrs_srv2_size)",
              "99",
              "# Index for size-age error matrix",
              "1",
              "# Observed longline survey size compositions (proportions at length): osc_srv2(1,nyrs_srv2_size,1,nlenbins)",
              paste(seq(1/nlenbins, 1/nlenbins, length.out=nlenbins), collapse=" "),
              "",
              "")
  }

  # size-age transition matrix ----
  sizeage <- c(sep,
               "# Size-age transition matrix: proportion at size given age: ",
               sep,
               collapse_row(dplyr::select(saa, -age)),
               "#",
               "",
               "")

  # age error matrix ----
  aa <- c(sep,
          "# age error transition matrix: ",
          sep,
          collapse_row(ae),
          "#",
          "",
          "")

  # eof ----
  eof <- c(sep,
           "# end of file marker",
           sep,
           "42",
           "#!")

  # Compile DAT file for ADMB ----

  if(is.null(maturity)){
    dat <- c(header,
             mipv,
             waa,
             fishery_catch,
             cpue,
             trawl_biomass,
             ll_biomass,
             fac,
             tsac,
             flc,
             tslc,
             llsc,
             sizeage,
             aa,
             eof)
  } else {
    dat <- c(header,
             mipv,
             mat,
             waa,
             fishery_catch,
             cpue,
             trawl_biomass,
             ll_biomass,
             fac,
             tsac,
             flc,
             tslc,
             llsc,
             sizeage,
             aa,
             eof)
  }

  if(is.null(retro)){
    write.table(dat, file = here::here(year, folder, paste0(dat_name, "_", year ,".dat")) ,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  } else {
    write.table(dat, file = here::here(year, folder, "retro", "model", ryear, paste0(dat_name, ".dat")) ,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}

data_years_fn <- function(yr, area = "goa") {

  fishery = grep("fish", list.files(here::here(yr, 'data', "output")), value=TRUE)
  survey = grep("bts_", list.files(here::here(yr, 'data', "output")), value=TRUE)

  catch = read.csv(here::here(yr, "data", "output", grep("catch", fishery, value=TRUE)))
  fishac = read.csv(here::here(yr, "data", "output", grep("age", fishery, value=TRUE)))
  fishlc = read.csv(here::here(yr, "data", "output", grep("length", fishery, value=TRUE)))
  bts_biom = read.csv(here::here(yr, "data", "output", grep("biomass_vast", survey, value=TRUE)))
  tsac = read.csv(here::here(yr, "data", "output", grep("age", survey, value=TRUE)))

  yr_tab <- data.frame(data = c("fishery", "fish_age", "fish_len", "bts_biom", "bts_age"),
                       first_year= c(min(catch$year), min(fishac$year), min(fishlc$year), min(bts_biom$year), min(tsac$year)),
                       end_year=  c(max(catch$year), max(fishac$year), max(fishlc$year), max(bts_biom$year), max(tsac$year)) )
  yr_tab
}

# Calculating discard rate
fish_discard_dusk <- function(year) {
  fish_dat <- read.csv(here::here(year, "data", "raw", "fish_catch_data.csv")) %>%
    dplyr::filter(year >= 2000,
                  agency_species_code %in% cas_code )

  fish_gear_disc_dat <- fish_dat %>%
    dplyr::mutate(gear = ifelse(fmp_gear == "TRW", "Trawl", "Fixed")) %>%
    dplyr::group_by(year, retained_or_discarded, gear) %>%
    dplyr::summarise(catch = sum(weight_posted, na.rm=T)) %>%
    dplyr::ungroup()

  fish_gear <- fish_gear_disc_dat %>%
    dplyr::group_by(year, gear) %>%
    dplyr::summarise(catch = sum(catch, na.rm=T)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(perc_catch = (catch/ sum(catch, na.rm=T) )*100 ) %>%
    dplyr::select(-catch) %>%
    tidyr::pivot_wider(names_from = gear, names_prefix = "catch_", values_from = perc_catch)

  fish_gear_discard <- fish_gear_disc_dat %>%
    tidyr::pivot_wider(names_from = retained_or_discarded, values_from = catch) %>%
    dplyr::mutate(discard = (D/ (D+R)) * 100) %>%
    dplyr::select(year, gear, discard) %>%
    tidyr::pivot_wider(names_from = gear, names_prefix = "discard_", values_from = discard)

  tot_discard <- fish_gear_disc_dat %>%
    dplyr::group_by(year, retained_or_discarded) %>%
    dplyr::summarise(catch = sum(catch, na.rm= T)) %>%
    tidyr::pivot_wider(names_from = retained_or_discarded, values_from = catch) %>%
    dplyr::mutate(tot_discard = (D/ (D+R)) * 100) %>%
    dplyr::select(year, tot_discard)

  fish_gear_discard_final <- dplyr::full_join(fish_gear, fish_gear_discard, by = join_by(year)) %>%
    dplyr::full_join(tot_discard, by = join_by(year)) %>%
    dplyr::select(Year = year, catch_Trawl, discard_Trawl, catch_Fixed, discard_Fixed, tot_discard)

  write.csv(fish_gear_discard_final, here::here(year, "data", "output", "fish_discard_gear.csv"), row.names = FALSE)

}

q_psc_dusk <- function(year, target, area, db, save = TRUE) {
  # globals
  area = toupper(area)
  area = if(isTRUE(area == "GOA")){
    area = c("WG", "CG", "WY", "EY", "SE")
  } else if(isTRUE(area=="BSAI")){
    area = c("BS", "AI")
  } else {
    area
  }

  target = toupper(target)
  yr = year

  # call table
  dplyr::tbl(db, dplyr::sql("council.comprehensive_psc")) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::select(year,
                  species = species_group_name,
                  psc = pscnq_estimate,
                  fmp_subarea, trip_target_code) %>%
    dplyr::filter(year >= yr-4, year <= yr,
                  fmp_subarea %in% area,
                  trip_target_code %in% target) %>%
    dplyr::group_by(year, species) %>%
    dplyr::summarise(psc = round(sum(psc, na.rm = T),3)) %>%
    dplyr::collect() %>%
    dplyr::mutate(year = as.integer(year)) %>%
    dplyr::arrange(year, species) %>%
    tidytable::pivot_wider(names_from = year, values_from = psc) -> psc

  if(isTRUE(save)){
    vroom::vroom_write(psc, here::here(year, "data", "output", "psc_catch.csv"),
                       delim = ",")
    message("PSC table written to data/output folder.")
  } else {
    psc
  }
}

q_fmp_target_dusk <- function(year, area= "GOA", target = "K", db, save = TRUE) {

  area = toupper(area)
  area = if(isTRUE(area == "GOA")){
    area = c("WG", "CG", "WY", "EY", "SE")
  } else if(isTRUE(area=="BSAI")){
    area = c("BS", "AI")
  } else {
    area
  }

  target = toupper(target)
  yr = year

  table_tmp <- dplyr::tbl(db, dplyr::sql("council.comprehensive_blend_ca")) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::filter(year >= yr-4, year <= yr,
                  fmp_subarea %in% area,
                  trip_target_code %in% target) %>%
    dplyr::mutate(conf = 1) %>%
    dplyr::group_by(year, species_group_name) %>%
    dplyr::summarise(weight = round(sum(weight_posted, na.rm = T),0),
                     conf = sum(conf)) %>%
    dplyr::collect() %>%
    dplyr::mutate(weight = ifelse(conf <=3 , "conf", weight)) %>%
    dplyr::select(year, species_group_name, weight) %>%
    dplyr::arrange(year, species_group_name) %>%
    tidytable::pivot_wider(names_from = year, values_from = weight)

  if(isTRUE(save)) {
    vroom::vroom_write(table_tmp, here::here(year, "data", "output", "fmp_catch.csv"),
                       delim = ",")
    message("FMP species group table written to data/output folder.")
  } else {
    table_tmp
  }

}

q_non_comm_dusk <- function(year, area= "GOA", species = species, db= akfin, save = TRUE) {

  area_name = toupper(area)
  area_sub = if(isTRUE(area == "GOA")){
    area_sub = c("WG", "CG", "WY", "EY", "SE")
  } else if(isTRUE(area=="BSAI")){
    area_sub = c("BS", "AI")
  } else {
    area
  }

  species_grp = toupper(species)

  yr = year
    non_fsh_catch_query <- paste0("select   *
                from     AKR.V_NONCOMMERCIAL_FISHERY_CATCH
                where    species_group_code = ('DUSK')
                         collection_year >= 2010" )

  query_non_fsh_catch_dat <- sql_run(akfin, non_fsh_catch_query) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(fmp_area_code %in% area_name )

  if(isTRUE(save)) {
    vroom::vroom_write(query_non_fsh_catch_dat, here::here(year, "data", "user_input", "non_comm_catch.csv"),
                       delim = ",")
    message("Non-commercial catch in user_input folder.")
  } else {
    table_tmp
  }

}
