## Table Functions
## October 2024
## K. Omori

# Age/ length comps with sample sizes
data_comps_table <- function(year, folder, comp_type ) {

  options(scipen = 999)
  if(comp_type == "fac") {
    comp_dat = read.csv(here::here(year, "data", "output", "fish_age_comp.csv"))
    eff_n_comp <- read.csv(here::here(year, folder, "processed","fac.csv"))
  }
  if(comp_type == "fsc") {
    comp_dat = read.csv(here::here(year, "data", "output", "fish_length_comp.csv"))
    eff_n_comp <- read.csv(here::here(year, folder, "processed","fsc.csv"))
  }
  if(comp_type == "sac") {
    comp_dat = read.csv(here::here(year, "data", "output", "goa_bts_age_comp.csv"))
    eff_n_comp <- read.csv(here::here(year, folder, "processed","sac.csv"))
  }
  if(comp_type == "ssc") {
    comp_dat = read.csv(here::here(year, "data", "output", "goa_bts_sizecomp.csv"))
    eff_n_comp <- read.csv(here::here(year, folder, "processed","ssc.csv"))
  }

  eff_n <- eff_n_comp %>%
    dplyr::select(year, n_ess = eff_N) %>% # for the effective sample size
    dplyr::distinct() %>%
    dplyr::mutate(n_ess = round(n_ess)) %>%
    tidyr::pivot_wider(names_from = year, values_from = n_ess )

  comp_dat %>%
    dplyr::select(n_s, n_h) %>%
    t(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") -> samps

  comp_dat %>%
    dplyr::select(-n_s, -n_h, -contains("_Index") ) %>%
    tidyr::pivot_longer(-c(year)) %>%
    tidyr::pivot_wider(names_from = year, values_from = value, names_prefix = "y") %>%
    as.data.frame() %>%
    dplyr::mutate_if(is.numeric, round, digits = 3) %>%
    dplyr::mutate(name = gsub("X", "", name),
                  name = ifelse(dplyr::row_number() == dplyr::n(), paste0(name, "+"), name )) %>%
    dplyr::rename_all(~stringr::str_replace(., "y", "")) -> comp

  names(samps) <- names(comp)
  if(comp_type %nin% c("ssc") ) {
    samps <- bind_rows(samps, eff_n %>% dplyr::mutate(name = "n_ess")) %>%
      dplyr::mutate(name = ifelse(name == "n_s", "n_iss", name))
  }

  dplyr::bind_rows(comp, samps) %>%
    write.csv(here::here(year, model, "tables", paste0("t_", comp_type, ".csv" ) ), row.names = FALSE)

} # end comp_table

# Model comparison for LL, priors/ penalties, param est
model_compare_ll_penalties <- function(year, base_pr_mod = "base", base_mod= "m22.3a_base", alt_mod = alt_mod) {

  folder_names <- c(base_pr_mod, base_mod, alt_mod)

  for(i in 1:length(folder_names)) {
    ll_val_tmp <- read.csv(here::here(year, folder_names[i], "processed", "ll_values.csv")) %>%
      dplyr::mutate(LL = round(LL, 2))
    colnames(ll_val_tmp) <- c("Likelihood", folder_names[i])

    prior_pen_tmp <- read.csv(here::here(year, folder_names[i], "processed", "prior_penalty_values.csv")) %>%
      dplyr::mutate(value = round(value, 2))
    colnames(prior_pen_tmp) <- c("Penalties/Priors", folder_names[i])

    par_est_tmp <- read.csv(here::here(year, folder_names[i], "processed", "parameter_est_values.csv"))
    colnames(par_est_tmp) <- c("Parameter Estimates", folder_names[i])

    if(i == 1) {
      ll_store <- ll_val_tmp
      prior_pen_store <- prior_pen_tmp
      par_est_store <- par_est_tmp
    }
    if(i>1) {
      ll_store <- dplyr::full_join(ll_store, ll_val_tmp)
      prior_pen_store <- dplyr::full_join(prior_pen_store, prior_pen_tmp)
      par_est_store <- dplyr::full_join(par_est_store, par_est_tmp)
    }
  } # end for loop

  write.csv(ll_store, here::here(year, "safe_input", "LL_compare.csv"), row.names = FALSE)
  write.csv(prior_pen_store, here::here(year, "safe_input", "prior_pen_compare.csv"), row.names = FALSE)
  write.csv(par_est_store, here::here(year, "safe_input", "par_est_compare.csv"), row.names = FALSE)

}

## MCMC (mceval) model results summarized
mcmc_eval_summary_dusk <- function(year, model){

  # get mcmc data - clean it, calculate annual uci and lci
  # next year include a50 (in the mcmc.csv and maybe deltas?)


  mcmc_sum <- read.csv(here::here(year, model, "mcmc" ,"processed", "mceval.csv")) %>%
    dplyr::select(sigr, q_srv1, log_mean_rec, F40, SSB_proj = spawn_biom_proj,
                tot_biom_proj = paste0("tot_biom_proj_", year), ABC) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(-id, names_to = "Parameter", values_to = "value") %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(mean_mcmc = mean(value),
                     median_mcmc = median(value),
                     sd_mcmc = sd(value),
                     lci_mcmc = quantile(value, 0.025),
                     uci_mcmc = quantile(value, 0.975))
  mle_sum <- read.delim(here::here(year, model, "base.std"), sep="", header = TRUE) %>%
    dplyr::filter(name %in% c("sigr", "q_srv1", "log_mean_rec", "F40", "spawn_biom_proj", "ABC", "tot_biom_proj")) %>%
    dplyr::filter(!(name == "spawn_biom_proj" & index > 343),
                  !(name == "tot_biom_proj" & index > 358))  %>%
    dplyr::mutate(Parameter = ifelse(name == "spawn_biom_proj", "SSB_proj", name )) %>% # note hard coded these values need to automate
    dplyr::select(Parameter, mean_mle = value, sd_mle = std.dev)

  mle_mcmc <- dplyr::full_join(mle_sum, mcmc_sum, by = join_by(Parameter) ) %>%
    dplyr::mutate(Parameter = factor(Parameter,
                     levels= c("sigr", "q_srv1", "log_mean_rec", "F40", "SSB_proj", "tot_biom_proj", "ABC"))) %>%
    write.csv(here::here(year, model, "tables", "t_pars_mle_mcmc.csv"), row.names = FALSE)
}

# Estimated numbers w/ selectivity, weight, mature
numbers_mat_selex_dusk <- function(year, model, model.name = "base", rec_age, plus_age) {

  age_dusk <- rec_age:plus_age

  selex <- read.csv(here::here(year, model, "processed", "selex.csv"))
  REP <- readLines(here::here(year, model, paste0(model_name, ".rep")))

  suppressWarnings(data.frame(age = unlist(strsplit(REP[grep("Age", REP)[1]]," "))) %>%
                     dplyr::mutate(age = as.numeric(age)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(age)) -> ages_long

  yrs_length <- length(1977:year)
  num_age <- data.frame(age_tmp = ages_long,
                        Num = as.numeric(unlist(strsplit(REP[grep("Numbers", REP) + 48]," "))[-(1:2)])) %>%
    dplyr::mutate(age = ifelse(age_tmp >30, 30, age_tmp)) %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(Num = round(sum(Num)*1000) )
  wt_age <- data.frame(age = age_dusk,
                       Weight = round(as.numeric(unlist(strsplit(REP[grep("Weight", REP)]," "))[-(1:2)])[1:length(age_dusk)] ) )

  Age_num_selex <- selex %>%
    dplyr::mutate(Mature = round(maturity*100),
                  Fishery = round(fish, 2),
                  Survey = round(srv1, 2)) %>%
    dplyr::filter(age <= 30) %>%
    dplyr::full_join(num_age, by = join_by(age)) %>%
    dplyr::full_join(wt_age,  by = join_by(age) ) %>%
    dplyr::select(Age = age, Abundance = Num, Mature, Weight, Fishery, Survey) %>%
    dplyr::mutate(Age = ifelse(Age == 30, "30+", Age))

  write.csv(Age_num_selex, here::here(year, model,"tables", "t_age_num_selex.csv"), row.names = FALSE)

}

# Time series of SB, 6+ biomass, catch rate, age 4+ recruits with last year's model
modout_compare_dusk <- function(year, prev_mod = "base", pref_mod) {

  base_dat <- read.csv(here::here(year, "base", "processed", "bio_rec_f.csv")) %>%
    dplyr::mutate(model = "Previous")
  pref_dat <- read.csv(here::here(year, pref_mod, "mcmc", "processed", "bio_rec_f.csv")) %>%
    dplyr::mutate(model = "Current")


}

recr_B_tbl_dusk <- function(year, model, model_name= "base", rec_age){

  # hard to filter year with year so change the name
  mod_year = year

  if (!dir.exists(here::here(year, model, "mcmc", "processed"))){
    stop("must run 'process_results' before creating tables")
  }

  # read in data
  REP <- readLines(here::here(year, model, "mcmc", paste0(model_name, ".rep")))
  STD <- read.delim(here::here(year, model, "mcmc", paste0(model_name, ".std")), sep="", header = TRUE)

  suppressWarnings(data.frame(year = unlist(strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     dplyr::mutate(year = as.numeric(year)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(year)) -> yrs

  f_bio_rec = read.csv(here::here(year, model, "processed", "bio_rec_f.csv")) %>%
    #dplyr::select(-F) %>%
    dplyr::mutate(recruits = round(recruits * 1000))

  f_bio_rec %>%
    dplyr::filter(year %in% (1977 + rec_age):(mod_year - rec_age)) %>%
    dplyr::summarise(recruits = round(mean(recruits))) -> pred_rec

  data.frame(year = mod_year + 1:2,
             tot_biom = STD$value[which(STD$name=="tot_biom_proj")][1:2],
             sp_biom = STD$value[which(STD$name=="spawn_biom_proj")][1:2],
             recruits = pred_rec$recruits) -> std_data

  values = dplyr::bind_rows(f_bio_rec, std_data) %>%
    dplyr::mutate("Frate" = round(F, 3),
                  recruits = round(recruits, 0),
                  sp_biom = round(sp_biom, 0),
                  tot_biom = round(tot_biom, 0))

  # mcmc data for fully selected F
  year_vec <- 1977:mod_year

  F_dat <- read.csv(here::here(year, model, "mcmc", "processed", "mcmc.csv")) %>%
    dplyr::select( "log_avg_F" , dplyr::starts_with("log_F_devs")) #%>%
  colnames(F_dat) <- c("log_avg_F", paste0("log_F_devs_", year_vec) )

  F_avg <- F_dat %>% select(log_avg_F) %>%
    dplyr::mutate( id = 1:dplyr::n())
 # (add these to the each F_dev_year, then exp() to get the F_devs to Frate)

  full_f_tmp <- F_dat %>%
    dplyr::select(!log_avg_F) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::mutate_if(is.character, dplyr::funs(as.numeric(gsub(",", "", .)))) %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::full_join(F_avg, by = join_by(id)) %>%
    dplyr::mutate(value = value + log_avg_F) %>%
    dplyr::mutate(value = exp(value),
                  year =  as.numeric(stringr::str_extract(name, "[[:digit:]]+")) ,
                  name ="Frate" ) %>%
    dplyr::select(-log_avg_F)

  # get mcmc data - clean it, calculate annual uci and lci
  read.csv(here::here(year, model, "mcmc", "processed", "mceval.csv")) %>%
    dplyr::select(dplyr::starts_with(c( "tot_biom", "spawn_biom", "log_rec_dev", "rec_proj"))) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::mutate_if(is.character, dplyr::funs(as.numeric(gsub(",", "", .)))) %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::mutate(value = ifelse(grepl("log", name), exp(value), value),
                  year =  as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  name = dplyr::case_when(stringr::str_detect(name, "tot_biom") ~ "tot_biom",
                                          stringr::str_detect(name, "spawn_biom") ~ "sp_biom",
                                          stringr::str_detect(name, "log_rec") ~ "recruits",
                                          stringr::str_detect(name, "rec_proj") ~ "recruits")) %>%
    dplyr::bind_rows(full_f_tmp) %>%
    group_by(year, name) %>%
    dplyr::summarise(lci = quantile(value, 0.025),
                     uci = quantile(value, 0.975)) %>%
    dplyr::left_join(values, .) %>%
    dplyr::filter(year >= 1977 & year <= mod_year + rec_age) %>%
    dplyr::mutate(tot_lci = ifelse(name == 'tot_biom', lci, NA),
                  tot_uci = ifelse(name == 'tot_biom', uci, NA),
                  sp_lci = ifelse(name == 'sp_biom', lci, NA),
                  sp_uci = ifelse(name == 'sp_biom', uci, NA),
                  rec_lci = ifelse(name == 'recruits', lci * 1000, NA),
                  rec_uci = ifelse(name == 'recruits', uci * 1000, NA),
                  frate_lci = ifelse(name == 'Frate', lci, NA),
                  frate_uci = ifelse(name == 'Frate', uci, NA) ) %>%
    group_by(year) %>%
    dplyr::summarise(Frate = mean(Frate, na.rm=T),
                     Frate_lci = mean(frate_lci, na.rm=T),
                     Frate_uci = mean(frate_uci, na.rm=T),
                     recruits = mean(recruits, na.rm = T),
                     rec_lci = mean(rec_lci, na.rm = T),
                     rec_uci = mean(rec_uci, na.rm = T),
                     tot_biom = mean(tot_biom, na.rm = T),
                     tot_lci = mean(tot_lci, na.rm = T),
                     tot_uci = mean(tot_uci, na.rm = T),
                     sp_biom = mean(sp_biom, na.rm = T),
                     sp_lci = mean(sp_lci, na.rm = T),
                     sp_uci = mean(sp_uci, na.rm = T)) %>%
    write.csv(here::here(year, model,"tables", "t_recr_biom_mcmc.csv"), row.names = FALSE)

}

## Apportionment tables

apport_table <- function() {
  ########## BELOW from POP
  # make apportionment tables
  rec_table <- vroom::vroom(here::here(year, 'mgmt', model, 'processed', 'exec_summ.csv'))

  # abc
  apport_out$proportion_biomass_by_strata %>%
    tidytable::filter(year == max(year)) %>%
    tidytable::select(-model_name) %>%
    tidytable::rename(cgoa = 'CENTRAL GOA',
                      egoa = 'EASTERN GOA',
                      wgoa = 'WESTERN GOA') %>%
    tidytable::pivot_longer(cols = c(cgoa, egoa, wgoa),
                            names_to = 'region',
                            values_to = 'apport') %>%
    tidytable::mutate(apport = round(apport, digits = 3),
                      diff = 1 - sum(apport)) %>%
    tidytable::mutate(apport_corr = case_when(max(diff) > 0 ~ case_when(region == 'cgoa' ~ apport + diff,
                                                                        region != 'cgoa' ~ apport),
                                              max(diff) == 0 ~ apport)) %>%  # if rounding error happens, add to cgoa
    tidytable::select(year, region, apport_corr) %>%
    tidytable::mutate(y1 = round(rec_table$abc[1] * apport_corr, digits = 0),
                      y2 = round(rec_table$abc[2] * apport_corr, digits = 0),
                      diff_y1 = rec_table$abc[1] - sum(y1),
                      diff_y2 = rec_table$abc[2] - sum(y2)) %>%
    tidytable::mutate(y1_corr = case_when(max(diff_y1) > 0 ~ case_when(region == 'cgoa' ~ y1 + diff_y1,
                                                                       region != 'cgoa' ~ y1),
                                          max(diff_y1) == 0 ~ y1),
                      y2_corr = case_when(max(diff_y2) > 0 ~ case_when(region == 'cgoa' ~ y2 + diff_y2,
                                                                       region != 'cgoa' ~ y2),
                                          max(diff_y2) == 0 ~ y2)) %>%  # if rounding error happens, add to cgoa
    tidytable::select(region, apport_corr, y1_corr, y2_corr) %>%
    tidytable::rename(apport = 'apport_corr',
                      y1 = 'y1_corr',
                      y2 = 'y2_corr') -> abc_apport1

  abc_apport1 %>%
    tidytable::pivot_longer(-1) %>%
    tidytable::pivot_wider(names_from = region, values_from = value) %>%
    tidytable::mutate(goa = c(1, rec_table$abc[1], rec_table$abc[2])) -> abc_apport

  abc_apport1 %>%
    tidytable::filter(region == 'egoa') %>%
    tidytable::select(-apport) %>%
    tidytable::pivot_longer(cols = c(y1, y2),
                            names_to = 'year',
                            values_to = 'abc') %>%
    tidytable::mutate(wyak = round(wyak_p$wyak * abc, digits = 0),
                      eyak_se = round((1 - wyak_p$wyak) * abc, digits = 0),
                      diff = abc - (wyak + eyak_se),
                      wyak_corr = wyak + diff) %>% # if difference in rounding, add (or take out) from wyak
    tidytable::select(year, abc, wyak_corr, eyak_se) %>%
    tidytable::rename(wyak = 'wyak_corr') %>%
    tidytable::mutate(wyak_p = wyak_p$wyak) -> abc_apport_wyak

  # ofl
  abc_apport1 %>%
    tidytable::select(region, apport) %>%
    tidytable::mutate(y1 = round(rec_table$ofl[1] * apport, digits = 0),
                      y2 = round(rec_table$ofl[2] * apport, digits = 0),
                      diff_y1 = rec_table$ofl[1] - sum(y1),
                      diff_y2 = rec_table$ofl[2] - sum(y2)) %>%
    tidytable::mutate(y1_corr = case_when(max(diff_y1) > 0 ~ case_when(region == 'cgoa' ~ y1 + diff_y1,
                                                                       region != 'cgoa' ~ y1),
                                          max(diff_y1) == 0 ~ y1),
                      y2_corr = case_when(max(diff_y2) > 0 ~ case_when(region == 'cgoa' ~ y2 + diff_y2,
                                                                       region != 'cgoa' ~ y2),
                                          max(diff_y2) == 0 ~ y2)) %>%  # if rounding error happens, add to cgoa
    tidytable::select(region, apport, y1_corr, y2_corr) %>%
    tidytable::rename(y1 = 'y1_corr',
                      y2 = 'y2_corr') -> ofl_apport1

  ofl_apport1 %>%
    tidytable::filter(region == 'egoa') %>%
    tidytable::select(-apport) %>%
    tidytable::pivot_longer(cols = c(y1, y2),
                            names_to = 'year',
                            values_to = 'ofl') %>%
    tidytable::mutate(wyak = round(wyak_p$wyak * ofl, digits = 0),
                      eyak_se = round((1 - wyak_p$wyak) * ofl, digits = 0),
                      diff = ofl - (wyak + eyak_se),
                      wyak_corr = wyak + diff) %>% # if difference in rounding, add (or take out) from wyak
    tidytable::select(year, ofl, wyak_corr, eyak_se) %>%
    tidytable::rename(wyak = 'wyak_corr') %>%
    tidytable::select(-ofl) %>%
    tidytable::bind_cols(ofl_apport1 %>%
                           tidytable::filter(region != 'egoa') %>%
                           tidytable::select(-apport) %>%
                           tidytable::pivot_longer(cols = c(y1, y2),
                                                   names_to = 'year',
                                                   values_to = 'ofl') %>%
                           tidytable::pivot_wider(names_from = region, values_from = ofl) %>%
                           tidytable::select(-year)) %>%
    tidytable::mutate(ofl_wcgoa_wyak = wyak + cgoa + wgoa,
                      ofl_eyak_se = eyak_se,
                      ofl = c(rec_table$ofl[1], rec_table$ofl[2])) %>%
    tidytable::select(year, ofl_wcgoa_wyak, ofl_eyak_se, ofl) -> ofl_apport

  write.csv(abc_apport, here::here(year, 'mgmt', model, 'processed', 'abc_apport.csv'))
  write.csv(abc_apport_wyak, here::here(year, 'mgmt', model, 'processed', 'abc_apport_wyak.csv'))
  write.csv(ofl_apport, here::here(year, 'mgmt', model, 'processed', 'ofl_apport.csv'))
  ################################################################################ END FROM POP
}


