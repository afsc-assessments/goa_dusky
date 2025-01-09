## Table Functions
## October 2024
## K. Omori

# functions
# data_comps_table
# model_compare_ll_penalties
# mcmc_eval_summary_dusk
# numbers_mat_selex_dusk
# mod_compare_dusk
# recr_B_tbl_dusk
# apport_table

# Discard table


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
    write.csv(here::here(year, folder, "tables", paste0("t_", comp_type, ".csv" ) ), row.names = FALSE)

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

  write.csv(ll_store, here::here(year, "safe_dev", "LL_compare.csv"), row.names = FALSE)
  write.csv(prior_pen_store, here::here(year, "safe_dev", "prior_pen_compare.csv"), row.names = FALSE)
  write.csv(par_est_store, here::here(year, "safe_dev", "par_est_compare.csv"), row.names = FALSE)

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
mod_compare_dusk <- function(year, prev_mod = "base", pref_mod, admb_name = "base") {

  ayr <- year

  base_dat <- read.csv(here::here(year, "base", "processed", "bio_rec_f.csv")) %>%
    dplyr::select(Year = year, tot_biom_base = tot_biom, sp_biom_base = sp_biom, F_base = "F", recruits_base = recruits)
  pref_dat <- read.csv(here::here(year, pref_mod, "processed", "bio_rec_f.csv")) %>%
    dplyr::select(Year = year, tot_biom_curr = tot_biom, sp_biom_curr = sp_biom, F_curr = "F", recruits_curr = recruits)

  compare_dat <- dplyr::full_join(base_dat, pref_dat, by = join_by( Year )) %>%
    dplyr::select(Year, sp_biom_base, sp_biom_curr, tot_biom_base, tot_biom_curr,
                  F_base, F_curr, recruits_base, recruits_curr)

  write.csv(compare_dat, here::here(year, "compare_models", "compare.csv"), row.names = FALSE, na = "")

  # For CI if needed can look at the figures (plot_compare_biom_F_dusk) function ~line 292
  # also not sure why the previous table had the age 6 and the catch/ age 6??


} # table


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

  # new for mean rec, total b, spawning biom

 ko2 <-  read.csv(here::here(year, model, "mcmc", "processed", "mceval.csv")) %>%
    dplyr::select(dplyr::starts_with(c( "tot_biom", "spawn_biom", "log_rec_dev", "rec_proj"))) %>%
    dplyr::mutate(id = 1:dplyr::n())

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
                  name ="frate" ) %>%
    dplyr::select(-log_avg_F)

  log_mean_rec_vec <- read.csv(here::here(year, model, "mcmc", "processed", "mceval.csv")) %>%
    dplyr::select("log_mean_rec") %>%
    dplyr::mutate(id = 1:dplyr::n())

  # get mcmc data - clean it, calculate annual uci and lci
  read.csv(here::here(year, model, "mcmc", "processed", "mceval.csv")) %>%
    dplyr::select(dplyr::starts_with(c( "tot_biom", "spawn_biom", "log_rec_dev", "rec_proj", "log_mean_rec"))) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::mutate_if(is.character, dplyr::funs(as.numeric(gsub(",", "", .)))) %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::full_join(log_mean_rec_vec) %>%
    dplyr::mutate(value = ifelse(grepl("log_rec_dev", name), (value + log_mean_rec), value ),
                  value = ifelse(grepl("log", name), exp(value), value),
                  year =  as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  name = dplyr::case_when(stringr::str_detect(name, "tot_biom") ~ "tot",
                                          stringr::str_detect(name, "spawn_biom") ~ "sp",
                                          stringr::str_detect(name, "log_rec") ~ "rec",
                                          stringr::str_detect(name, "rec_proj") ~ "rec")) %>%
    dplyr::select(-log_mean_rec) %>%
    dplyr::bind_rows(full_f_tmp) %>%
    dplyr::filter(year >= 1977 & year <= mod_year + 2) %>%
    dplyr::group_by(year, name) %>%
    dplyr::summarise( "mean" = mean(value, na.rm=T),
                      "lci" = quantile(value, 0.025),
                      "uci" = quantile(value, 0.975) ) %>%
    tidyr::pivot_wider(names_from = name, values_from = c("mean", "lci", "uci"), names_glue = "{name}_{.value}") %>%
    dplyr::mutate(rec_mean = rec_mean*1000,
                  rec_lci = rec_lci*1000,
                  rec_uci = rec_uci*1000) %>%
    dplyr::select(Year = year, frate_mean, frate_lci, frate_uci, rec_mean, rec_lci, rec_uci,
                  tot_mean, tot_lci, tot_uci, sp_mean, sp_lci, sp_uci)%>%
    write.csv(here::here(year, model,"tables", "t_recr_biom_mcmc.csv"), row.names = FALSE)

}

## Apportionment tables

apport_table <- function(year = year, model = alt_mod, biom_name = "vast") {

  comma <- function(x) format(round(as.numeric(x), digits = 0), big.mark = ",")
  #perc_fn <- function(x) paste0(as.character(round(x*100, 1)), "%")
  perc_fn <- function(x) round(x*100, 1)
  # make apportionment tables
  rec_table <- read.csv(here::here(year, model, 'processed', 'exec_summ.csv'))
  ABC_yr1 <- rec_table %>%
    dplyr::filter(item == "ABC (t)") %>%
    dplyr::mutate(y3 = as.numeric(y3)) %>%
    pull(y3)
  ABC_yr2 <- rec_table %>%
    dplyr::filter(item == "ABC (t)") %>%
    dplyr::mutate(y4 = as.numeric(y4)) %>%
    pull(y4)
  OFL_totals <- data.frame(Year = c(year+1, year+2),
                           mgmt = "OFL",
                           Area = "Total",
                           value = c(as.numeric(rec_table[11, 4]), as.numeric(rec_table[11, 5]) ) )

  apport_out <- read.csv(here::here(year, model, "apport", paste0("apport_prop_", biom_name,".csv")) )

  abc_apport_goa_dat <- apport_out %>%
    dplyr::filter(year == max(year)) %>%
    dplyr::select(Year = year, Western = western, Central = central, Eastern = eastern) %>%
    tidyr::pivot_longer(cols = c(Western, Central, Eastern), names_to = "area", values_to = "apport") %>%
    dplyr::mutate(ABC_sub1 = ABC_yr1*apport,
                  diff_calc1 = ABC_yr1 - sum(ABC_sub1),
                  ABC_sub2 = ABC_yr2*apport,
                  diff_calc2 = ABC_yr2 - sum(ABC_sub2)) %>%
    # adding rounding errors to Central
    dplyr::mutate(ABC_sub1_corr = case_when(max(diff_calc1) > 0 ~ case_when(area == "Central" ~ ABC_sub1 + diff_calc1,
                                                                        area != "Central" ~ ABC_sub1),
                                         max(diff_calc1) == 0 ~ ABC_sub1)) %>%
    dplyr::mutate(ABC_sub2_corr = case_when(max(diff_calc2) > 0 ~ case_when(area == "Central" ~ ABC_sub2 + diff_calc2,
                                                                            area != "Central" ~ ABC_sub2),
                                            max(diff_calc2) == 0 ~ ABC_sub2)) %>%
    dplyr::select(Area = area, apport, ABC1 = ABC_sub1_corr, ABC2 = ABC_sub2_corr)

  # abc and ofl for goa
  abc_apport_sub <- abc_apport_goa_dat %>%
    dplyr::select(Area, value = apport) %>%
    dplyr::mutate(mgmt = "Area Apportionment" ) %>%
    dplyr::mutate(value = perc_fn(value)) %>%
    tidyr::pivot_wider(names_from = Area, values_from = value) %>%
    dplyr::mutate(Total = 100)

  abc_apport_goa_tmp <- abc_apport_goa_dat %>%
    dplyr::mutate(Year = year +1) %>%
    dplyr::bind_rows(abc_apport_goa_dat %>%
                       dplyr::mutate(Year = year+2) ) %>%
    dplyr::mutate(ABC = case_when(Year == (year+1) ~ ABC1,
                                  Year == (year+2) ~ ABC2)) %>%
    dplyr::mutate(value = ABC,
                  mgmt = "ABC") %>%
    dplyr::select(Year, mgmt, Area, value)

  abc_apport_goa_tmp2 <- abc_apport_goa_tmp %>%
    dplyr::bind_rows(OFL_totals %>%
                       dplyr::mutate(value = as.numeric(value)) ) %>%
    dplyr::bind_rows(data.frame(Year = c(year+1, year+2),
                                mgmt = "ABC",
                                Area = "Total",
                                value = c(ABC_yr1, ABC_yr2))) %>%
    #dplyr::mutate(value = comma(value)) %>%
     tidyr::pivot_wider(names_from = Area, values_from = value) %>%
    dplyr::arrange(Year, mgmt)

  abc_apport_goa <- bind_rows(abc_apport_sub, abc_apport_goa_tmp2) %>%
    dplyr::select(Year, mgmt, Western, Central, Eastern, Total)

  # abc for egoa with split
  apport_egoa_dat <- apport_out %>%
    dplyr::filter(year == max(year)) %>%
    dplyr::select( wyak, seo)

  abc_appport_egoa  <- abc_apport_goa_tmp %>%
    dplyr::filter(Area == "Eastern") %>%
    dplyr::select(Year, mgmt, value) %>%
    dplyr::bind_cols(apport_egoa_dat) %>%
    dplyr::mutate(wyak = round(wyak*round(value), 0),
                  seo = round(value)-wyak) %>%
    dplyr::select(Year, mgmt, WYAK = wyak, SEO = seo)


  write.csv(abc_apport_goa, here::here(year, model, 'apport', 'abc_apport.csv'), row.names = FALSE)
  write.csv(abc_appport_egoa, here::here(year, model, 'apport', 'abc_apport_egoa.csv'), row.names = FALSE)

}

#apport_dusk(year= year, model = alt_mod, biom_name = "vast",abc = NULL)
