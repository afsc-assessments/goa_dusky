# functions for GOA dusky analysis
# process_results_dusk
# process_retro_dusk
# run_retro_dusk
# report_admb_dusk


process_results_dusk <- function(year, folder, admb_name = "base", dat_name,
                            rec_age, plus_age, len_bins,
                            mcmc= FALSE, n_mcmc = 100000, mcsave= 2000,
                            proj = FALSE, retro= FALSE, retro_mcmc= FALSE){

  # for mcmc runs: add extra 'mcmc folder' --> model > mcmc
  #  change the model to: 'model/mcmc
  # setup
  if (!dir.exists(here::here(year, folder, "processed"))){
    dir.create(here::here(year, folder, "processed"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, folder, "figs"))){
    dir.create(here::here(year, folder, "figs"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, folder, "tables"))){
    dir.create(here::here(year, folder, "tables"), recursive=TRUE)
  }

  # helper functions
  rep_item <- function(name){
    t <- strsplit(REP[grep(name, REP)]," ")
    t <- subset(t[[1]], t[[1]]!="")
    if(t[[1]][1] == "TWL"){
      as.numeric(t[3:length(t)])
    } else {
      as.numeric(t[2:length(t)])
    }
  }

  # For retrospective
  get_colnames <- function(ret_yr){
    eval_names = c("sigr", "q_srv1", "q_srv2", "F40", "natmort",
                   "ABC", "obj_fun",
                   paste0("tot_biom_", yrs[which(yrs <= ret_yr)]),
                   paste0("log_rec_dev_", seq(styr_rec, yrs[length(yrs[which(yrs <= ret_yr)])])),
                   paste0("spawn_biom_", yrs[which(yrs <= ret_yr)]),
                   "log_mean_rec",
                   paste0("spawn_biom_proj_", max(yrs[which(yrs <= ret_yr)]) + 1:15),
                   paste0("pred_catch_proj_", max(yrs[which(yrs <= ret_yr)]) + 1:15),
                   paste0("rec_proj_", max(yrs[which(yrs <= ret_yr)]) + 1:10),
                   paste0("tot_biom_proj_", max(yrs[which(yrs <= ret_yr)]) +1:2),
                   paste0("srv1_biom_",yrs_srv1_biom$srv_yr[which(yrs_srv1_biom$srv_yr <= ret_yr)]))
    eval_names
  }

  # not sure if it works
  get_retro_res <- function(index_yr){
    ret <- read.delim(list.files(paste0(model_dir, "/retro/results"), pattern="*.prj", full.names = TRUE)[index_yr], sep="", header = FALSE)
    ret_yr <- as.numeric(strsplit(strsplit(list.files(paste0(model_dir, "/retro/results"), pattern="*.prj", full.names = TRUE)[index_yr], split = c("mcmc_"))[[1]][2], split = ".prj")[[1]])
    ret %>%
      tidytable::rename(!!!setNames(colnames(ret), get_colnames(ret_yr))) %>%
      tidytable::select(contains("spawn_biom")) %>%
      tidytable::select(contains(as.character(yrs[1]:yrs[length(yrs[which(yrs <= ret_yr)])]))) -> ret_res
    ret_res = ret_res[(0.1 * mcmcruns_ret / mcmcsave_ret + 1):nrow(ret_res),]
    ret_res
  }



  # read in rep and ctl files
  REP <- readLines(here::here(year, folder, paste0(admb_name, ".rep")))
  CTL <- readLines(here::here(year, folder, paste0(dat_name, ".ctl")))
  STD <- read.delim(here::here(year, folder, paste0(admb_name, ".std")), sep="", header = TRUE)
  DAT <- readLines(here::here(year, folder, paste0(dat_name, "_", year, ".dat")))
  PAR <- readLines(here::here(year, folder, paste0(admb_name, ".par")))

  if(mcmc == TRUE) {
    PSV <- file(here::here(year, folder, paste0(admb_name, ".psv")), "rb")
    mceval <- read.delim(here::here(year, folder, "evalout.prj"), sep="", header=FALSE)
  }

  # clean rep file
  suppressWarnings(data.frame(year = unlist(strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     dplyr::mutate(year = as.numeric(year)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(year)) -> yrs

  suppressWarnings(data.frame(age = unlist(strsplit(REP[grep("Age", REP)[1]]," "))) %>%
                     dplyr::mutate(age = as.numeric(age)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(age)) -> ages

  #styr_rec <- yrs[1] - length(ages) + -2 #rec_age #
  styr_rec <- yrs[1] -length(ages) + 2 # "length(ages) --> length(rec_age:plus_grp) - 2 (for plus group and first age)
  # not estimating recruit devs in the first year for the plus group and initial recruitment (initial recruitment already has recruit_dev)

  suppressWarnings(as.data.frame(cbind(yrs = yrs, ages = ages, styr_rec = styr_rec)) %>%
                     dplyr::mutate(ages = replace(ages, duplicated(ages), NA),
                                   styr_rec = replace(styr_rec, duplicated(styr_rec), NA))) %>%
    write.csv(here::here(year, folder, "processed", "ages_yrs.csv"), row.names = FALSE)

  # For MCMC = TRUE

  if(isTRUE(mcmc)) {
    # MCMC parameters ----

    npar = readBin(PSV, what = integer(), n=1)
    mcmcs = readBin(PSV, what = numeric(), n = (npar * n_mcmc / mcsave))
    close(PSV)
    mcmc_params = matrix(mcmcs, byrow=TRUE, ncol=npar)
    mcmc_params = mcmc_params[501:nrow(mcmc_params),]
    colnames(mcmc_params) = STD$name[1:ncol(mcmc_params)]
    write.csv(mcmc_params, here::here(year, folder, "processed", "mcmc.csv"), row.names = FALSE)

    # mceval phase output ----

    #Curry's Change (adding burn in)
    mceval = mceval[501:nrow(mceval),]

    # not sure why the mceval is only 500 (do I need to change the eval something with number)

    dim(mceval) # has 222 columns
    # below has 220

    #1-8: Through objective function value

    colnames(mceval) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
                         "ABC", "obj_fun",
                         paste0("tot_biom_", yrs),
                         paste0("log_rec_dev_", seq(styr_rec, yrs[length(yrs)])),
                         paste0("spawn_biom_", yrs),
                         "log_mean_rec",
                         paste0("spawn_biom_proj_", max(yrs) + 1:15),
                         paste0("pred_catch_proj_", max(yrs) + 1:15),
                         paste0("rec_proj_", max(yrs) + 1:10),
                         paste0("tot_biom_proj_", max(yrs)))
    write.csv(mceval, here::here(year, folder, "processed", "mceval.csv"), row.names = FALSE)
  }


  # catch data ----

  pred = strsplit(REP[grep("Pred_Catch", REP)], " ")
  r1 = which(pred[[1]] == "Pred_Catch")
  r2 = which(pred[[1]] == "Pred_catch_later")
  r3 = which(pred[[1]] == "")
  pred = as.numeric(pred[[1]][-c(r1, r2, r3)])

  obs = strsplit(REP[grep("Obs_Catch", REP)], " ")
  r1 = which(obs[[1]] == "Obs_Catch")
  r2 = which(obs[[1]] == "Obs_Catch_Later")
  r3 = which(obs[[1]] == "")
  obs = as.numeric(obs[[1]][-c(r1, r2, r3)])

  data.frame(obs = obs, pred = pred) %>%
    write.csv(here::here(year, folder, "processed", "catch.csv"))


  # survey data ----
  syr = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][2]
  syr = strsplit(syr," ")
  syr = subset(syr[[1]], syr[[1]]!="")
  syr = as.numeric(syr[2:length(syr)])

  obs = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][4]
  obs = strsplit(obs," ")
  obs = subset(obs[[1]], obs[[1]]!="")
  obs = as.numeric(obs[2:length(obs)])

  se = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][6]
  se = strsplit(se," ")
  se = subset(se[[1]], se[[1]]!="")
  se = as.numeric(se[2:length(se)])

  pred = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][3]
  pred = strsplit(pred," ")
  pred = subset(pred[[1]], pred[[1]]!="")
  pred = as.numeric(pred[2:length(pred)])

  data.frame(year = syr, biomass = obs, se = se) %>%
    dplyr::mutate(lci = biomass - 1.96 * se,
                  uci = biomass + 1.96 * se ) %>%
    dplyr::bind_cols(pred = pred) %>%
    write.csv(here::here(year, folder, "processed", "survey.csv"), row.names = FALSE)


  # recruitment ----

  N = REP[grep("Numbers", REP):(grep("Obs_P_fish_age", REP)-2)]
  t = NA
  for(i in 1:length(yrs)){
    ts = as.numeric(strsplit(N[i+1]," ")[[1]][3])
    t = c(t, ts)}
  pred_rec = t[!is.na(t)]

  # biomass & F & recruits ----
  data.frame(year = yrs,
             tot_biom = rep_item("Tot_biom"),
             sp_biom = rep_item("SpBiom"),
             F = rep_item("Fully_selected_F"),
             recruits = pred_rec) %>%
    write.csv(here::here(year, folder, "processed", "bio_rec_f.csv"), row.names = FALSE)


  # selectivity ----
  data.frame(age = ages,
             fish = rep_item("Fishery_Selectivity"),
             srv1 = rep_item("TWL Survey_Selectivity"),
             maturity = rep_item("Maturity")) %>%
    write.csv(here::here(year, folder, "processed", "selex.csv"), row.names = FALSE)

  # yield ratio B40 & B35----

  data.frame(B40 = STD$value[which(STD$name=="B40")],
             B35 = as.numeric(REP[(grep("B_35",REP)+1):(grep("F_40",REP)[1]-1)]),
             yld_rat = as.numeric(unlist(strsplit(CTL[grep("yieldratio", CTL)], "\t"))[1])) %>%
    write.csv(here::here(year, folder, "processed", "b35_b40_yld.csv"), row.names = FALSE)

  # size comps ----

  #! this will need a switch for multiple surveys

  obs_f_a = REP[grep("Obs_P_fish_age",REP):(grep("Pred_P_fish_age",REP)-2)]
  pred_f_a = REP[grep("Pred_P_fish_age",REP):(grep("Obs_P_fish_size",REP)-2)]

  obs_f_l = REP[grep("Obs_P_fish_size",REP):(grep("Pred_P_fish_size",REP)-2)]
  pred_f_l = REP[grep("Pred_P_fish_size",REP):(grep("Obs_P_srv1_age",REP)-2)]

  obs_s_a = REP[grep("Obs_P_srv1_age",REP):(grep("Pred_P_srv1_age",REP)-2)]
  pred_s_a = REP[grep("Pred_P_srv1_age",REP):(grep("Obs_P_srv1_size",REP)-2)]

  obs_s_l = REP[grep("Obs_P_srv1_size",REP):(grep("Pred_P_srv1_size",REP)-2)]

  purrit_dusk(obs = obs_f_a, pred = pred_f_a, rec_age = rec_age, plus_age = plus_age, comp = "age", lenbins = len_bins) %>%
    write.csv(here::here(year, folder, "processed", "fac.csv"))

  purrit_dusk(obs= obs_f_l, pred= pred_f_l, rec_age= rec_age, plus_age= plus_age, comp = "length", lenbins = len_bins) %>%
    write.csv(here::here(year, folder, "processed", "fsc.csv"))

  purrit_dusk(obs= obs_s_a, pred= pred_s_a, rec_age= rec_age, plus_age= plus_age, comp = "age", lenbins = len_bins) %>%
    write.csv(here::here(year, folder, "processed", "sac.csv"))

  purrit_dusk(obs= obs_s_l, pred = NULL, rec_age= rec_age, plus_age= plus_age, comp = "length", lenbins = len_bins) %>%
    write.csv(here::here(year, folder, "processed", "ssc.csv"))

  # Likelihoods saves (KO added)

  # Likelihoods
  tot_ll <-  data.frame(name = "total_LL",
                        LL=as.numeric(unlist(strsplit(REP[grep("Total likelihood",REP)+1],split=" "))))
  data_ll <-  data.frame(name = "data_LL" ,
                         LL=as.numeric(unlist(strsplit(REP[grep("Data likelihood",REP)+1],split=" "))))

  wts_ll_tmp <- data.frame(tmp= REP[(grep("Wts_and_Likelihoods", REP)+1):(grep("SigmaR:", REP)-1)] )
  ll_df_tmp <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "_Likelihood")) %>%
    tidyr::separate(tmp, into= c("wt", "LL", "name1", "name2"), sep = "\\s") %>%
    tidyr::unite("name", name1, name2, na.rm=T) %>%
    dplyr::mutate(LL = as.numeric(LL)) %>%
    dplyr::filter(!LL==0)

  if(dim(ll_df_tmp)[1] == 6 ) {
    mat_ll_tmp <- data.frame(wt = "1", LL = 65.00, name = "Maturity_Likelihood")
    ll_df_tmp <- bind_rows(ll_df_tmp, mat_ll_tmp)

  }
  #ll_names <- c("Catch", "Trawl survey", "Fishery ages", "Survey ages", "Fishery lengths", "Recruitment devs")
  #ll_df <- bind_cols(ll_df, ll_names = ll_names)

  # not sure where maturity ll is


  # Create likelihood large table with names and values from model name

  # after rewriting will need to include the maturity LL in there
  ll_names <- c("Catch", "Trawl survey", "Fishery ages", "Survey ages", "Fishery lengths",
                "Recruitment devs", "Maturity", "Data LL", "Total LL")
  ll_df <- bind_rows(ll_df_tmp, data_ll, tot_ll) %>% select(-wt)
  ll_df <- bind_cols(ll_df, Likelihood = ll_names) %>%
    dplyr::select(Likelihood, LL)

  write.csv(ll_df, here::here(year, folder, "processed", "ll_values.csv"), row.names = FALSE)

  # Penalties/ Priors (KO added)

  priors_tmp <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "priors") ) %>%
    tidyr::separate(tmp, into= c( "wt" ,"LL", "name1", "name2", "name3"), sep = "\\s") %>%
    tidyr::unite("name", name1, name2, name3, sep= " ", na.rm=TRUE) %>%
    dplyr::select(name, value = LL) %>%
    dplyr::filter(value > 0)

  penalty_tmp <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "_Penalty")) %>%
    tidyr::separate(tmp, into= c("wt", "LL", "name1"), sep = "\\s") %>%
    dplyr::filter(!LL==0) %>%
    dplyr::select(name = name1, value= LL)

  obj_fn <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "obj_fun")) %>%
    tidyr::separate(tmp, into= c("wt",  "value" , "name"), sep = "\\s") %>%
    dplyr::select(name, value)

  prior_penalty_df <- dplyr::bind_rows(priors_tmp, penalty_tmp, obj_fn)

  write.csv(prior_penalty_df, here::here(year, folder, "processed", "prior_penalty_values.csv"), row.names = FALSE)

  # Parameter estimates (KO added)
 # NEXT YEAR need to add the sd with the par est
  # single values
  num_par <- as.numeric(unlist(strsplit(REP[grep("Number parameters estimated",REP)+1],split=" ")))
  B100 <- round(as.numeric(unlist(strsplit(REP[grep("B_zero",REP)+1],split=" "))),digits=0) # B0 = B100
  B40 <- round(as.numeric(unlist(strsplit(REP[grep("B_40",REP)+1],split=" "))),digits=2) #Get B40 from report file
  F40 <- round(as.numeric(unlist(strsplit(REP[grep("F_40",REP)+1],split=" "))), digits=3) [1]
  F.ABC <- round(as.numeric(unlist(strsplit(REP[grep("F_ABC",REP)+1],split=" "))), digits=3)[1]
  ABC <- round(as.numeric(unlist(strsplit(REP[grep("ABC for",REP)+1],split=" "))), digits=0)[3]   # first 2 are the F_ABC values, this grabs ABC for first projection year
  q <- as.numeric(unlist(strsplit(REP[grep("q_trawl",REP)+1],split=" ")))[1]
  avg_log_rec <- round(as.numeric(unlist(strsplit(REP[grep("log_mean_rec",REP)+1],split=" "))), digits=3)[1]
  total_B <- round(as.numeric(unlist(strsplit(REP[grep("tot_biom",REP)+1],split=" "))), digits=0)[1]
  spawn_B <- round(as.numeric(unlist(strsplit(REP[grep("spawn_biom",REP)+1],split=" "))), digits=0)[1] # SSB for model_year+1
  sigmar <- as.numeric(unlist(strsplit(REP[grep("SigmaR: ",REP)],split=" ")))[2]
  a50 <- PAR[(grep("a50:", PAR)+1)]

  param_names <- c("# parameters", "sigmaR" , "q", "avg rec", "a50", "F40", "Total Biomass", "SSB", "B100", "B40", "ABC")
  par_est_values <- c(num_par, sigmar, q, exp(avg_log_rec), a50, F40, total_B, spawn_B, B100, B40, ABC)
  par_est_df <- bind_cols(Parameter = param_names, Estimates = par_est_values)

  write.csv(par_est_df, here::here(year, folder, "processed", "parameter_est_values.csv"), row.names = FALSE)

}

print("process_results_dusk")

# copied from afscaccess -> utils (but changed the length bin read in)
purrit_dusk <- function(obs, pred = NULL, rec_age, plus_age, comp = "length", lenbins = len_bins){

  if(is.null(lenbins)){
    lenbins = read.csv(here::here(year, "data", "user_input", "len_bin.csv"))$len_bins
  } else {
    lenbins = len_bins
  }

  obs_tmp = stringr::str_split(obs, " ")
  obs_tmp_samp_size <- unlist(obs_tmp)
  eff_N <- as.numeric(obs_tmp_samp_size[grep("eff_N",obs_tmp_samp_size)+1])
  N_tmp <- as.numeric(obs_tmp_samp_size[grep("eff_N",obs_tmp_samp_size)+3]) # "N" searches for all "N's"

  purrr::map_if(obs_tmp[-1], is.character, as.numeric) %>%
    purrr::map_df(., ~as.data.frame(t(.))) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) -> obs_tmp


  if(comp == "age" & !is.null(pred)){
    obs_tmp %>%
      dplyr::select(1:(plus_age - rec_age + 2)) -> obs_tmp

    pred_tmp = stringr::str_split(pred, " ")
    purrr::map_if(pred_tmp[-1], is.character, as.numeric) %>%
      purrr::map_df(., ~as.data.frame(t(.))) %>%
      dplyr::select_if(~sum(!is.na(.)) > 0) %>%
      dplyr::select(1:(plus_age - rec_age + 2)) -> pred_tmp

    names(pred_tmp) <- names(obs_tmp) <- c("year", rec_age:plus_age)

    obs_tmp %>%
      tidyr::pivot_longer(!year, names_to = "age") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::bind_rows(pred_tmp %>%
                         tidyr::pivot_longer(!year, names_to= "age") %>%
                         dplyr::mutate(groups = "pred")) %>%
      dplyr::mutate(age = as.integer(age),
                    Age = factor(age),
                    Year = factor(year)) -> dat_tmp

  } else if(comp == "length" & !is.null(pred)){

    obs_tmp %>%
      dplyr::select(1:(length(lenbins) + 1)) -> obs_tmp

    pred_tmp = stringr::str_split(pred, " ")
    purrr::map_if(pred_tmp[-1], is.character, as.numeric) %>%
      purrr::map_df(., ~as.data.frame(t(.))) %>%
      dplyr::select_if(~sum(!is.na(.)) > 0) %>%
      dplyr::select(1:(length(lenbins) + 1)) -> pred_tmp

    names(pred_tmp) <- names(obs_tmp) <- c("year", lenbins)

    obs_tmp %>%
      tidyr::pivot_longer(!year, names_to = "length") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::bind_rows(pred_tmp %>%
                         tidyr::pivot_longer(!year, names_to= "length") %>%
                         dplyr::mutate(groups = "pred")) %>%
      dplyr::mutate(length = as.integer(length),
                    Length = factor(length),
                    Year = factor(year)) -> dat_tmp

  } else if(comp == "age" & is.null(pred)){
    obs_tmp %>%
      dplyr::select(1:(plus_age - rec_age + 2)) -> obs_tmp

    names(obs_tmp) <- c("year", rec_age:plus_age)

    obs_tmp %>%
      tidyr::pivot_longer(!year, names_to= "age") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::mutate(age = as.integer(age),
                    Age = factor(age),
                    Year = factor(year)) -> dat_tmp

  } else if(comp == "length" & is.null(pred)){

    obs_tmp %>%
      dplyr::select(1:(length(lenbins) + 1)) -> obs_tmp

    names(obs_tmp) <- c("year", lenbins)

    obs_tmp %>%
      tidyr::pivot_longer(!year, names_to= "length") %>%
      dplyr::mutate(groups = "obs") %>%
      dplyr::mutate(length = as.integer(length),
                    Length = factor(length),
                    Year = factor(year)) -> dat_tmp
  }


  N_Neff <- data.frame(Year = unique(dat_tmp$Year) , eff_N = eff_N, N_tmp = N_tmp )
  dat <- dat_tmp %>% dplyr::full_join(N_Neff, by = join_by(Year))
}



process_retro_dusk <- function(year, model, model_name = alt, dat_name,
                          rec_age, plus_age, mcmc = 100000, mcsave = 100, folder = '2013'){

  # setup
  if (!dir.exists(here::here(year, model, "retro", "model", folder, "processed"))){
    dir.create(here::here(year, model, "retro", "model", folder, "processed"), recursive=TRUE)
  }


  # helper functions
  rep_item <- function(name){
    t <- strsplit(REP[grep(name, REP)]," ")
    t <- subset(t[[1]], t[[1]]!="")
    if(t[[1]][1] == "TWL"){
      as.numeric(t[3:length(t)])
    } else {
      as.numeric(t[2:length(t)])
    }
  }


  # read in rep and ctl files
  REP <- readLines(here::here(year, model, "retro", "model", folder, paste0(model_name, ".rep")))
  CTL <- readLines(here::here(year, model, "retro", "model", folder,  paste0(dat_name, ".ctl")))
  PSV <- file(here::here(year, model,  "retro", "model", folder, paste0(model_name, ".psv")), "rb")
  STD <- read.delim(here::here(year, model,  "retro", "model", folder, paste0(model_name, ".std")), sep="", header = TRUE)
  mceval <- read.delim(here::here(year, model,  "retro", "model", folder, "evalout.prj"), sep="", header=FALSE)

  # clean rep file
  suppressWarnings(data.frame(year = unlist(strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     dplyr::mutate(year = as.numeric(year)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(year)) -> yrs

  suppressWarnings(data.frame(age = unlist(strsplit(REP[grep("Age", REP)[1]]," "))) %>%
                     dplyr::mutate(age = as.numeric(age)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(age)) -> ages

  styr_rec <- yrs[1] - length(ages) + 2

  suppressWarnings(as.data.frame(cbind(yrs = yrs, ages = ages, styr_rec = styr_rec)) %>%
                     dplyr::mutate(ages = replace(ages, duplicated(ages), NA),
                                   styr_rec = replace(styr_rec, duplicated(styr_rec), NA))) %>%
    write.csv(here::here(year, model, "retro", "model", folder, "processed", "ages_yrs.csv"), row.names = FALSE)

  # MCMC parameters ----

  npar = readBin(PSV, what = integer(), n=1)
  mcmcs = readBin(PSV, what = numeric(), n = (npar * mcmc / mcsave))
  close(PSV)
  mcmc_params = matrix(mcmcs, byrow=TRUE, ncol=npar)
  mcmc_params = mcmc_params[501:nrow(mcmc_params),]
  colnames(mcmc_params) = STD$name[1:ncol(mcmc_params)]
  write.csv(mcmc_params, here::here(year, model,  "retro", "model", folder, "processed", "mcmc.csv"), row.names = FALSE)

  # mceval phase output ----

  #Curry's Change
  mceval = mceval[501:nrow(mceval),]

  #Length colnames = 286
  # columns mcmc_other = 271

  #1-8: Through objective function value

  colnames(mceval) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
                       "ABC", "obj_fun",
                       paste0("tot_biom_", yrs),
                       paste0("log_rec_dev_", seq(styr_rec, yrs[length(yrs)])),
                       paste0("spawn_biom_", yrs),
                       "log_mean_rec",
                       paste0("spawn_biom_proj_", max(yrs) + 1:15),
                       paste0("pred_catch_proj_", max(yrs) + 1:15),
                       paste0("rec_proj_", max(yrs) + 1:10),
                       paste0("tot_biom_proj_", max(yrs)))
  write.csv(mceval, here::here(year, model, "retro", "model", folder, "processed", "mceval.csv"), row.names = FALSE)

  # catch data ----

  pred = strsplit(REP[grep("Pred_Catch", REP)], " ")
  r1 = which(pred[[1]] == "Pred_Catch")
  r2 = which(pred[[1]] == "Pred_catch_later")
  r3 = which(pred[[1]] == "")
  pred = as.numeric(pred[[1]][-c(r1, r2, r3)])

  obs = strsplit(REP[grep("Obs_Catch", REP)], " ")
  r1 = which(obs[[1]] == "Obs_Catch")
  r2 = which(obs[[1]] == "Obs_Catch_Later")
  r3 = which(obs[[1]] == "")
  obs = as.numeric(obs[[1]][-c(r1, r2, r3)])

  data.frame(obs = obs, pred = pred) %>%
    write.csv(here::here(year, model, "retro", "model", folder, "processed", "catch.csv"))



  # survey data ----
  syr = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][2]
  syr = strsplit(syr," ")
  syr = subset(syr[[1]], syr[[1]]!="")
  syr = as.numeric(syr[2:length(syr)])

  obs = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][4]
  obs = strsplit(obs," ")
  obs = subset(obs[[1]], obs[[1]]!="")
  obs = as.numeric(obs[2:length(obs)])

  se = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][6]
  se = strsplit(se," ")
  se = subset(se[[1]], se[[1]]!="")
  se = as.numeric(se[2:length(se)])

  pred = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][3]
  pred = strsplit(pred," ")
  pred = subset(pred[[1]], pred[[1]]!="")
  pred = as.numeric(pred[2:length(pred)])


  data.frame(year = syr, biomass = obs, se = se) %>%
    dplyr::mutate(lci = biomass - 1.96 * se,
                  uci = biomass + 1.96 * se ) %>%
    dplyr::bind_cols(pred = pred) %>%
    write.csv(here::here(year, model, "retro", "model", folder, "processed", "survey.csv"), row.names = FALSE)


  # recruitment ----

  N = REP[grep("Numbers", REP):(grep("Obs_P_fish_age", REP)-2)]
  t = NA
  for(i in 1:length(yrs)){
    ts = as.numeric(strsplit(N[i+1]," ")[[1]][3])
    t = c(t, ts)}
  pred_rec = t[!is.na(t)]

  # biomass & F & recruits ----
  data.frame(year = yrs,
             tot_biom = rep_item("Tot_biom"),
             sp_biom = rep_item("SpBiom"),
             F = rep_item("Fully_selected_F"),
             recruits = pred_rec) %>%
    write.csv(here::here(year, model,  "retro", "model", folder, "processed", "bio_rec_f.csv"), row.names = FALSE)


  # selectivity ----
  data.frame(age = ages,
             fish = rep_item("Fishery_Selectivity"),
             srv1 = rep_item("TWL Survey_Selectivity")) %>%
    write.csv(here::here(year, model,  "retro", "model", folder, "processed", "selex.csv"), row.names = FALSE)

  # yield ratio B40 & B35----

  data.frame(B40 = STD$value[which(STD$name=="B40")],
             B35 = as.numeric(REP[(grep("B_35",REP)+1):(grep("F_40",REP)[1]-1)]),
             yld_rat = as.numeric(unlist(strsplit(CTL[grep("yieldratio", CTL)], "\t"))[1])) %>%
    write.csv(here::here(year, model,  "retro", "model", folder, "processed", "b35_b40_yld.csv"), row.names = FALSE)

}

print("process_retro_dusk")


run_retro_dusk <- function(year, model="m22.3a", tpl_name = 'base', dat_name = "goa_dusk",
                           n_retro = 10, mcmc_on = FALSE, mcmc_num = 100000, mcsave_num = 100){

  if (!dir.exists(here::here(year, model, "retro"))){
    dir.create(here::here(year, model, "retro", "model"), recursive=TRUE)
    dir.create(here::here(year, model, "retro", "results"), recursive=TRUE)
  }

  file.copy(here::here(year, model, paste0(tpl_name, ".tpl")),
            here::here(year, model, "retro", "model"),
            overwrite = TRUE)

  file.copy(here::here(year, model, "MAT.dat"),
            here::here(year, model, "retro", "model"))

  file.copy(here::here(year, model, paste0(tpl_name, ".exe")),
            here::here(year, model, "retro", "model"),
            overwrite = TRUE)

  if(file.exists(here::here(year, model, paste0(tpl_name, ".pin"))) ) {
    file.copy(here::here(year, model, paste0(tpl_name, ".pin")),
              here::here(year, model, "retro", "model", paste0(tpl_name, "_og.pin")),
              overwrite = TRUE)
    PIN_og = readLines(here::here(year, model, "retro", "model", paste0(tpl_name, "_og.pin")), warn = FALSE)
  }


  # model .dat, ctl, pin files
  CTL = read.delim(here::here(year, model, paste0(dat_name, ".ctl") ), header = FALSE )
  DAT = readLines(here::here(year, model, paste0(dat_name, "_", year, ".dat")), warn = FALSE)

  # define .dat file breaks
  st_end = data.frame(st = grep("#-", DAT),
                      end = grep("#!", DAT))

  # model dims
  styr = as.numeric(DAT[st_end$st[2] - 3]) # start of model (example 1961 for POP)
  nages = as.numeric(DAT[st_end$st[2] + 3]) # number of age bins
  nlens = as.numeric(DAT[st_end$st[2] + 7]) # number of length bins
  yrs_retro = seq(year - n_retro + 1, year) # end years for retro run

  # Set up some results files
  res_sb <- matrix(nrow = length(seq(styr, year)), ncol = n_retro)
  rownames(res_sb) <- seq(styr, year)
  colnames(res_sb) <- seq( year - n_retro + 1, year)
  res_rec <- matrix(nrow = length(seq(styr, year)), ncol = n_retro)
  rownames(res_rec) <- seq(styr, year)
  colnames(res_rec) <- seq(year - n_retro + 1, year)
  res_tb <-  matrix(nrow = length(seq(styr, year)), ncol = n_retro)
  rownames(res_tb) <- seq(styr, year)
  colnames(res_tb) <- seq( year - n_retro + 1, year)

  # retro data loop
  T_start <- Sys.time() #Timer start

  for(y in 1:n_retro){
    # Set endyr
    endyr = yrs_retro[y]
    nyrs = endyr - styr + 1

    dat_retro = c(DAT[st_end$st[1]:st_end$end[1]],
                  # end year of model
                  as.character(endyr), # end year
                  DAT[st_end$st[2]:st_end$end[2]],
                  # fishery catch
                  paste(scan(text = DAT[st_end$st[3] - 1])[1:nyrs],collapse=" "),
                  DAT[st_end$st[3]:st_end$end[3]],
                  # trawl survey biomass
                  as.character(length(which(scan(text = DAT[st_end$st[5]-1]) <= endyr))),
                  DAT[st_end$st[4]:st_end$end[4]],
                  paste(scan(text = DAT[st_end$st[5] - 1])[1:length(which(scan(text = DAT[st_end$st[5]-1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[5]:st_end$end[5]],
                  paste(scan(text = DAT[st_end$st[6] - 1])[1:length(which(scan(text = DAT[st_end$st[5]-1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[6]:st_end$end[6]],
                  paste(scan(text = DAT[st_end$st[7] - 1])[1:length(which(scan(text = DAT[st_end$st[5]-1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[7]:st_end$end[7]],
                  paste(scan(text = DAT[st_end$st[8] - 1])[1:length(which(scan(text = DAT[st_end$st[5]-1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[8]:st_end$end[8]],
                  paste(scan(text = DAT[st_end$st[9] - 1])[1:length(which(scan(text = DAT[st_end$st[5]-1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[9]:st_end$end[9]],
                  # fishery age comp (note: will never have ages in the final year of model)
                  as.character(length(which(scan(text = DAT[st_end$st[11]-1]) < (endyr - 1)))),
                  DAT[st_end$st[10]:st_end$end[10]],
                  paste(scan(text = DAT[st_end$st[11] - 1])[1:length(which(scan(text = DAT[st_end$st[11]-1]) < (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[11]:st_end$end[11]],
                  paste(scan(text = DAT[st_end$st[12] - 1])[1:length(which(scan(text = DAT[st_end$st[11]-1]) < (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[12]:st_end$end[12]],
                  paste(scan(text = DAT[st_end$st[13] - 1])[1:length(which(scan(text = DAT[st_end$st[11]-1]) < (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[13]:st_end$end[13]],
                  paste(scan(text = DAT[st_end$st[14] - 1])[1:length(which(scan(text = DAT[st_end$st[11]-1]) < (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[14]:st_end$end[14]],
                  DAT[(st_end$end[14] + 1):(st_end$end[14] + length(which(scan(text = DAT[st_end$st[11] - 1]) < (endyr - 1))))],
                  DAT[st_end$st[15]:st_end$end[15]],
                  # trawl survey age comp (note: will never have ages in the final year of model)
                  as.character(length(which(scan(text = DAT[st_end$st[17] - 1]) <= (endyr - 1)))),
                  DAT[st_end$st[16]:st_end$end[16]],
                  paste(scan(text = DAT[st_end$st[17] - 1])[1:length(which(scan(text = DAT[st_end$st[17] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[17]:st_end$end[17]],
                  paste(scan(text = DAT[st_end$st[18] - 1])[1:length(which(scan(text = DAT[st_end$st[17] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[18]:st_end$end[18]],
                  paste(scan(text = DAT[st_end$st[19] - 1])[1:length(which(scan(text = DAT[st_end$st[17] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[19]:st_end$end[19]],
                  paste(scan(text = DAT[st_end$st[20] - 1])[1:length(which(scan(text = DAT[st_end$st[17] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[20]:st_end$end[20]],
                  DAT[(st_end$end[20] + 1):(st_end$end[20] + length(which(scan(text = DAT[st_end$st[17] - 1]) <= (endyr - 1))))],
                  DAT[st_end$st[21]:st_end$end[21]],
                  # fishery size comp
                  as.character(length(which(scan(text = DAT[st_end$st[23] - 1]) <= (endyr - 1)))),
                  DAT[st_end$st[22]:st_end$end[22]],
                  paste(scan(text = DAT[st_end$st[23] - 1])[1:length(which(scan(text = DAT[st_end$st[23] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[23]:st_end$end[23]],
                  paste(scan(text = DAT[st_end$st[24] - 1])[1:length(which(scan(text = DAT[st_end$st[23] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[24]:st_end$end[24]],
                  paste(scan(text = DAT[st_end$st[25] - 1])[1:length(which(scan(text = DAT[st_end$st[23] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[25]:st_end$end[25]],
                  paste(scan(text = DAT[st_end$st[26] - 1])[1:length(which(scan(text = DAT[st_end$st[23] - 1]) <= (endyr - 1)))],collapse=" "),
                  DAT[st_end$st[26]:st_end$end[26]],
                  DAT[(st_end$end[26] + 1):(st_end$end[26] + length(which(scan(text = DAT[st_end$st[23] - 1]) <= (endyr - 1))))],
                  DAT[st_end$st[27]:st_end$end[27]],
                  # survey size comp (will have size comps in final year of model)
                  as.character(length(which(scan(text = DAT[st_end$st[29] - 1]) <= endyr))),
                  DAT[st_end$st[28]:st_end$end[28]],
                  paste(scan(text = DAT[st_end$st[29] - 1])[1:length(which(scan(text = DAT[st_end$st[29] - 1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[29]:st_end$end[29]],
                  paste(scan(text = DAT[st_end$st[30] - 1])[1:length(which(scan(text = DAT[st_end$st[29] - 1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[30]:st_end$end[30]],
                  paste(scan(text = DAT[st_end$st[31] - 1])[1:length(which(scan(text = DAT[st_end$st[29] - 1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[31]:st_end$end[31]],
                  paste(scan(text = DAT[st_end$st[32] - 1])[1:length(which(scan(text = DAT[st_end$st[29] - 1]) <= endyr))],collapse=" "),
                  DAT[st_end$st[32]:st_end$end[32]],
                  DAT[(st_end$end[32] + 1):(st_end$end[32] + length(which(scan(text = DAT[st_end$st[29] - 1]) <= endyr)))],
                  DAT[st_end$st[33]:st_end$end[33]]
                  )


    # Write data and control file
    write.table(dat_retro,
                file = here::here(year, model, "retro", "model", paste0(dat_name, "_", endyr, ".dat")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    CTL_retro = as.matrix(CTL)
    CTL_retro[2,1] = paste0(dat_name, "_", endyr, ".dat")
    CTL_retro[5,1] = as.character(endyr)

    #rewrites CTL file with retro year data name and end year
    write.table(CTL_retro,
                file = here::here(year, model, "retro", "model", paste0(dat_name, ".ctl")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    # update pin file and save if exists
    if(file.exists(here::here(year, model, "retro", "model", paste0(tpl_name, "_og.pin") ) ) ) {

      PIN_retro = as.matrix(PIN_og)
      PIN_retro[35] = paste("",scan(text = PIN_retro[35])[1:nyrs],collapse=" ")
      pin_log_rec_dev = scan(text = PIN_retro[37])
      PIN_retro[37] = paste("",scan(text = PIN_retro[37])[1:(length(pin_log_rec_dev)-(year-endyr))],collapse=" ")

      write.table(PIN_retro,
                  file= here::here(year, model, "retro", "model", paste0(tpl_name, ".pin")),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }

    # run retro models

    ## set folder directory
    setwd(here::here(year, model, "retro", "model"))

    ## to run with mcmc ON
    if(isTRUE(mcmc_on)) {
      system(paste0(tpl_name,'.exe', ' -mcmc ', mcmc_num, ' -mcsave ', mcsave_num))
      system(paste0(tpl_name,'.exe',' -mceval'))

      file.copy("evalout.prj",
                here::here(year, model, "retro", "results", paste0("mcmc_", endyr, ".std")), overwrite = TRUE)
    }

    ## to run with mcmc OFF
    if(isFALSE(mcmc_on)) {
      system(paste0(tpl_name,'.exe'))
    }

    file.copy(paste0(tpl_name, ".std"),
              here::here(year, model, "retro", "results", paste0("std_", endyr, ".std")), overwrite = TRUE)

    file.copy(here::here(year, model, "retro", "model", paste0(tpl_name, ".rep")),
              here::here(year, model, "retro", "results", paste0("rep_", endyr, ".rep")), overwrite = TRUE)

    file.copy(here::here(year, model, "retro", "model", paste0(tpl_name, ".par")),
              here::here(year, model, "retro", "results", paste0("rep_", endyr, ".par")), overwrite = TRUE)

    # # sometimes I get weird output from admb it drops the model name and produces these files

    file.copy(here::here(year, model, "retro", "model", "update~1.rep"),
              here::here(year, model, "retro", "results", paste0("rep_", endyr, ".rep")), overwrite = TRUE)

    file.copy("update~1.std",
              here::here(year, model, "retro", "results", paste0("std_", endyr, ".std")), overwrite = TRUE)

    # compile ssb/recruitment results
    std = read.delim(here::here(year, model, "retro", "model", paste0(tpl_name, ".std")),
                     header = T, sep = "")
    res_sb[1:nyrs, y] = std$value[which(std$name=="spawn_biom")]
    res_rec[1:nyrs, y] = std$value[which(std$name=="pred_rec")]
    res_tb[1:nyrs, y] = std$value[which(std$name=="tot_biom")]

  } # end retro for loop

  write.csv(res_sb, here::here(year, model, 'processed', 'retro_ssb.csv'))
  write.csv(res_rec, here::here(year, model, 'processed', 'retro_rec.csv'))
  write.csv(res_tb, here::here(year, model, 'processed', 'retro_totbiom.csv'))

  # print total time
  T_end<-Sys.time()
  T_end-T_start

} # end run_retro_dusk

print("run_retro_dusk")

report_admb_dusk <- function(report_name, std_name , model_year = year_mod, year = year, model) {

  # create directory
  if (!dir.exists(here::here(year, model, "processed"))){
    dir.create(here::here(year, model, "processed"), recursive=TRUE)
  }

  temp.rep <- readLines(here::here(year, model, 'base.rep'))
  temp.std <- readLines(here::here(year, model, 'base.std'))[-1]

  # single values
  num_par <- as.numeric(unlist(strsplit(temp.rep[grep("Number parameters estimated",temp.rep)+1],split=" ")))
  B100 <- round(as.numeric(unlist(strsplit(temp.rep[grep("B_zero",temp.rep)+1],split=" "))),digits=0) # B0 = B100
  B40 <- round(as.numeric(unlist(strsplit(temp.rep[grep("B_40",temp.rep)+1],split=" "))),digits=2) #Get B40 from report file
  F40 <- round(as.numeric(unlist(strsplit(temp.rep[grep("F_40",temp.rep)+1],split=" "))), digits=3) [1]
  F.ABC <- round(as.numeric(unlist(strsplit(temp.rep[grep("F_ABC",temp.rep)+1],split=" "))), digits=3)[1]
  ABC <- round(as.numeric(unlist(strsplit(temp.rep[grep("ABC for",temp.rep)+1],split=" "))), digits=0)[3]   # first 2 are the F_ABC values, this grabs ABC for first projection year
  q <- as.numeric(unlist(strsplit(temp.rep[grep("q_trawl",temp.rep)+1],split=" ")))[1]
  avg_log_rec <- round(as.numeric(unlist(strsplit(temp.rep[grep("log_mean_rec",temp.rep)+1],split=" "))), digits=3)[1]
  total_B <- round(as.numeric(unlist(strsplit(temp.rep[grep("tot_biom",temp.rep)+1],split=" "))), digits=0)[1]
  spawn_B <- round(as.numeric(unlist(strsplit(temp.rep[grep("spawn_biom",temp.rep)+1],split=" "))), digits=0)[1] # SSB for model_year+1
  sigmar <- as.numeric(unlist(strsplit(temp.rep[grep("SigmaR: ",temp.rep)],split=" ")))[2]

  # vectors
  yrs_vec <- 1977:model_year

  SSB_vec <- na.omit(round(as.numeric(unlist(strsplit(temp.rep[grep("SpBiom",temp.rep)],split=" "))), 0 ) )
  SSB_df <- data.frame(year= yrs_vec, SSB = SSB_vec)

  TotB_vec <- na.omit(round(as.numeric(unlist(strsplit(temp.rep[grep("Tot_biom",temp.rep)],split=" "))), 0) )
  TotB_df <- data.frame(year= yrs_vec, Tot_B = TotB_vec)

  #TotB_vec <-

  surv_yrs <- na.omit(as.numeric(unlist(strsplit(temp.rep[grep("Survey Biomass",temp.rep)+1],split=" "))))
  surv_yrs <- surv_yrs[1:length(surv_yrs)-1] # weirdly repeating the 1999 (need to go back and look at ADMB report)
  surveyB_pred_vec <- na.omit(as.numeric(unlist(strsplit(temp.rep[grep("Survey Biomass",temp.rep)+2],split=" "))))
  surveyB_pred_vec <- surveyB_pred_vec[1:length(surveyB_pred_vec)-1] # weirdly repeating the 1999 (need to go back and look at ADMB report)
  surveyB_pred_df <- data.frame(year = surv_yrs, survB_pred = surveyB_pred_vec)

  surveyB_obs_vec <- na.omit(as.numeric(unlist(strsplit(temp.rep[grep("Survey Biomass",temp.rep)+3],split=" "))))
  surveyB_obs_vec <- surveyB_obs_vec[1:length(surveyB_obs_vec)-1] # weirdly repeating the 1999 (need to go back and look at ADMB report)
  surveyB_obs_df <- data.frame(year = surv_yrs, survB_obs = surveyB_obs_vec)

  surveyB_df <- full_join(surveyB_pred_df, surveyB_obs_df)

  # recruitment
  # Need sd?? neet to get them from teh .std file
  N_tmp <- temp.rep[grep("Numbers", temp.rep):(grep("Obs_P_fish_age", temp.rep)-2)]
  t_tmp = NA
  for(i in 1:length(yrs_vec)){
    ts = as.numeric(base::strsplit(N_tmp[i+1]," ")[[1]][3])
    t_tmp = c(t_tmp, ts)
  }
  pred_rec <- data.frame(year = yrs_vec, recruit= t_tmp[!is.na(t_tmp)])

  # fully selected f
  rep_item <- function(name){
    t <- strsplit(temp.rep[grep(name, temp.rep)]," ")
    t <- subset(t[[1]], t[[1]]!="")
    if(t[[1]][1] == "TWL"){
      as.numeric(t[3:length(t)])
    } else {
      as.numeric(t[2:length(t)])
    }
  }

  F_df <- data.frame( year = yrs_vec,  F_vec =rep_item("Fully_selected_F"))

  #data.frame(year = yrs,
  #           tot_biom = rep_item("Tot_biom"),
  #           sp_biom = rep_item("SpBiom"),
  #           F = rep_item("Fully_selected_F"),
  #           recruits = pred_rec) %>%
  #  write.csv(here::here(year, model, "processed", "bio_rec_f.csv"), row.names = FALSE)

  # Likelihoods
  tot_ll <-  data.frame(name = "total_LL",
                        LL=as.numeric(unlist(strsplit(temp.rep[grep("Total likelihood",temp.rep)+1],split=" "))))
  data_ll <-  data.frame(name = "data_LL" ,
                         LL=as.numeric(unlist(strsplit(temp.rep[grep("Data likelihood",temp.rep)+1],split=" "))))

  wts_ll_tmp <- data.frame(tmp= temp.rep[(grep("Wts_and_Likelihoods", temp.rep)+1):(grep("SigmaR:", temp.rep)-1)] )
  ll_df_tmp <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "_Likelihood")) %>%
    tidyr::separate(tmp, into= c("wt", "LL", "name1", "name2"), sep = "\\s") %>%
    tidyr::unite("name", name1, name2, na.rm=T) %>%
    dplyr::mutate(LL = as.numeric(LL)) %>%
    dplyr::filter(!LL==0)
  #ll_names <- c("Catch", "Trawl survey", "Fishery ages", "Survey ages", "Fishery lengths", "Recruitment devs")
  #ll_df <- bind_cols(ll_df, ll_names = ll_names)

  sigmr_prior <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "priors sigr")) %>%
    tidyr::separate(tmp, into= c( "wt" ,"LL", "name1", "name2"), sep = "\\s") %>%
    tidyr::unite("name", name1, name2, sep = " ") %>%
    dplyr::mutate(LL = as.numeric(LL))
  q_prior <- wts_ll_tmp %>%
    dplyr::filter(str_detect(tmp, "priors q TWL")) %>%
    tidyr::separate(tmp, into= c( "wt" ,"LL", "name1", "name2"), sep = "\\s") %>%
    tidyr::unite("name", name1, name2, sep = " ") %>%
    dplyr::mutate(LL = as.numeric(LL))

  # Create likelihood large table with names and values from model name

  ll_names <- c("Catch", "Trawl survey", "Fishery ages", "Survey ages", "Fishery lengths",
                "Recruitment devs", "sigmaR", "q prior", "Data LL", "Total LL")
  ll_df <- bind_rows(ll_df_tmp, sigmr_prior, q_prior, data_ll, tot_ll) %>% select(-wt)
  ll_df <- bind_cols(ll_df, Likelihood = ll_names) %>%
    dplyr::select(Likelihood, LL) %>%
    dplyr::mutate(model_name = paste0(model))

  write.csv(ll_df, here::here(year, model, "processed", "Table_likelihood_values.csv"), row.names = FALSE)

  # Create parameter estimate table

  param_names <- c("sigmaR" , "q", "avg rec", "F40", "Total Biomass", "SSB", "B100", "B40", "ABC")
  par_est_values <- c(sigmar, q, exp(avg_log_rec), F40, total_B, spawn_B, B100, B40, ABC)
  par_est_df <- bind_cols(Parameter = param_names, Estimates = par_est_values) %>%
    dplyr::mutate(model_name = paste0(model))

  write.csv(par_est_df, here::here(year, model, "processed", "Table_parameter_est.csv"), row.names = FALSE)

  # format std to readable form

  std_base <- data.frame(text = temp.std) %>%
    dplyr::mutate(text = stringr::str_squish(text)) %>%
    tidyr::separate(text, c("index", "name", "value", "std_dev"), sep = " ")

  write.csv(std_base, here::here(year, model, "processed", "std_output.csv"), row.names = FALSE )

  # Vectors to save
  write.csv(SSB_df, here::here(year, model, "processed", "df_SSB.csv"), row.names = FALSE)
  write.csv(TotB_df, here::here(year, model, "processed", "df_totalB.csv"), row.names = FALSE)
  write.csv(surveyB_df, here::here(year, model, "processed", "df_surveyB.csv"), row.names = FALSE)
  write.csv(pred_rec, here::here(year, model, "processed", "df_pred_recruit.csv"), row.names = FALSE)
  write.csv(F_df, here::here(year, model, "processed", "df_fully_selected_F.csv"), row.names = FALSE)

}

print("report_admb_dusk")

# apportionment split

### still need to redo for dusky

run_apport_dusk <- function(year, model, biom_dat, biom_name = "vast", frac_dat){

  if (!dir.exists(here::here(year, model, "apport"))){
    dir.create(here::here(year, model, "apport"), recursive=TRUE)
  }

  # run rema model for e/c/wgoa, save output and plots into 'apport' folder
  biomass_dat <- biom_dat %>% # use area-specific dataframe
    dplyr::mutate(cv = se / biomass,
                  lci = (biomass- (se*1.96)),
                  uci =  (biomass+ (se*1.96))) %>%
    dplyr::filter(year >= 1990)

  if(biom_name %in% c("GAP", "gap")) {
    # db_gap_index
    db_index <- biomass_dat %>%
      dplyr::select(year, strata= area, biomass, cv, lci, uci) %>%
      dplyr::mutate(strata = factor(strata, levels = c("western", "central", "eastern"))) %>%
      dplyr::ungroup()

    input_db_rema <- prepare_rema_input( model_name = 'db_rema',
                                         biomass_dat = db_index,
                                         end_year = max(db_index$year, na.rm=T),
                                         zeros = list(assumption = 'NA' ) )
    mod_db_rema <- fit_rema(input_db_rema)
    out_db_rema <- tidy_rema(mod_db_rema)

    prop_out <- out_db_rema$proportion_biomass_by_strata
  } # end gap proportions

  # for vast apportionment
  if(sum(stringr::str_detect(biom_name, pattern = c("VAST", "vast")))>0 ) {

    vast_prop <- biomass_dat %>%
      dplyr::select(year, area, biomass) %>%
      dplyr::group_by(year) %>%
      dplyr::mutate(proporions = biomass/ sum(biomass))

    prop_out <- vast_prop %>%
      tidyr::pivot_wider(id_cols = -biomass, names_from = area, values_from = proporions) %>%
      dplyr::mutate(model_name = biom_name)

  } # end vast proportions

  # compute wyak/eyak-se split

  # function to calculate wyak fractions each survey year
  split_function <- function(data_input, year_input) {

    data_input_sub <- data_input %>%
      dplyr::filter(year %in% c( (max(year_input) - 5):max(year_input)) )

    wyak_split_tmp <- data_input_sub %>%
      dplyr::mutate(wt = c(4,6,9),
                    wtd_var = (wt / sum(wt))^2 * var(wyak_fract) ) %>%
      dplyr::summarise(sum_wtd_var = sum(wtd_var),
                       wtd_avg = sum(wyak_fract * wt) / sum(wt)) %>%
      dplyr::mutate(wyak = round(wtd_avg + 2 * sqrt(sum_wtd_var), digits = 2),
                    year = max(data_input_sub$year) )
  }

  # running function to compute wyak fractions and storing all results
  wyak_split_store <- data.frame()
  yr_vec <- seq(2008, year, by= 2)

  for(y in 1:length(yr_vec)) {
    year_input_num <- yr_vec[y]
    wyak_split_tmp_store <- split_function(data_input = frac_dat, year_input = year_input_num)
    wyak_split_store <- dplyr::bind_rows(wyak_split_store, wyak_split_tmp_store)
  }

  #wyak_split_store

  prop_egoa_out <- prop_out %>%
    dplyr::full_join(wyak_split_store %>%
                       dplyr::select(year, wyak)) %>%
    dplyr::mutate(seo = 1-wyak)

  write.csv(prop_egoa_out , here::here(year, model, "apport", paste0("apport_prop_", biom_name, ".csv")), row.names = FALSE)
}

print("run_apport_dusk")
