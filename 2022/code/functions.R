process_results <- function(year, model, model_name = alt, dat_name,
                            rec_age, plus_age, mcmc, mcsave, len_bins){

  # setup
  if (!dir.exists(here::here(year, model, "processed"))){
    dir.create(here::here(year, model, "processed"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, model, "figs"))){
    dir.create(here::here(year, model, "figs"), recursive=TRUE)
  }

  if (!dir.exists(here::here(year, model, "tables"))){
    dir.create(here::here(year, model, "tables"), recursive=TRUE)
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
  REP <- readLines(here::here(year, model, paste0(model_name, ".rep")))
  CTL <- readLines(here::here(year, model, paste0(dat_name, ".ctl")))
  PSV <- file(here::here(year, model, paste0(model_name, ".psv")), "rb")
  STD <- read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)
  mceval <- read.delim(here::here(year, model, "evalout.prj"), sep="", header=FALSE)

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
    write.csv(here::here(year, model, "processed", "ages_yrs.csv"), row.names = FALSE)

  # MCMC parameters ----

  npar = readBin(PSV, what = integer(), n=1)
  mcmcs = readBin(PSV, what = numeric(), n = (npar * mcmc / mcsave))
  close(PSV)
  mcmc_params = matrix(mcmcs, byrow=TRUE, ncol=npar)
  mcmc_params = mcmc_params[501:nrow(mcmc_params),]
  colnames(mcmc_params) = STD$name[1:ncol(mcmc_params)]
  write.csv(mcmc_params, here::here(year, model, "processed", "mcmc.csv"), row.names = FALSE)

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
  write.csv(mceval, here::here(year, model, "processed", "mceval.csv"), row.names = FALSE)

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
    write.csv(here::here(year, model, "processed", "catch.csv"))



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
    write.csv(here::here(year, model, "processed", "survey.csv"), row.names = FALSE)


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
    write.csv(here::here(year, model, "processed", "bio_rec_f.csv"), row.names = FALSE)


  # selectivity ----
  data.frame(age = ages,
             fish = rep_item("Fishery_Selectivity"),
             srv1 = rep_item("TWL Survey_Selectivity"),
             maturity = rep_item("Maturity")) %>%
    write.csv(here::here(year, model, "processed", "selex.csv"), row.names = FALSE)

  # yield ratio B40 & B35----

  data.frame(B40 = STD$value[which(STD$name=="B40")],
             B35 = as.numeric(REP[(grep("B_35",REP)+1):(grep("F_40",REP)[1]-1)]),
             yld_rat = as.numeric(unlist(strsplit(CTL[grep("yieldratio", CTL)], "\t"))[1])) %>%
    write.csv(here::here(year, model, "processed", "b35_b40_yld.csv"), row.names = FALSE)

  # size comps ----

  #! this will need a switch for multiple surveys

  obs = REP[grep("Obs_P_fish_age",REP):(grep("Pred_P_fish_age",REP)-2)]
  pred = REP[grep("Pred_P_fish_age",REP):(grep("Obs_P_fish_size",REP)-2)]

  obs_l = REP[grep("Obs_P_fish_size",REP):(grep("Pred_P_fish_size",REP)-2)]
  pred_l = REP[grep("Pred_P_fish_size",REP):(grep("Obs_P_srv1_age",REP)-2)]

  s_obs = REP[grep("Obs_P_srv1_age",REP):(grep("Pred_P_srv1_age",REP)-2)]
  s_pred = REP[grep("Pred_P_srv1_age",REP):(grep("Obs_P_srv1_size",REP)-2)]

  s_obs_l = REP[grep("Obs_P_srv1_size",REP):(grep("Pred_P_srv1_size",REP)-2)]

  rockfishr::purrit(obs, pred, rec_age, plus_age, comp = "age", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "fac.csv"))

  rockfishr::purrit(obs_l, pred_l, rec_age, plus_age, comp = "length", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "fsc.csv"))

  rockfishr::purrit(s_obs, s_pred, rec_age, plus_age, comp = "age", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "sac.csv"))

  rockfishr::purrit(s_obs_l, pred = NULL, rec_age, plus_age, comp = "length", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "ssc.csv"))


}
process_retro <- function(year, model, model_name = alt, dat_name,
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


plot_compare_survey <- function(year, models = c('2022, default', '2022, lognormal')) {

  if (!dir.exists(here::here(year, "compare_models"))){
    dir.create(here::here(year, "compare_models"))
  }

  dat = data.frame()
  m = list(rep(NA, length(models)))
  for(i in 1:length(models)) {
    m[[i]] = scan(text = models[i], sep = ",", what = "")
    m[[i]][2] = gsub(" ", "", m[[i]][2])
    # m[[i]][3] = gsub(" ", "", m[[i]][3])

    year = m[[i]][1]
    # folder = m[[i]][2]
    model = m[[i]][2]

    id = model
    # if(model=="db"){
    #   id = "design-based"
    # } else if(model=="m15.5a"){
    #   id = "A"
    # } else if(model=="pois_gamma_750"){
    #   id = "B"
    # } else if(model=="log_1000"){
    #   id = "C"
    # } else if(model=="log_750"){
    #   id = "D"
    # } else if(model=="pois_log_500"){
    #   id = "E"
    # } else {
    #   id = "F"
    # }

    dat %>%
      dplyr::bind_rows(
        read.csv(here::here(year,  model, "processed", "survey.csv")) %>%
          dplyr::rename_all(tolower) %>%
          dplyr::select(year = starts_with("y"),
                        Observed = starts_with("bio"),
                        Predicted = pred,
                        se, lci, uci) %>%
          tidyr::pivot_longer(-c(year, se, uci, lci)) %>%
          dplyr::mutate(value = value / 1000,
                        uci = uci / 1000,
                        lci = lci / 1000,
                        model = id)) -> dat
  }

  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = model)) +
    ggplot2::geom_point(data = dplyr::filter(dat, name == "Observed"), position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_errorbar(data = dplyr::filter(dat, name == "Observed"),
                           ggplot2::aes(ymin = lci, ymax = uci, color = model), width = 0.4, position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_line(data = dplyr::filter(dat, name == "Predicted"),
                       ggplot2::aes(color = model)) +
    scico::scale_color_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year)$breaks,
                                labels = funcr::tickr(dat, year)$labels) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    funcr::theme_report() +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.8,0.70))

  # ggplot2::ggsave(here::here(year, "compare_models", "figs", "compare_est_survey.png"),
  # width = 6.5, height = 6.5, units = "in", dpi = 200)
}
plot_compare_biomass <- function(year, models = c('2022, lognormal')) {

  if (!dir.exists(here::here(year, "compare_models"))){
    dir.create(here::here(year, "compare_models"))
  }

  dat = data.frame()
  m = list(rep(NA, length(models)))
  for(i in 1:length(models)) {
    m[[i]] = scan(text = models[i], sep = ",", what = "")
    m[[i]][2] = gsub(" ", "", m[[i]][2])
    # m[[i]][3] = gsub(" ", "", m[[i]][3])

    year = m[[i]][1]
    # folder = m[[i]][2]
    model = m[[i]][2]
    id = model
    # if(model=="db"){
    #   id = "design-based"
    # } else if(model=="m15.5a"){
    #   id = "A"
    # } else if(model=="pois_gamma_750"){
    #   id = "B"
    # } else if(model=="log_1000"){
    #   id = "C"
    # } else if(model=="log_750"){
    #   id = "D"
    # } else if(model=="pois_log_500"){
    #   id = "E"
    # } else {
    #   id = "F"
    # }

    yrs = read.csv(here::here(year, model, "processed", "ages_yrs.csv"))$yrs
    bio = read.csv(here::here(year, model, "processed", "bio_rec_f.csv"))

    dat %>%
      dplyr::bind_rows(
        read.csv(here::here(year, model, "processed", "mceval.csv"))  %>%
          dplyr::select(paste0("tot_biom_", yrs)) %>%
          dplyr::mutate(group = 1:dplyr::n()) %>%
          tidyr::pivot_longer(-group) %>%
          dplyr::mutate(year = as.numeric(gsub("tot_biom_", "", name)),
                        name = "Total biomass") %>%
          dplyr::bind_rows( read.csv(here::here(year,model, "processed", "mceval.csv")) %>%
                              dplyr::select(paste0("spawn_biom_", yrs)) %>%
                              dplyr::mutate(group = 1) %>%
                              tidyr::pivot_longer(-group) %>%
                              dplyr::mutate(year = as.numeric(gsub("spawn_biom_", "", name)),
                                            name = "Spawning biomass")) %>%
          dplyr::mutate(name = factor(name, levels = c("Total biomass", "Spawning biomass"))) %>%
          dplyr::group_by(year, name) %>%
          dplyr::summarise(median = median(value) / 1000,
                           lci = quantile(value, 0.025) / 1000,
                           uci = quantile(value, 0.975) / 1000) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(data.frame(year = yrs,
                                      tot = bio$tot_biom / 1000,
                                      bio = bio$sp_biom / 1000)) %>%
          dplyr::mutate(biomass = ifelse(name == "Total biomass", tot, bio),
                        model = id) %>%
          dplyr::select(-tot, -bio)) -> dat
  }

  dummy = data.frame(year = rep(unique(dat$year),4),
                     name = rep(c("Total biomass", "Spawning biomass"), each = 2 * length(unique(dat$year))),
                     biomass = c(rep(0, length(unique(dat$year))), rep(160, length(unique(dat$year))),
                                 rep(0, length(unique(dat$year))), rep(60, length(unique(dat$year)))),
                     model = NA)

  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, biomass, color = model, fill = model)) +
    ggplot2::geom_blank(data = dummy) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
    ggplot2::facet_wrap(~name, dir = "v", scales = "free_y") +
    ggplot2::scale_y_continuous(name = "Biomass (kt)", labels = scales::comma) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_x_continuous(name = "Year",
                                breaks = funcr::tickr(dat, year, 10, start = 1960)$breaks,
                                labels = funcr::tickr(dat, year, 10, start = 1960)$labels) +
    scico::scale_color_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    scico::scale_fill_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    funcr::theme_report() +
    ggplot2::theme(legend.position = c(0.2, .8))

  # ggplot2::ggsave(here::here(year, "compare_models", "figs", "compare_est_biomass.png"),
  # width = 6.5, height = 6.5, units = "in", dpi = 200)
}

run_retro <- function(year, model="m22.3a", tpl_name = 'base', n_retro = 10, admb_home = NULL, mcon = NULL, mcmc = 100000, mcsave = 100){

  if (!dir.exists(here::here(year, model, "retro"))){
    dir.create(here::here(year, model, "retro", "model"), recursive=TRUE)
    dir.create(here::here(year, model, "retro", "results"), recursive=TRUE)
  }

  file.copy(here::here(year, model, paste0(tpl_name, ".tpl")),
            here::here(year, model, "retro", "model"),
            overwrite = TRUE)

  file.copy(here::here(year, model, "mat.dat"),
            here::here(year, model, "retro", "model"))

  if(is.null(admb_home)){
    R2admb::setup_admb()
  } else {
    R2admb::setup_admb(admb_home)
  }

  setwd(here::here(year, model, "retro", "model"))

  #Compile the Model
  R2admb::compile_admb(tpl_name)

  setwd(here::here())

  # model .dat and ctl files
  CTL = read.delim(here::here(year, model, paste0("goa_dr_", year, ".ctl")), header = FALSE)
  DAT = readLines(here::here(year, model, paste0("goa_dr_", year, ".dat")), warn = FALSE)


  # define .dat file breaks
  st_end = data.frame(sec_st = grep("#!", DAT),
                      sec_end = grep("#-", DAT))


  # model dims
  styr = as.numeric(DAT[st_end$sec_st[1] - 1]) # start of model (example 1961 for POP)
  nages = as.numeric(DAT[st_end$sec_st[1] + 5]) # number of age bins
  nlens = as.numeric(DAT[st_end$sec_st[1] + 9]) # number of length bins

  # retro data loop

  for(y in 1:n_retro){
    # Set endyr
    yrs_retro = seq(year - n_retro + 1, year)
    endyr = yrs_retro[y]
    nyrs = endyr - styr + 1

    DAT_retro = c(DAT[st_end[1,1]:st_end[1,2]], as.character(endyr), DAT[st_end[2,1]:st_end[2,2]])

    # Fishery catch
    DAT_retro = c(DAT_retro,
                  paste(scan(text=DAT[st_end$sec_st[2] + 1])[1:nyrs], collapse=" "),
                  DAT[48])


    # Trawl survey biomass
    BTSb_yrs = length(which(scan(text=DAT[st_end$sec_st[5] - 1]) <= endyr))
    DAT_retro<-c(
      DAT_retro,
      as.character(BTSb_yrs),
      DAT[st_end[4,1]:st_end[4,2]],
      paste(scan(text=DAT[st_end$sec_st[5]-1])[1:BTSb_yrs],collapse=" "),
      DAT[st_end[5,1]:st_end[5,2]],
      paste(scan(text=DAT[st_end$sec_st[6]-1])[1:BTSb_yrs],collapse=" "),
      DAT[st_end[6,1]:st_end[6,2]],
      paste(scan(text=DAT[st_end$sec_st[7]-1])[1:BTSb_yrs],collapse=" "),
      DAT[st_end[7,1]:st_end[7,2]],
      paste(scan(text=DAT[st_end$sec_st[8]-1])[1:BTSb_yrs],collapse=" "),
      DAT[st_end[8,1]:st_end[8,2]],
      paste(scan(text=DAT[st_end$sec_st[9]-1])[1:BTSb_yrs],collapse=" "),
      DAT[st_end[9,1]:st_end[9,2]])

    # Fish age comp
    FAC_yrs<-length(which(scan(text=DAT[st_end$sec_st[11]-1])<(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(FAC_yrs),
                 DAT[st_end[10,1]:st_end[10,2]],
                 paste(scan(text=DAT[st_end$sec_st[11]-1])[1:FAC_yrs],collapse=" "),
                 DAT[st_end[11,1]:st_end[11,2]],
                 paste(scan(text=DAT[st_end$sec_st[12]-1])[1:FAC_yrs],collapse=" "),
                 DAT[st_end[12,1]:st_end[12,2]],
                 paste(scan(text=DAT[st_end$sec_st[13]-1])[1:FAC_yrs],collapse=" "),
                 DAT[st_end[13,1]:st_end[13,2]],
                 paste(scan(text=DAT[st_end$sec_st[14]-1])[1:FAC_yrs],collapse=" "),
                 DAT[st_end[14,1]:st_end[14,2]])

    for(i in 1:FAC_yrs){DAT_retro = c(DAT_retro,
                                      paste(scan(text=DAT[st_end$sec_st[15] - FAC_yrs - 1 + i]), collapse = " "))
    }

    DAT_retro = c(DAT_retro,
                  DAT[st_end[15,1]:st_end[15,2]])

    # Survey age comp
    SAC_yrs = length(which(scan(text=DAT[st_end$sec_st[17] - 1]) <= (endyr - 1)))
    DAT_retro<-c(DAT_retro,
                 as.character(SAC_yrs),
                 DAT[st_end[16,1]:st_end[16,2]],
                 paste(scan(text = DAT[st_end$sec_st[17] - 1])[1:SAC_yrs], collapse=" "),
                 DAT[st_end[17,1]:st_end[17,2]],
                 paste(scan(text = DAT[st_end$sec_st[18] - 1])[1:SAC_yrs], collapse=" "),
                 DAT[st_end[18,1]:st_end[18,2]],
                 paste(scan(text = DAT[st_end$sec_st[19] - 1])[1:SAC_yrs], collapse=" "),
                 DAT[st_end[19,1]:st_end[19,2]],
                 paste(scan(text = DAT[st_end$sec_st[20] - 1])[1:SAC_yrs], collapse=" "),
                 DAT[st_end[20,1]:st_end[20,2]])
    for(i in 1:SAC_yrs) {
      DAT_retro = c(DAT_retro, paste(scan(text = DAT[st_end$sec_st[21] - SAC_yrs - 1 + i]), collapse = " "))
    }
    DAT_retro = c(DAT_retro, DAT[st_end[21,1]:st_end[21,2]])

    # Fish size comp
    FSC_yrs = length(which(scan(text = DAT[st_end$sec_st[23] - 1]) <= (endyr - 1)))

    DAT_retro<-c(DAT_retro,
                 as.character(FSC_yrs),
                 DAT[st_end[22,1]:st_end[22,2]],
                 paste(scan(text = DAT[st_end$sec_st[23] - 1])[1:FSC_yrs], collapse=" "),
                 DAT[st_end[23,1]:st_end[23,2]],
                 paste(scan(text = DAT[st_end$sec_st[24] - 1])[1:FSC_yrs], collapse=" "),
                 DAT[st_end[24,1]:st_end[24,2]],
                 paste(scan(text = DAT[st_end$sec_st[25] - 1])[1:FSC_yrs], collapse=" "),
                 DAT[st_end[25,1]:st_end[25,2]],
                 paste(scan(text = DAT[st_end$sec_st[26] - 1])[1:FSC_yrs], collapse=" "),
                 DAT[st_end[26,1]:st_end[26,2]])
    for(i in 1:FSC_yrs){
      DAT_retro = c(DAT_retro,
                    paste(scan(text = DAT[st_end$sec_st[27] - FSC_yrs - 1 + i]), collapse = " "))
    }
    DAT_retro = c(DAT_retro,
                  DAT[st_end[27,1]:st_end[27,2]])

    # Survey size comp
    SSC_yrs = length(which(scan(text = DAT[st_end$sec_st[29] - 1]) <= endyr))
    DAT_retro<-c(DAT_retro,
                 as.character(SSC_yrs),
                 DAT[st_end[28,1]:st_end[28,2]],
                 paste(scan(text = DAT[st_end$sec_st[29] - 1])[1:SSC_yrs], collapse=" "),
                 DAT[st_end[29,1]:st_end[29,2]],
                 paste(scan(text = DAT[st_end$sec_st[30] - 1])[1:SSC_yrs], collapse=" "),
                 DAT[st_end[30,1]:st_end[30,2]],
                 paste(scan(text = DAT[st_end$sec_st[31] - 1])[1:SSC_yrs], collapse=" "),
                 DAT[st_end[31,1]:st_end[31,2]],
                 paste(scan(text = DAT[st_end$sec_st[32] - 1])[1:SSC_yrs], collapse=" "),
                 DAT[st_end[32,1]:st_end[32,2]])
    for(i in 1:SSC_yrs)
    {DAT_retro = c(DAT_retro,
                   paste(scan(text = DAT[st_end$sec_st[33] - SSC_yrs - 1 + i]), collapse = " "))
    }
    DAT_retro = c(DAT_retro,
                  DAT[st_end[33,1]:st_end[33,2]])

    # Write data and control file
    write.table(DAT_retro,
                file = here::here(year, model, "retro", "model", paste0("goa_nr_", endyr, ".dat")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    CTL_retro = as.matrix(CTL)
    CTL_retro[2,1] = paste0("goa_nr_", endyr, ".dat")
    CTL_retro[5,1] = as.character(endyr)

    #Updated to account for fact that .tpl is looking for current model year
    write.table(CTL_retro,
                file = here::here(year, model, "retro", "model", paste0("goa_nr_", year, ".ctl")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    # run retro models

    ## set your number of MCMC runs at the top of the program...
    setwd(here::here(year, model, "retro", "model"))

    #Determine Operating system

    system(paste0(tpl_name,'.exe', ' -mcmc ', mcmc, ' -mcsave ', mcsave))
    system(paste0(tpl_name,'.exe',' -mceval'))

    file.copy("evalout.prj",
              here::here(year, model, "retro", "results", paste0("mcmc_", endyr, ".std")), overwrite = TRUE)

    file.copy(paste0(tpl_name, ".std"),
              here::here(year, model, "retro", "results", paste0("std_", endyr, ".std")), overwrite = TRUE)

    file.copy(here::here(year, model, "retro", "model", paste0(tpl_name, ".rep")),
              here::here(year, model, "retro", "results", paste0("rep_", endyr, ".rep")), overwrite = TRUE)

    # # sometimes I get weird output from admb it drops the model name and produces these files

    file.copy(here::here(year, model, "retro", "model", "update~1.rep"),
              here::here(year, model, "retro", "results", paste0("rep_", endyr, ".rep")), overwrite = TRUE)

    file.copy("update~1.std",
              here::here(year, model, "retro", "results", paste0("std_", endyr, ".std")), overwrite = TRUE)

  }
}




plot_compare_survey(year, models = c('2022, db', '2022, db.1', '2022, db.2', '2022, db.3', '2022, m22'))
plot_compare_biomass(year, models = c('2022, db', '2022, db.1', '2022, db.2', '2022, db.3', '2022, m22'))

concat_dat <- function(year, species, area = "goa", model="db", dat_name="goa_nr_2022", tsb_id =NULL, rec_age, plus_age, spawn_mo = 5,
                       maturity = NULL, n_ageage = 1, n_sizeage = 1, retro = NULL, ryear = NULL, n_fleets = 1, n_ts = NULL, n_lls = NULL){

  # create directory
  if (!dir.exists(here::here(year, model))){
    dir.create(here::here(year, model), recursive=TRUE)
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

  fishery = grep("fsh", list.files(here::here(year, 'data', "output")), value=TRUE)
  survey = grep("ts_", list.files(here::here(year, 'data', "output")), value=TRUE)
  ll_survey = grep("lls_", list.files(here::here(year, 'data', "output")), value=TRUE)

  catch = read.csv(here::here(year, "data", "output", grep("catch", fishery, value=TRUE)))
  waa = read.csv(here::here(year, "data", "output", "waa.csv"))
  saa = read.csv(here::here(year, "data", "output", "saa.csv"))
  ae = read.csv(here::here(year, "data", "output", "ae_model.csv"))
  fishac = read.csv(here::here(year, "data", "output", grep("age", fishery, value=TRUE)))
  fishlc = read.csv(here::here(year, "data", "output", grep("length", fishery, value=TRUE)))
  tsac = read.csv(here::here(year, "data", "output", grep("age", survey, value=TRUE)))
  tslc = read.csv(here::here(year, "data", "output", grep("length", survey, value=TRUE)))
  if(is.null(tsb_id)) {
    tsb = read.csv(here::here(year, "data", "output",
                              grep("biomass", survey, value=TRUE)))
  } else {
    tsb = read.csv(here::here(year, "data", "output",
                              grep(paste0("goa_ts_biomass_", tsb_id, ".csv"), survey, value=TRUE)))
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
              as.character(min(catch$Year)),
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
              as.character(min(catch$Year)),
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
                    paste0("#! ", paste(min(catch$Year):year, collapse=" ")),
                    paste(catch$Catch, collapse=" "),
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
    write.table(dat, file = here::here(year, model, paste0(dat_name, ".dat")) ,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  } else {
    write.table(dat, file = here::here(year, model, "retro", "model", ryear, paste0(dat_name, ".dat")) ,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}
