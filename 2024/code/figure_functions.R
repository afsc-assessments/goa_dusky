## Figure Functions
## October 2024
## K. Omori
# plot_compare_survey_dusk
# plot_compare_biomass_dusk
# plot_compare_biom_F_dusk
# plot_catch_dusk
# plot_catch_pres
# plot_catch_rate
# plot_recr_devs_dusk
# plot_retro_dusk
# plot_final_biomass
# plot_params_dusk (histograms of param est from mcmc)
# apport_catch_plot
# pearson_resid_plot

plot_compare_survey_dusk <- function(ayr, final= TRUE, models = c("m22.5a", "m22.3a_base", "base"),
                                     add_rema = TRUE, legend_x = .15, legend_y = .80, save = TRUE) {

  if (!dir.exists(here::here(ayr, "compare_models"))){
    dir.create(here::here(ayr, "compare_models"))
  }

  colors_mod <- c("darkgoldenrod1", "green4", "blue4")[c(1:length(models))]
  color_obs <- c("grey20", "grey50")
  mod_order <- c(models[1:(length(models)-1)], paste0(ayr-2) )

  dat = data.frame()

  for(i in 1:length(models)) {

    mod_name <- models[i]

    dat %>%
      dplyr::bind_rows(
        read.csv(here::here(ayr,  mod_name, "processed", "survey.csv")) %>%
          dplyr::rename_all(tolower) %>%
          dplyr::select(year = starts_with("y"),
                        Observed = starts_with("bio"),
                        Predicted = pred,
                        se, lci, uci) %>%
          tidyr::pivot_longer(-c(year, se, uci, lci)) %>%
          dplyr::mutate(value = value / 1000,
                        uci = uci / 1000,
                        lci = lci / 1000,
                        model = ifelse(mod_name == "base", paste0(ayr-2), mod_name)) ) -> dat
  }


  year_plot <- unique(dat$year)
  dat_new <- dat %>%
    dplyr::filter(name == "Predicted") %>%
    dplyr::select(year, name, value, model, uci, lci) %>%
    dplyr::mutate(model = factor(model, levels = c(rev(mod_order)))) %>%
    dplyr::bind_rows( dat %>%
                        dplyr::filter(name == "Observed",
                                      model %in% c(models[1])) %>%
                        dplyr::mutate(model = "observed (VAST)")
    )

  if(isTRUE(add_rema)) {
    rema_dat <- read.csv(here::here(ayr, "data", "output", "goa_biom_rema.csv")) %>%
      dplyr::mutate(model = "observed (db-rema)",
                    name = "Observed",
                    pred = pred/1000,
                    pred_lci = pred_lci/1000,
                    pred_uci = pred_uci/1000) %>%
      dplyr::select(year, name, value = pred, model, lci = pred_lci, uci = pred_uci)

    dat_new <- dplyr::bind_rows(dat_new, rema_dat)
  }

  if(isTRUE(add_rema)) {
    color_obs <- color_obs
  }
  if(isFALSE(add_rema)) {
    color_obs <- color_obs[1]
  }
  surv_plot1 <- ggplot2::ggplot() + #data = dat_new, ggplot2::aes(year, value, color = model)
    # Model points
    ggplot2::geom_point(data = dplyr::filter(dat_new, name == "Observed" ),
                        aes(x= year, y= value, color= model) , alpha = 0.8,
                        position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_errorbar(data = dplyr::filter(dat_new, name == "Observed" ),
                           ggplot2::aes(x= year, ymin = lci, ymax = uci, color = model, linetype = model), width = 0.5, alpha = 0.8,
                           position=ggplot2::position_dodge(width=0.7)) +
    scale_color_manual(name = "Observed", values = c(color_obs[2], color_obs[1])) +
    scale_linetype_manual(name = "Observed", values = c(2, 1)) +
    # adding predicted survey
    ggnewscale::new_scale_color() +
    ggplot2::geom_line(data = dplyr::filter(dat_new, name == "Predicted"),
                       ggplot2::aes(x= year, y= value, color = model), linewidth= 1) +
    ggplot2::geom_point(data = dplyr::filter(dat_new, name == "Predicted"), size = 1.5,
                        ggplot2::aes(x= year, y= value, color = model)) +
    scale_color_manual(name = "Model", values = colors_mod , guide = "none") +

    #scico::scale_color_scico_d(palette = 'navia', begin = 0.2, end = 0.8) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    theme_bw() +
    #funcr::theme_report() +
    ggplot2::theme(legend.position = c(legend_x, legend_y))


  surv_plot2 <- ggplot2::ggplot() + #data = dat_new, ggplot2::aes(year, value, color = model)
    ggplot2::geom_point(data = dplyr::filter(dat_new, name == "Observed" & model == "observed (VAST)" ),
                        aes(x= year, y= value, color= "Observed (VAST)") , alpha = 0.8, color = color_obs[1],
                        #ggplot2::aes(color = "Observed (VAST)") , alpha = 0.8,
                        #ggplot2::aes(color = model) , alpha = 0.8,
                        position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_errorbar(data = dplyr::filter(dat_new, name == "Observed" & model == "observed (VAST)" ),
                           ggplot2::aes(x= year, ymin = lci, ymax = uci, color = "Observed (VAST)"), width = 0.5, alpha = 0.8,
                           color = color_obs[1], linetype = 1,
                           #ggplot2::aes(ymin = lci, ymax = uci, color = "Observed (VAST)"), width = 0.5, alpha = 0.8,
                           #ggplot2::aes(ymin= lci, ymax = uci, color = model), width = 0.5, alpha = 0.8,
                           position=ggplot2::position_dodge(width=0.7)) +
    scale_linetype_manual(name = "Observed", values = c(2, 1) ) +
    ggplot2::geom_line(data = dplyr::filter(dat_new, name == "Predicted"),
                       ggplot2::aes(x= year, y= value, color = model), linewidth= 1) +
    ggplot2::geom_point(data = dplyr::filter(dat_new, name == "Predicted"), size = 1.5,
                        ggplot2::aes(x= year, y= value, color = model)) +
    scale_color_manual(name = "Model", values = colors_mod ) +
    #scico::scale_color_scico_d(palette = 'navia', begin = 0.2, end = 0.8) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    theme_bw() +
    #funcr::theme_report() +
    ggplot2::theme(legend.position = c(legend_x, legend_y))
  #surv_plot2

  survey_fig <- ggpubr::ggarrange(surv_plot1, surv_plot2,
                                      labels= c("A", "B"),
                                      nrow= 2)

  if(isFALSE(add_rema)) {
    pres_surv <- surv_plot2 +
      ggplot2::theme(legend.justification=c(1,0),
                     legend.position=c(0.25,0.68),
                     axis.text=element_text(size=10),
                     axis.title=element_text(size=12,face="bold"),
                     legend.text=element_text(size=10))

    ggplot2::ggsave(plot = pres_surv,
                    here::here(ayr, "compare_models", "survey_pres.png"),
                    width = 6, height = 4, units = "in", dpi = 200)
  }


  if(isTRUE(save)) {
    if(isTRUE(final)) {
      ggplot2::ggsave(plot = survey_fig,
                      here::here(ayr, "compare_models", "compare_est_survey.png"),
                      width = 6.5, height = 7, units = "in", dpi = 200)
    }
    if(isFALSE(final)) {
      ggplot2::ggsave(plot = survey_fig,
                      here::here(ayr, "compare_models", "compare_est_survey.png"),
                      width = 6.5, height = 6.5, units = "in", dpi = 200)
    }
  }

}

print("plot_compare_survey_dusk")



plot_compare_biom_F_dusk <- function(year, pref_mod, base_current,   base_pr_mod = "base", admb_name = "base") {

  ayr <- year
  pref_mod_std <- read.delim(here::here(year, pref_mod, paste0(admb_name, ".std")), sep="", header = TRUE) %>%
    dplyr::mutate(model = pref_mod)
  base_new_mod_dat <- read.delim(here::here(year, base_current, paste0(admb_name, ".std")), sep="", header = TRUE) %>%
    dplyr::mutate(model = base_current )
  base_old_mod_dat <- read.delim(here::here(year, base_pr_mod, paste0(admb_name, ".std")), sep="", header = TRUE) %>%
    dplyr::mutate(model = as.character(paste0(ayr-2, "_base")) )

  ### SOMETHING IS INCORRECT with the LCI AND UCI for plots F
  ### WILL NEED TO CHANGE THE F DEVS INTO EXP FOR THE VALUE AND FOR CI
  F_dat_base_old <- read.csv(here::here(year, base_pr_mod, "tables", "t_recr_biom_mcmc.csv")) %>%
    dplyr::select(Year = year, value = Frate, lci = Frate_lci, uci = Frate_uci) %>%
    dplyr::mutate(model = as.character(paste0(ayr-2, "_base"))) %>%
    dplyr::filter(!is.na(value))
  F_dat_pref_mod <- read.csv(here::here(year, pref_mod, "tables", "t_recr_biom_mcmc.csv")) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::select(Year = year, value = frate_mean, lci = frate_lci, uci = frate_uci) %>%
    dplyr::mutate(model = as.character(paste0(ayr-2, "_base")))
  F_dat_base_new <- read.csv(here::here(year, base_current, "processed", "bio_rec_f.csv")) %>%
    dplyr::select(Year = year, value = F) %>%
    dplyr::mutate(model = base_current)

  full_f_dat <- dplyr::bind_rows(F_dat_pref_mod, F_dat_base_new, F_dat_base_old) %>%
    dplyr::filter(Year <= ayr)

  mod_dat <- dplyr::bind_rows(pref_mod_std, base_new_mod_dat, base_old_mod_dat)

  sb_dat <- data.frame()
  tb_dat <- data.frame()
  recruit_dat <- data.frame()
  #full_f_dat <- data.frame()
  model_all <- c(pref_mod , base_current, as.character(paste0(ayr-2, "_base")))

  for(m in 1:length(model_all)) {

    dat_tmp <- mod_dat %>%
      dplyr::filter(model %in% c(model_all[m] ))

    if(m<3) {years_vec <- (1977:ayr) }
    if(m==3) {years_vec <- 1977:(ayr-2)}

    sb_dat <- dat_tmp %>%
      dplyr::filter(name == "spawn_biom") %>%
      dplyr::mutate(Year = years_vec) %>%
      dplyr::bind_rows(sb_dat)
    tb_dat <- dat_tmp %>%
      dplyr::filter(name == "tot_biom") %>%
      dplyr::mutate(Year = years_vec) %>%
      dplyr::bind_rows(tb_dat)
    recruit_dat <- dat_tmp %>%
      dplyr::filter(name == "pred_rec") %>%
      dplyr::mutate(Year = years_vec) %>%
      dplyr::bind_rows(recruit_dat)
    log_avg_f <- dat_tmp %>%
      dplyr::filter(name == "log_avg_F")
    #full_f_dat <- dat_tmp %>%
    #  dplyr::filter(name == "log_F_devs") %>%
    #  dplyr::mutate(value = exp(value+log_avg_f$value),
    #                std.dev = std.dev  ) %>%
    #  dplyr::mutate(Year = years_vec) %>%
    #  dplyr::bind_rows(full_f_dat)

  } # end model loop

  # Individual plots

  ssb_fig <- ggplot(sb_dat, aes(x= Year, y= value/1000, color= model, fill= model)) +
    geom_line(linewidth= 1.15, aes(color = model)) +
    geom_ribbon(aes(ymin= (value-std.dev)/1000, ymax= (value+std.dev)/1000 ) , alpha= 0.1, linetype=0) +
    scale_color_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
    scale_fill_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
    labs(title = "Spawning biomass" , x = NULL, y = "Biomass (kt)") +
    expand_limits(y=0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5) )

  ggsave(here::here( year, 'compare_models', 'fig_ssb_mod_compare.png'), plot= ssb_fig,
         width = 7, height = 4)

  tb_fig <- ggplot(tb_dat, aes(x= Year, y= value/1000, color= model, fill= model)) +
    geom_line(linewidth= 1.15, aes(color = model)) +
    geom_ribbon(aes(ymin= (value-std.dev)/1000, ymax= (value+std.dev)/1000 ) , alpha= 0.1, linetype=0) +
    scale_color_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
    scale_fill_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
    labs(title = "Total biomass" , x = NULL, y = "Biomass (kt)") +
    expand_limits(y=0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5) )

  ggsave(here::here( year, 'compare_models', 'fig_tb_mod_compare.png'), plot= ssb_fig,
         width = 7, height = 4)

  recruit_fig <- ggplot(recruit_dat, aes(x= Year, y= value, color= model, fill= model)) +
    geom_line(linewidth= 1.15, aes(color = model)) +
    geom_ribbon(aes(ymin= (value-std.dev), ymax= (value+std.dev) ) , alpha= 0.1, linetype=0) +
    scale_color_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
    scale_fill_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
    labs(title = "Age-4 Recruits" , x = "Year", y = "Numbers (millions)") +
    expand_limits(y=0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5) )

  ggsave(here::here( year, 'compare_models', 'fig_recruit_mod_compare.png'), plot= recruit_fig,
         width = 7, height = 4)

  f_fig <- ggplot(full_f_dat, aes(x= Year, y= value, color= model, fill= model)) +
    #geom_line(linewidth= 1.15, aes(color = model)) +
    #geom_ribbon(aes(ymin= (value-std.dev), ymax= (value+std.dev) ) , alpha= 0.1, linetype=0) +
    geom_ribbon(aes(ymin= lci, ymax= uci ) , alpha= 0.1, linetype=0) +
    geom_line(linewidth= 1.15, aes(color = model)) +
    scale_color_manual(values = c("blue4", "green4")) +
    scale_fill_manual(values = c("blue4", "green4")) +
    expand_limits(y=0) +
    labs(title = "Fishing Mortality" , x = NULL, y = "F") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5) )

  ggsave(here::here( year, 'compare_models', 'fig_f_mod_compare.png'), plot= f_fig,
         width = 7, height = 4)

  # Combine graphics into one panel
  ts_compare_fig <- ggpubr::ggarrange(ssb_fig, tb_fig, f_fig, recruit_fig,
                                      labels= c("A", "B", "C", "D"),
                                      common.legend = TRUE,
                                      legend= "right",
                                      nrow= 4)
  ggsave(here::here( year, 'compare_models', 'fig_ts_mod_compare.png'), plot= ts_compare_fig,
         width = 7, height = 9)

}

# catch
#' Plot catch
#'
#' @param year model year
#' @param folder folder name model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_catch
#'
#' @examples
plot_catch_dusk <- function(year, folder, save=TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs
  vroom::vroom(here::here(year, folder, "processed", "catch.csv")) %>%
    tidytable::bind_cols(year = yrs) %>%
    dplyr::mutate(Observed = obs/1000 ,
                  Estimated = pred/ 1000 ,
                  years = "All years") %>%
    dplyr::select(year, Observed, Estimated, years) -> dat


  tidytable::filter(dat, year %in% (max(yrs) - 20):max(yrs)) %>%
    tidytable::mutate(years = "Recent years") %>%
    tidytable::bind_rows(dat) %>%
    tidytable::pivot_longer(c(-year, -years)) %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name, lty = name)) +
    ggplot2::geom_line() +
    scico::scale_color_scico_d(name = "", palette = "roma") +
    ggplot2::scale_linetype_manual(name = "",
                                   values = c(1,1)) +
    ggplot2::facet_wrap(~years, scales = "free",
                        dir = "v") +
    ggplot2::ylab("Catch (kt)") +
    ggplot2::xlab("Year") +
    ggplot2::expand_limits(y = 0) +
    afscassess::scale_x_tickr(data=dat, var=year, start=1960) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.15,0.9)) -> plot_catch

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "catch.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }

}

plot_catch_pres <- function(year, folder, save=TRUE){

  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results' before creating figures")
  }
  # set view
  #ggplot2::theme_set(afscassess::theme_report())

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs
  vroom::vroom(here::here(year, folder, "processed", "catch.csv")) %>%
    tidytable::bind_cols(year = yrs) %>%
    dplyr::mutate(Observed = obs/1000 ,
                  Estimated = pred/ 1000 ,
                  years = "All years") %>%
    dplyr::select(year, Observed, Estimated, years) %>%
    tidyr::pivot_longer(cols = c(Observed, Estimated), names_to = "name", values_to = "value") -> dat

catch_plot <-  ggplot2::ggplot(dat, ggplot2::aes(year, value, color = name, lty = name)) +
    ggplot2::geom_line(linewidth = 2) +
    scico::scale_color_scico_d(name = "", palette = "roma") +
    ggplot2::scale_linetype_manual(name = "",
                                   values = c(1,1)) +
    ggplot2::ylab("Catch (kt)") +
    ggplot2::xlab("Year") +
    ggplot2::expand_limits(y = 0) +
    #afscassess::scale_x_tickr(data=dat, var=year, start=1960) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.2,0.8),
                   axis.text=element_text(size=12),
                   axis.title=element_text(size=14,face="bold"),
                   legend.text=element_text(size=12))

    ggplot2::ggsave(here::here(year, folder, "figs", "catch_pres.png"),
                    width = 7, height = 5, units = "in", dpi = 200)


}

## Plot Catch/ Biomass (catch rate)

plot_catch_rate <- function(year, admb_name = "base", model) {
  ## from 2023_analysis need to turn into function
  # plot catch/biomass
  ayr = year

  std <- read.delim(here::here(year, model, paste0(admb_name, ".std")), sep="", header = TRUE)
  catch <- vroom::vroom(here::here(year, "data", "output", "fish_catch.csv"))

  filter(catch, year == max(ayr)) %>%
    left_join(vroom::vroom(here::here(year, model, "proj", "author_f", "bigsum.csv")) %>%
                filter(Year == year, Alt == 2) %>%
                dplyr::select(year=Year, value = Total_Biom)) -> final_yr_catch

  std %>%
    filter(name=="tot_biom") %>%
    #filter(row_number() < n() ) %>% # removing last row for "current year"
    #mutate(year = 1977:(ayr-1) ) %>%
    bind_cols(filter(catch, year <= max(ayr))) %>%
    filter(year >= 1991) %>%
    dplyr::select(year, catch, value, std.dev) %>%
    bind_rows(filter(catch, year == max(ayr)) %>%
                left_join(vroom::vroom(here::here(year, model, "proj", "author_f", "bigsum.csv")) %>%
                            filter(Year == year, Alt == 2) %>%
                            mutate(value = Total_Biom * 1000) %>%
                            dplyr::select(year=Year, value))) %>%
    mutate(std.dev = ifelse(is.na(std.dev), std.dev[year==max(ayr)-1], std.dev)) %>%
    mutate(lci = value - std.dev * 1.96,
           uci = value + std.dev * 1.96) %>%
    mutate(perc = catch / value,
           lci = catch / lci,
           uci = catch / uci,
           mean = mean(perc)) %>%
    dplyr::select(year, value, mean, perc, lci, uci) -> df_tmp

   catch_rate_plot <-  ggplot2::ggplot(df_tmp, ggplot2::aes(year, perc)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size= 1.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.2) +
    ggplot2::geom_hline(yintercept = df_tmp$mean, lty = 3) +
    ggplot2::expand_limits(y = c(0, 0.08)) +
    afscassess::scale_x_tickr(data=df_tmp, var=year, start = 1990) +
    afscassess::theme_report() +
    ggplot2::xlab("\nYear") +
    ggplot2::ylab("Catch/Biomass\n")

  ggsave(here::here( year, pref_mod, 'figs' , 'f_catchrate.png'), plot= catch_rate_plot,
         width = 6, height = 4)

}

## Plot of recruitment deviations ----

plot_recr_devs_dusk <- function(year, pref_mod, save = TRUE) {

  recr_dev_mcmc <- read.csv(here::here(year, pref_mod, "mcmc", "processed", "mceval.csv")) %>%
    dplyr::select(contains("log_rec_dev")) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(-id, names_to = "Parameter", values_to = "value") %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(mean_mcmc = mean(value),
                     median_mcmc = median(value),
                     sd_mcmc = sd(value),
                     lci_mcmc = quantile(value, 0.025),
                     uci_mcmc = quantile(value, 0.975)) %>%
    dplyr::mutate(Year = as.integer(str_extract(Parameter, "[0-9]+")),
                  model = pref_mod)

  plot_log_rec_devs <- ggplot(recr_dev_mcmc, aes(x= Year, y= median_mcmc)) +
    geom_linerange(aes(x= Year, ymin=lci_mcmc, ymax=uci_mcmc), linewidth = 0.5, colour = "grey55") +
    geom_point(size = 1.5) +
    geom_hline(yintercept= 0, color= "grey50", linetype = 2) +
    labs(x = "Year", y = "Log recruit deviations")

  if(isTRUE(save)) {
    ggsave(here::here( year, pref_mod, "mcmc", 'figs' , 'fig_recruit_devs.png'), plot= plot_log_rec_devs,
           width = 6, height = 4)
  }

} # end plot recruit devs

## Plot retrospective plots (and save Mohn's rho) ----

plot_retro_dusk <- function(year, folder, n_retro=10, save_indiv=TRUE, rhos = c("ssb", "totbiom")) {

  peels = n_retro - 1
  max_year = year
  # loop through mcmc output
  age_yr = read.csv(here::here(year, folder, "processed", "ages_yrs.csv"))
  yrs = age_yr %>%
    dplyr::select(yrs) %>%
    tidyr::drop_na() %>%
    dplyr::pull(yrs)
  styr_rec = age_yr[1,3]
  retro_yrs = (year - n_retro + 1):year

  dat = list()

  for(i in 1:n_retro) {

    read.delim(here::here(year, folder, "retro", "results",
                          paste0("mcmc_", retro_yrs[i], ".std")),
               sep = "",  header = FALSE) -> df

    df = df[(0.2 * nrow(df)):nrow(df),] # drop burn in

    colnames(df) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
                     "ABC", "obj_fun",                                                    # 8
                     paste0("tot_biom_", yrs[1]:retro_yrs[i]),                            #39
                     paste0("log_rec_dev_", (styr_rec-2):retro_yrs[i]),  # hardcoded because styr_rec is incorrect
                     paste0("spawn_biom_", yrs[1]:retro_yrs[i]),
                     "log_mean_rec",
                     paste0("spawn_biom_proj_", max(retro_yrs[i]) + 1:15),
                     paste0("pred_catch_proj_", max(retro_yrs[i]) + 1:15),
                     paste0("rec_proj_", max(retro_yrs[i]) + 1:10),
                     paste0("tot_biom_proj_", max(retro_yrs[i]) + 1))

    dat[[i]] = df %>% dplyr::mutate(retro_year = retro_yrs[i])

  }

  # clean up columns so can bind all together
  col_name_tmp <- unique(unlist(sapply(dat, names)))
  dat <- lapply(dat, function(df) {
    df[, setdiff(col_name_tmp, names(df))] <- NA
    df
  })

  do.call(rbind, dat)  -> retro_mc

  # save output
  write.csv(retro_mc, here::here(year, folder, "processed", "retro_mcmc.csv"), row.names = FALSE)

  #retro_mc <- read.csv(here::here(year, folder, "processed", "retro_mcmc.csv"))

  # functions for quantiles
  q_name <- purrr::map_chr(c(.025,.975), ~ paste0("q", .x*100))
  q_fun <- purrr::map(c(.025,.975), ~ purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>%
    purrr::set_names(nm = q_name)

  # rho and plots for spawning biomass and total biomass
  retro_mc %>%
    dplyr::select(paste0("spawn_biom_", yrs), paste0("tot_biom_", yrs), retro_year) %>%
    tidyr::pivot_longer(c(-retro_year), values_to = "biomass") %>%
    dplyr::mutate(biom_var = paste0(stringr::str_extract(name, "[a-z]+"), "_biom" ),
                  year = as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  biomass = biomass / 1000) %>%
    dplyr::group_by(biom_var, year, retro_year) %>%
    dplyr::summarise_at(dplyr::vars(biomass), tibble::lst(!!!q_fun, median)) %>%
    dplyr::mutate(Retro = factor(retro_year)) %>%
    dplyr::ungroup() -> dat_retro_quant


  dat_retro_quant %>%
    dplyr::select(biom_var, year, retro_year, median) %>%
    dplyr::group_by(biom_var, year) %>%
    dplyr::mutate(pdiff = (median - median[retro_year==max_year]) /
                    median[retro_year==max_year]) %>%
    tidyr::drop_na() %>%
    dplyr::filter(year %in% (max_year-peels):max_year) %>%
    dplyr::ungroup() %>%
    dplyr::filter(year == retro_year, year !=max_year) %>%
    dplyr::group_by(biom_var) %>%
    dplyr::summarise(rho = mean(pdiff)) -> ssb_tb_rho

  write.csv(ssb_tb_rho, here::here(year, folder, "processed", "rho_ssb_tb.csv"), row.names = FALSE)

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  ssb_plot <- ggplot(dat_retro_quant %>% dplyr::filter(biom_var == "spawn_biom"),
                     aes(x= year, y= median, color = Retro, group = Retro)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5, fill = Retro),
                         alpha = .1, color = NA) +
    ggplot2::ylab("Spawning biomass (kt)\n") +
    ggplot2::xlab("\nYear") +
    ggplot2::expand_limits(y = 0) +
    scico::scale_fill_scico_d(palette = "roma") +
    scico::scale_color_scico_d(palette = "roma") + #'batlow
    afscassess::theme_report() +
    ggplot2::scale_x_continuous(breaks = afscassess::tickr(dat_retro_quant, year, 10, start = 1977)$breaks,
                                labels = afscassess::tickr(dat_retro_quant, year, 10, start = 1977)$labels) +
    ggplot2::annotate(geom = "text", x=1978, y=Inf, hjust = -0.05, vjust = 2,
                      label = paste0("Mohn's rho = ",
                                     round(ssb_tb_rho %>% filter(biom_var == "spawn_biom") %>% pull(rho), 3)),
                      family = "Times") +
    ggplot2::theme(legend.position = "none")

  tb_plot <- ggplot(dat_retro_quant %>% dplyr::filter(biom_var == "tot_biom"),
                     aes(x= year, y= median, color = Retro, group = Retro)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5, fill = Retro),
                         alpha = .1, color = NA) +
    ggplot2::ylab("Total biomass (kt)\n") +
    ggplot2::xlab("\nYear") +
    ggplot2::ylim(0,250) +
    #ggplot2::expand_limits(y = 0 ) +
    scico::scale_fill_scico_d(palette = "roma") +
    scico::scale_color_scico_d(palette = "roma") + #'batlow
    afscassess::theme_report() +
    ggplot2::scale_x_continuous(breaks = afscassess::tickr(dat_retro_quant, year, 10, start = 1977)$breaks,
                                labels = afscassess::tickr(dat_retro_quant, year, 10, start = 1977)$labels) +
    ggplot2::annotate(geom = "text", x=1978, y=Inf, hjust = -0.05, vjust = 2,
                      label = paste0("Mohn's rho = ",
                                     round(ssb_tb_rho %>% filter(biom_var == "tot_biom") %>% pull(rho), 3)),
                      family = "Times") #+
    ggplot2::theme(legend.position = "none")

  if(isTRUE(save_indiv)) {

    ggplot2::ggsave(plot = ssb_plot, here::here(year, folder, "figs", "retro_ssb.png"),
                    width = 6.5, height = 4, units = "in", dpi = 200)
    ggplot2::ggsave(plot = tb_plot, here::here(year, folder, "figs", "retro_tb.png"),
                    width = 6.5, height = 4, units = "in", dpi = 200)

    ssb_plot2 <- ssb_plot +
      theme(axis.title.x = element_blank())

    biom_retro <- cowplot::plot_grid(ssb_plot2, tb_plot, ncol = 1)
    ggplot2::ggsave(plot = biom_retro, here::here(year, folder, "figs", "retro.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)

  }
} # end retro plots

# Final SSB & total biomass plot

plot_final_biomass <- function(year, model, save = TRUE ) {

  ayr = year

  dat_tmp <- read.csv(here::here(year, model, "tables", "t_recr_biom_mcmc.csv")) %>%
    dplyr::filter(Year <= ayr ) %>%
    dplyr::select(Year, tot_mean, tot_lci, tot_uci, sp_mean, sp_lci, sp_uci) %>%
    tidyr::pivot_longer(cols = -Year, names_to= "name", values_to = "value") %>%
    tidyr::separate(name, into = c("biom_var", "type"), sep = "_") %>%
    dplyr::mutate(biom_var = case_when(biom_var == "tot"~ "Total Biomass",
                                   biom_var == "sp"~ "Spawning Biomass")) %>%
    tidyr::pivot_wider(names_from = type, values_from = value) %>%
    dplyr::rename(biom = mean)


  plot_final_tb <- ggplot(dat_tmp, aes(x= Year, y= biom/1000)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci/1000, ymax = uci/1000 ), alpha =  0.15) +
    ggplot2::ylab("Biomass (kt)\n") +
    ggplot2::xlab("\nYear") +
    expand_limits(y = 0) +
    ggplot2::facet_wrap(~biom_var, ncol = 1, scales = "free_y")

  ggsave(here::here( year, model, "figs", 'f_final_biom.png'), plot= plot_final_tb,
         width = 7, height = 6)
}

# Historgram of parameter estimates
' Parameter histogram plots
#'
#' @param year assessment year
#' @param folder   the folder with the model in it
#' @param model_name the name of your .tpl file
#' @param pars parameter names
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @return
#' @export plot_params
#'
#' @examples plot_swath(year, folder)
plot_params_dusk <- function(year, folder, model_name, pars = c("q_srv1", "ABC", "nattymort", "log_mean_rec", "tot_biom", "F40","spawn_biom"), save=TRUE) {
  if (!dir.exists(here::here(year, folder, "processed"))){
    stop("must run 'process_results()' before creating figures")
  }

  # set view
  ggplot2::theme_set(afscassess::theme_report())

  vroom::vroom(here::here(year, folder, "processed", "ages_yrs.csv"))$yrs -> yrs

  read.delim(here::here(year, folder, paste0(model_name, ".std")), sep="", header = TRUE) %>%
    tidytable::filter(name %in% pars) %>%
    tidytable::slice(tidytable::n(),
                     .by = name) %>%
    tidytable::mutate(name = tidytable::case_when(name =="q_srv1" ~ "q_srv",
                                                  name =="nattymort" ~ "natmort",
                                                  name == "log_mean_rec" ~ "log_mean_rec",
                                                  name == "F40" ~ "F",
                                                  name == "spawn_biom" ~ "spawn_biom_",
                                                  name == "tot_biom" ~ "tot_biom_",
                                                  name == "pred_recr" ~ "pred_recr",
                                                  TRUE ~ name),
                      value = tidytable::case_when(name == "ABC" ~ value / 1000,
                                                   name == "spawn_biom_" ~ value / 1000,
                                                   name == "tot_biom_" ~ value / 1000,
                                                   TRUE ~ value),
                      name = factor(name, levels = c("q_srv", "natmort", "log_mean_rec", "F", "ABC", "tot_biom_", "spawn_biom_"))) -> fits

  vroom::vroom(here::here(year, folder, "processed", "mceval.csv"))  %>%
    tidytable::select(q_srv1, ABC, natmort, log_mean_rec, paste0("tot_biom_", yrs),
                      F40, paste0("spawn_biom_", yrs)) %>%
    tidytable::mutate(group = 1:dplyr::n()) %>%
    tidytable::pivot_longer(-group) %>%
    tidytable::mutate(years = as.numeric(gsub('\\D+','', name)),
                      name = gsub('[[:digit:]]+', '', name),
                      value = tidytable::case_when(name=="spawn_biom_" ~ value / 1000,
                                                   name=="tot_biom_" ~ value / 1000,
                                                   name=="ABC" ~ value / 1000,
                                                   TRUE ~ value),
                      name = factor(name, levels = c("q_srv", "natmort", "log_mean_rec", "F", "ABC", "tot_biom_", "spawn_biom_"))) -> dat

  p1 = dat %>%
    tidytable::filter(name == "q_srv") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_vline(data = tidytable::filter(fits, name == "q_srv"), aes(xintercept = value),
                        linewidth = 1.5 , color = "black", linetype = 2) +
    #ggplot2::geom_segment(data = tidytable::filter(fits, name == "q_srv"),
    #                      mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
    #                      size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab(expression("Trawl survey catchability ("*italic(q)*")"))

  p2_b = dat %>%
    tidytable::filter(name == "natmort") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = tidytable::filter(fits, name == "natmort"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), size = 2, color = "black") +
    ggplot2::ylab("Probability density") +
    ggplot2::xlab(expression("Natural mortality ("*italic(M)*")")) +
    ggplot2::theme(legend.position = "none")
# do log_mean_rec
  p2 = dat %>%
    tidytable::filter(name == "log_mean_rec") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    #ggplot2::geom_segment(data = tidytable::filter(fits, name == "log_mean_rec"),
    #                      mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), size = 2, color = "black") +
    ggplot2::geom_vline(data = tidytable::filter(fits, name == "log_mean_rec"), aes(xintercept = value),
                        linewidth = 1.5 , color = "black", linetype = 2) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::ylab("Probability density") +
    ggplot2::xlab(expression("Log mean recruitment")) +
    ggplot2::theme(legend.position = "none")

  p3 = dat %>%
    tidytable::filter(name == "F") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_vline(data = tidytable::filter(fits, name == "F"), aes(xintercept = value),
                        linewidth = 1.5 , color = "black", linetype = 2) +
    #ggplot2::geom_segment(data = tidytable::filter(fits, name == "F"),
    #                      mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
    #                      size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab(expression(italic(F)["40%"])) +
    ggplot2::theme(legend.position = "none")

  p4 = dat %>%
    tidytable::filter(name == "ABC") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_vline(data = tidytable::filter(fits, name == "ABC"), aes(xintercept = value),
                        linewidth = 1.5 , color = "black", linetype = 2) +
    #ggplot2::geom_segment(data = tidytable::filter(fits, name == "ABC"),
    #                      mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
    #                      size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("ABC (kt)") +
    ggplot2::theme(legend.position = "none")


  p5 = dat %>%
    tidytable::filter(name == "tot_biom_", years == year) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_vline(data = tidytable::filter(fits, name == "tot_biom_"), aes(xintercept = value),
                        linewidth = 1.5 , color = "black", linetype = 2) +
    #ggplot2::geom_segment(data = tidytable::filter(fits, name == "tot_biom_"),
    #                      mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
    #                      size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("Current total biomass (kt)") +
    ggplot2::theme(legend.position = "none")

  p6 = dat %>%
    tidytable::filter(name == "spawn_biom_", years == year) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_vline(data = tidytable::filter(fits, name == "spawn_biom_"), aes(xintercept = value),
                        linewidth = 1.5 , color = "black", linetype = 2) +
    #ggplot2::geom_segment(data = tidytable::filter(fits, name == "spawn_biom_"),
    #                      mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05),
    #                      size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("Current spawning biomass (kt)") +
    ggplot2::theme(legend.position = "none")

  #print((cowplot::plot_grid(p1, p4, p2,  p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))
  # removing natural mort plot
  print((cowplot::plot_grid(p1, p4, p2, p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))

  if(isTRUE(save)) {
    png(filename=here::here(year, folder, "figs", "hists.png"), width = 6.5, height = 6.5,
        units = "in", type ="cairo", res = 200)

    print((cowplot::plot_grid(p1, p4, p2, p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))

    dev.off()

  }

} # end histogram of param est

## Apportionment: catch with subarea ABCs + new one
apport_catch_plot <- function(year, model, model_name= "vast", save_area= TRUE) {

  ayr <- year

  # fishery catch by area
  fish_catch_area <- read.csv(here::here(ayr, "data", "raw", "fish_catch_data.csv")) %>%
    dplyr::filter(agency_species_code %in% c(154, 172),
                  year %in% c(2012:ayr)) %>%
    dplyr::group_by(year, fmp_subarea) %>%
    dplyr::summarise(weight = sum(weight_posted, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Year = year, Area = fmp_subarea, Catch = weight) %>%
    dplyr::mutate(Area = case_when(Area == "WG" ~ "Western",
                                   Area == "CG" ~ "Central",
                                   Area == "WY" ~ "West Yakutat",
                                   Area == "SE" ~ "Southeast")) %>%
    dplyr::mutate(Area = factor(Area, levels = c("Western", "Central", "West Yakutat", "Southeast")))

  if(isTRUE(save_area)) {
    write.csv(fish_catch_area, here::here(year, "data", "output", "fish_catch_area.csv"), row.names = FALSE)
  }

  #fishery catch total goa
  fish_catch_goa <- fish_catch_area %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(Catch = sum(Catch, na.rm= T)) %>%
    dplyr::mutate(area = "GOA")

  # goa harvest limits
  goa_hl <- read.csv(here::here(ayr, "data", "user_input", "abc_tac_goa.csv")) %>%
    dplyr::filter(Year %in% c(2012:ayr)) %>%
    dplyr::select(Year, ABC, OFL) %>%
    tidyr::pivot_longer(cols = c("ABC", "OFL"), names_to = "mgmt", values_to = "value") %>%
    dplyr::mutate(harv_specs = "Historical")

  # read in the subarea ABC
  abc_subarea <- read.csv(here::here(year, "data", "user_input", "goa_dusk_harvest_specs.csv")) %>%
    dplyr::filter(Area == "ABC") %>%
    dplyr::select(Year, mgmt = Area, Western = W, Central = C, "West Yakutat" = WYK, Southeast = SEO) %>%
    tidyr::pivot_longer(cols = -c(Year, mgmt), names_to = "Area", values_to = "value") %>%
    dplyr::mutate(harv_specs = "Historical")

  # Proposed harvest recommendations

  proposed_goa_abc <- read.csv(here::here(year, model, "processed", "exec_summ.csv")) %>%
    dplyr::filter(item %in% c("ABC (t)", "OFL (t)"  ) ) %>%
    dplyr::select(mgmt = item, value = y3) %>%
    dplyr::mutate(harv_specs = "Proposed",
                  mgmt = case_when(mgmt == "ABC (t)" ~ "ABC",
                                   mgmt == "OFL (t)" ~ "OFL"),
                  Year = (ayr +1),
                  value = as.numeric(value))

  proposed_subarea <- read.csv(here::here(year, model, "apport", paste0("apport_prop_", model_name, ".csv")) ) %>%
    dplyr::slice( n() ) %>% # takes final year
    dplyr::mutate(wyak = wyak*eastern,
                  seo = seo * eastern) %>%
    dplyr::select(Western = western, Central = central,
                  "West Yakutat" = wyak, Southeast = seo, model_name) %>%
    tidyr::pivot_longer(cols = -model_name, names_to = "Area", values_to = "proportion" ) %>%
    dplyr::mutate(value = proportion * (proposed_goa_abc[proposed_goa_abc$mgmt == "ABC", 2]),
                  harv_specs = "Proposed",
                  Year = (ayr +1),
                  mgmt = "ABC")

  # Final data for figures
  # goa
  goa_hl_dat <- dplyr::full_join(goa_hl, proposed_goa_abc, by = join_by(Year, mgmt, harv_specs, value))

  # subarea
  subarea_hl_dat <- dplyr::full_join(abc_subarea, proposed_subarea, by = join_by(Year, mgmt, Area, harv_specs, value)) %>%
    dplyr::mutate(Area = factor(Area, levels = c("Western", "Central", "West Yakutat", "Southeast")))

  # GOA total ABC and OFL with catch
  goa_hl_fig <- ggplot() +
    geom_bar(data= fish_catch_goa, stat = "identity", alpha = 0.5, aes(x = Year, y = Catch), color = "grey") +
    geom_line(data=goa_hl_dat, aes(x= Year, y= value, color= mgmt, linetype = harv_specs), linewidth= 1) +
    geom_line(data=goa_hl_dat %>% filter(Year %in% c(ayr: (ayr+1))), aes(x = Year, y = value, color = mgmt), linetype = 2) +
    geom_point(data= goa_hl_dat %>% filter(Year %in% (ayr+1)), aes(x= Year, y= value, color = mgmt), shape= 8) +
    scale_color_manual(values = c("orange2", "red3"), name = "Harvest specification") +
    ylab("Catch (t)") +
    xlab(NULL) +
    guides(linetype = guide_legend(title = NULL)) +
    ggtitle("Gulf of Alaska") +
    theme(plot.title = element_text(hjust = 0.5))

  # Subarea with ABC
  subarea_hl_fig <- ggplot() +
    geom_bar(data= fish_catch_area, stat = "identity", alpha = 0.5, aes(x = Year, y = Catch), color = "grey") +
    geom_line(data=subarea_hl_dat, aes(x= Year, y= value, linetype = harv_specs), color = "orange2", linewidth= 1) +
    geom_line(data=subarea_hl_dat %>% filter(Year %in% c(ayr: (ayr+1))), aes(x = Year, y = value), color = "orange2", linetype = 2) +
    geom_point(data= subarea_hl_dat %>% filter(Year %in% (ayr+1)), aes(x= Year, y= value), color = "orange2", shape= 8) +
    #scale_color_manual(values = c("orange2", "red3"), name = "Harvest specification") +
    facet_wrap(~Area, nrow =2, scales = "free_y") +
    ylab("Catch (t)") +
    #guides(linetype = guide_legend(title = NULL)) +
    theme(legend.position = "none")

  # combine into single plot (GOA total on top)

  harv_specs_catch <- ggpubr::ggarrange(goa_hl_fig, subarea_hl_fig,
                                      nrow= 2, heights = c(1, 1.5))
  ggsave(here::here( year, model, 'apport', 'fig_catch_abc.png'), plot= harv_specs_catch,
         width = 8, height = 8)

}

## spatial distribution of survey catches ----
#sp_dist_bts_plot(years = c(2019, 2021, 2023))
sp_dist_bts_plot <- function(years , save = TRUE ) {
  libs_figs <- c("ggplot2", "mapdata", "mapproj")
  lapply(libs_figs, library, character.only = TRUE)

  bts_vast_dat <- readRDS(here::here(year, "data" , "user_input" , "Data_Geostat_Sebastes_variabilis.rds"))

  bts_vast_dat_sub <- bts_vast_dat %>%
    dplyr::filter(Year %in% c(years)) %>%
    dplyr::mutate(Catch = Catch_KG/1000)

  # stopped here for map
  ak <- map_data('worldHires','USA:Alaska')
  ak <- map_data('worldHires')

  ak <- subset(ak, long<0) #drop the end of the Aleutian Islands, or use world2Hires
  bts_catch_map <- ggplot() +
    geom_polygon(data=ak, aes(long, lat, group=group), fill=8, color="black") +
    #theme(panel.background = element_rect(fill = 'aliceblue')) +
    xlab(expression(paste(Longitude^o,~'W'))) +
    ylab(expression(paste(Latitude^o,~'N'))) +
    coord_map(xlim = c(-180, -130), ylim = c(50, 65)) +
    geom_point(data= bts_vast_dat_sub, aes(x= Lon, y= Lat, size = Catch) ,
               color = "blue", shape= 16, alpha = 0.2) +
    #geom_point(data= bts_vast_dat_sub, aes(x= Lon, y= Lat, color = Year),
    #           size= .5, shape= 16, alpha = 0.2) +
    facet_wrap(~ Year, nrow = 3) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

 # bts_catch_map
  if(isTRUE(save)) {
    ggsave(here::here( year, model, 'figs', 'f_btscatch.png'), plot= bts_catch_map,
           width = 5, height =8)
  }

} # end sp_dist_bts_plot

### Plot for pearson residuals

pearson_resid_plot <- function(comp_type = c("fac", "fsc", "sac", "ssc")[2] ,
                               label = 'Age', outlier=3, title_name ="Fishery") {

  lookup <- c(ind_lab = "Age", ind_lab = "Length")

  iss_tmp <- read.csv(here::here(year, folder, "tables", paste0("t_",comp_type,".csv") ) ) %>%
    dplyr::rename_if(stringr::str_detect(names(.), "X"), ~str_replace_all(., "X", "")) %>%
    dplyr::filter(name == "n_iss") %>%
    tidyr::pivot_longer(cols = -name, names_to = "year", values_to = "iss") %>%
    dplyr::mutate(year = as.integer(year)) %>%
    dplyr::select(year, iss)

  ind_lab <-  read.csv(here::here(year, folder, "tables", paste0("t_",comp_type,".csv") ) ) %>%
    dplyr::select(name) %>%
    dplyr::filter(name %nin% c("n_iss", "n_h", "n_ess")) %>% pull(name)

  dat_tmp <- read.csv(here::here(year, folder, "processed", paste0(comp_type, ".csv")) ) %>%
    dplyr::left_join(iss_tmp, by = join_by(year)) %>%
    dplyr::select(-X) %>%
    tidyr::pivot_wider( names_from = groups, values_from = value ) %>%
    dplyr::mutate(Residual = iss * (obs-pred) / sqrt(iss*pred)) %>%
    dplyr::mutate(id = ifelse(Residual<0, 'a', 'b'),
                  Outlier = ifelse(abs(Residual) >= outlier, "Yes", "No"),
                  Outlier = factor(Outlier, levels = c("No", "Yes")) ) %>%
    dplyr::rename(any_of(lookup))

  dat_tmp %>%
    ggplot(aes(year, ind_lab, color = Residual, size = Residual, shape = Outlier) ) +
    geom_point(show.legend=TRUE) +
    scale_size_area(guide=F, max_size = 5) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = dat_tmp, var = ind_lab, start=0) +
    afscassess::scale_x_tickr(data = dat_tmp, var = year) +
    scale_shape_manual(values = c(19,8), drop = FALSE) +
    ylab(label) +
    xlab('Year') +
    ggtitle(paste('Pearson -',title_name, label, sep = " ") ) +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )

} # end pearson resid plot
