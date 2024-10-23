# ADMB Model comparisons 2024
# Using 2022 model and data (For Sept PT)
# kristen.omori@noaa.gov
# updated date: August 12, 2024

## libraries ----
library(tidyverse)
library(afscassess)
library(R2admb)
library(here)
library(scales)
library(ggplot2)

## globals ----
year_mod = 2022
year = 2024
species = "DUSK"
area = "goa"
rec_age = 4
plus_age = 30
len_bins = read.csv(here::here(year, 'data', 'user_input', 'lbins.csv'))

base_mod = "m1.base"
srv_mod = "m2.srv"
srvproj_mod = "m3.srvproj"

# setting folder location
folder_temp <- paste0(getwd(), "/", year,"/")

# In each model folder have:
# .tpl, .cpp, .ctl, .dat
# Data files= goa_dr_yyyy.dat AND MAT.DAT (static)

## To create new .exe file for each model ----

# Compiling model to get .exe file isn't working in R
#R2admb::compile_admb('base', admb_errors = "warn") # just creates .cout

# To compile in ADMB command line
#1. cd C:\Users\Kristen.Omori\Work\Assessments\dusk_test_2\dusk_test_2\2022\m22.3a
#       this will change directory/ path
#2. type: admb base.tpl
#       this will make the .exe file
#3. type: base.exe
#       this will run the admb model (make sure it runs!)

## Run all models ----

## Base (previously accepted model)
setwd(paste0(folder_temp, base_mod)) ; getwd()

R2admb::run_admb('base', verbose= TRUE)

rep_base <- readLines(paste0('base.rep'))
std_base <- readLines(paste0('base.std'))[-1]

report_admb(report_name= rep_base, std_name = std_base, model_year = year_mod, year = year, model_name = base_mod)

std_base <- data.frame(text = std_base) %>%
  dplyr::mutate(text = stringr::str_squish(text)) %>%
  tidyr::separate(text, c("index", "name", "value", "std_dev"), sep = " ")

write.csv(std_base, here::here(year, base_mod, "processed", "std_output.csv"), row.names = FALSE )

## Model 2: updated lognormal survey error structure
# Make sure model has a new .exe
setwd(paste0(folder_temp, srv_mod)); getwd()

R2admb::run_admb('base', verbose= TRUE)

rep_srvmod <- readLines(paste0('base.rep'))
std_srvmod <- readLines(paste0('base.std'))[-1]

report_admb(report_name= rep_srvmod, std_name = std_base, model_year = year_mod, year = year, model_name = srv_mod)

std_srvmod <- data.frame(text = std_srvmod) %>%
  dplyr::mutate(text = stringr::str_squish(text)) %>%
  tidyr::separate(text, c("index", "name", "value", "std_dev"), sep = " ")

write.csv(std_srvmod, here::here(year, srv_mod, "processed", "std_output.csv"), row.names = FALSE )

## Model 3: lognormal survey error + corrected projection of the recruitment age starting year
# Make sure model has a new .exe
setwd(paste0(folder_temp, srvproj_mod)); getwd()

R2admb::run_admb('base', verbose= TRUE)

rep_srvprojmod <- readLines(paste0('base.rep'))
std_srvprojmod <- readLines(paste0('base.std'))[-1]

report_admb(report_name= rep_srvprojmod, std_name = std_base, model_year = year_mod, year = year, model_name = srvproj_mod)

std_srvprojmod <- data.frame(text = std_srvprojmod) %>%
  dplyr::mutate(text = stringr::str_squish(text)) %>%
  tidyr::separate(text, c("index", "name", "value", "std_dev"), sep = " ")

write.csv(std_srvprojmod, here::here(year, srvproj_mod, "processed", "std_output.csv"), row.names = FALSE )

## Read in input/ output for model comparison ----

# Use if not running above part of script
std_base <- read.csv(here::here(year, "alt_models_september", "m1.base", "processed", "std_output.csv"))
std_srvmod <- read.csv(here::here(year, "alt_models_september", "m2.srv_v2", "processed", "std_output.csv"))
std_srvprojmod <- read.csv(here::here(year, "alt_models_september", "m3.srvproj_v3", "processed", "std_output.csv"))

# Likelihood values
LL_vals <- read.csv(here::here(year, "alt_models_september", base_mod, "processed", "Table_likelihood_values.csv")) %>%
  dplyr::full_join(read.csv(here::here(year, "alt_models_september", srv_mod, "processed", "Table_likelihood_values.csv"))) %>%
  dplyr::full_join(read.csv(here::here(year, "alt_models_september", srvproj_mod, "processed", "Table_likelihood_values.csv")))

# Parameter estimate values
param_est <- read.csv(here::here(year, base_mod, "processed", "Table_parameter_est.csv")) %>%
  dplyr::full_join(read.csv(here::here(year, srv_mod, "processed", "Table_parameter_est.csv"))) %>%
  dplyr::full_join(read.csv(here::here(year, srvproj_mod, "processed", "Table_parameter_est.csv")))

# survey biomass predicted and observed with std

vast_index <- read.csv(here::here(year, "data", "user_input", "vast_lognormal_2022.csv")) %>%
  dplyr::mutate(model = "vast_index (observed)")

survB_dat <-  read.csv(here::here(year, "alt_models_september", base_mod, "processed", "df_surveyB.csv")) %>%
  dplyr::mutate(model = base_mod) %>%
  dplyr::full_join(read.csv(here::here(year, "alt_models_september", srv_mod, "processed", "df_surveyB.csv")) %>%
                     dplyr::mutate(model = srv_mod)) %>%
  dplyr::mutate(survB_obs = NA) %>%
  dplyr::bind_rows(vast_index %>%
                     dplyr::mutate(survB_obs = biomass) %>%
                     dplyr::select(year, model, survB_obs, lci, uci) ) %>%
  dplyr::mutate(model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj",
                                  model == "vast_index (observed)" ~ "vast_index (observed)"))


# spawning biomass estimates with std
mod_yrs <- 1977:2022
ssb_dat <- std_base %>%
  dplyr::filter(name == "spawn_biom") %>%
  dplyr::mutate(year = mod_yrs,
                model = base_mod) %>%
  dplyr::full_join( std_srvmod %>%
                      dplyr::filter(name == "spawn_biom") %>%
                      dplyr::mutate(year = mod_yrs,
                                    model = srv_mod) ) %>%
  dplyr::full_join( std_srvprojmod %>%
                      dplyr::filter(name == "spawn_biom") %>%
                      dplyr::mutate(year = mod_yrs,
                                    model = srvproj_mod) ) %>%
  dplyr::select(-index) %>%
  dplyr::mutate(model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj"))


# Avg predicted recruitment with std

pred_rec_dat <- std_base %>%
  dplyr::filter(name == "pred_rec") %>%
  dplyr::mutate(year = mod_yrs,
                model = base_mod) %>%
  dplyr::full_join( std_srvmod %>%
                      dplyr::filter(name == "pred_rec") %>%
                      dplyr::mutate(year = mod_yrs,
                                    model = srv_mod) ) %>%
  dplyr::full_join( std_srvprojmod %>%
                      dplyr::filter(name == "pred_rec") %>%
                      dplyr::mutate(year = mod_yrs,
                                    model = srvproj_mod) ) %>%
  dplyr::mutate(model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj"))

avg_rec_dat <- std_base %>%
  dplyr::filter(name == "avg_rec") %>%
  dplyr::mutate(model = base_mod) %>%
  dplyr::full_join( std_srvmod %>%
                      dplyr::filter(name == "avg_rec") %>%
                      dplyr::mutate(model = srv_mod) ) %>%
  dplyr::full_join( std_srvprojmod %>%
                      dplyr::filter(name == "avg_rec") %>%
                      dplyr::mutate(model = srvproj_mod) ) %>%
  dplyr::mutate(model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj"))

# Fully selected F

F_table <- read.csv(here::here(year, "alt_models_september", base_mod, "processed", "df_fully_selected_F.csv")) %>%
  dplyr::mutate(model = base_mod) %>%
  dplyr::full_join(read.csv(here::here(year,"alt_models_september", srv_mod, "processed", "df_fully_selected_F.csv")) %>%
                     dplyr::mutate(model = srv_mod) ) %>%
  dplyr::full_join(read.csv(here::here(year, "alt_models_september", srvproj_mod, "processed", "df_fully_selected_F.csv")) %>%
                     dplyr::mutate(model = srvproj_mod) ) %>%
  dplyr::mutate(model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj"))

## Model comparison ----
# 1. Likelihood values for: (LL_vals)
#       catch, survey biomass, fishery ages, survey ages,
#       fishery lengths, recruitment devs, F regularity, q prior
# 2. Parameter estimates for: (param_est_tab)
#       q, avg rec, F40, total B, spawning B, B100, B40, ABC
# 3. Fig. survey fit with normal and lognormal comparison (survB_dat)
# 4. Fig. SSB fits with all models (ssb_dat)
# 5. Fig. recruitment (pred_rec_dat)
# 6. Fig. fully selected F (F_table)

# Likelihood values
LL_vals_table <- LL_vals %>%
  dplyr::mutate(LL = round(LL, 3),
                model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj") ) %>%
  tidyr::pivot_wider(names_from = model, values_from = LL)

write.csv(LL_vals_table, here::here(year, "sept_pt", "model_comparison", "Likelihood_table.csv"), row.names = F)

# Parameter estimates
param_est_table <- param_est %>%
  dplyr::mutate(Estimates = ifelse(Parameter %in% c("sigmaR", "q", "avg rec", "F40"),as.character(round(Estimates, 3)),
                                  ifelse(Parameter %in% c("Total Biomass", "SSB", "B100", "B40", "ABC"),
                                  scales::comma(round(Estimates, 0)), NA)),
                model = case_when(model == "m1.base" ~ "m22.3a_base",
                                  model == "m2.srv" ~ "m22.4a_srv",
                                  model == "m3.srvproj" ~ "m22.5a_srvproj")) %>%
  tidyr::pivot_wider(names_from = model, values_from = Estimates)

write.csv(param_est_table, here::here(year, "sept_pt", "model_comparison", "Parameter_est_table.csv"), row.names = F)

# Figure: survey fit
head(survB_dat)

surv_B_fig <- ggplot(survB_dat, aes(x= year, y= survB_pred/1000, color= model)) +
  geom_line(linewidth= 1.15) +
  scale_color_manual(values = c("blue4", "darkgoldenrod1", "grey30")) +
  geom_point(aes(x=year, y= survB_obs/1000), color = "grey30", size =1.2) +
  geom_errorbar(aes(ymin= lci/1000, ymax= uci/1000, width = 0.5), color = "grey30") +
  ylab("Survey biomass (kt)") +
  theme_bw()

ggsave(here::here( year, 'sept_pt', 'model_comparison', 'fig_survey_fit.png'), plot= surv_B_fig,
       width = 7, height = 4)

survB_dat %>% filter(!is.na(survB_pred)) %>%
  select(year, survB_pred, model) %>%
  pivot_wider(names_from = model, values_from = survB_pred) %>%
  mutate(perc_diff = ((m22.4a_srv-m22.3a_base)/m22.3a_base)*100 )

# Figure: ssb
head(ssb_dat)

ssb_fig <- ggplot(ssb_dat, aes(x= year, y= value/1000, color= model, fill= model)) +
  geom_line(linewidth= 1.15) +
  geom_ribbon(aes(ymin= (value-std_dev)/1000, ymax= (value+std_dev)/1000 ) , alpha= 0.1, linetype=0) +
  scale_color_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
  scale_fill_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
  ylab("Spawning biomass (kt)") +
  xlab(NULL) +
  theme_bw()

ggsave(here::here( year, 'sept_pt', 'model_comparison', 'fig_ssb.png'), plot= ssb_fig,
       width = 7, height = 4)

tmp_ssb_diff <- ssb_dat %>%
  select(year, value, model) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(perc_diff_srv = ((m22.4a_srv-m22.3a_base)/m22.3a_base)*100,
         perc_diff_base_srvproj = ((m22.5a_srvproj-m22.3a_base)/m22.3a_base)*100,
         perc_diff_srv_proj = ((m22.5a_srvproj-m22.4a_srv)/m22.4a_srv)*100)

# Figure: recruitment
head(pred_rec_dat)

recruit_fig <- ggplot(pred_rec_dat, aes(x= year, y= value, color= model, fill= model)) +
  geom_line(linewidth= 1.15) +
  geom_ribbon(aes(ymin= (value-std_dev), ymax= (value+std_dev) ) , alpha= 0.1, color= NA) +
  scale_color_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
  scale_fill_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
  ylab("Recruitment (millions)") +
  #xlab(NULL) +
  theme_bw()

ggsave(here::here( year, 'sept_pt', 'model_comparison', 'fig_recruit.png'), plot= recruit_fig,
       width = 7, height = 4)

# Figure: fully selected F
head(F_table)

sel_F_fig <- ggplot(F_table, aes(x= year, y= F_vec, color= model)) +
  geom_line(linewidth= 1.15) +
  scale_color_manual(values = c("blue4", "darkgoldenrod1", "green4")) +
  ylab("Fishing mortality rate (F)") +
  xlab(NULL) +
  theme_bw()

ggsave(here::here( year, 'sept_pt', 'model_comparison', 'fig_F_rate.png'), plot= sel_F_fig,
       width = 7, height = 4)

# Combine graphics into one panel
sb_F_rec_fig <- ggpubr::ggarrange(ssb_fig, sel_F_fig, recruit_fig,
                                  labels= c("A", "B", "C"),
                                  common.legend = TRUE,
                                  legend= "right",
                                  nrow= 3)
ggsave(here::here( year, 'sept_pt', 'model_comparison', 'fig_sb_F_rec.png'), plot= sb_F_rec_fig,
       width = 7, height = 9)
