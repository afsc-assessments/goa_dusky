# Apportionment Update 2024
# VAST vs. rema
# kristen.omori@noaa.gov
# vast code originally from Cecilia O'Leary
# updated date: July 16, 2024

## load ----
library(tidyverse)
library(rema)
library(here)
library(ggplot2)
library(mapdata)
library(mapproj)
library(afscassess)

## data and switches ----

yr <- 2024
species_code_race <- c(30152, 30150) # dusky (and dark_rockfish) & dusky and dark rockfish
species_name <- "Dusky"
sci_name <- "Sebastes_variabilis"

# data for map
data_geostat <- readRDS(here::here(yr, 'data', 'user_input', 'Data_Geostat_Sebastes_variabilis.rds'))

db_gap_index <- read.csv(here::here(yr, 'data','raw','goa_area_bts_biomass_data.csv')) %>%
  dplyr::mutate(type = "db_gap") %>%
  dplyr::filter(year >= 1990) %>%
  dplyr::select(year, area = regulatory_area_name, species_code, biomass = area_biomass, var= biomass_var) %>%
  dplyr::group_by(year, area ) %>%
  dplyr::summarise(biomass = sum(biomass, na.rm=T),
                   var = sum(var, na.rm=T))

# vast geostat data
# Cecilia O'Leary ran the VAST model with updated survey data last year
# lognomal error structure (the accepted model)

strata_limits_tmp <- data.frame('strata' = as.factor(c("Total","Western","Central","Eastern")),
                            'west_border' = c(-Inf,-Inf, -159,-147),
                            'east_border' = c(Inf,-159,-147,Inf ),
                            'Stratum' = c("Stratum_1", "Stratum_2", "Stratum_3", "Stratum_4"))

#vast_fit <- readRDS(here::here(yr, 'data', 'user_input', 'Sebastes_variabilisVASTfit.rds'))

vast_index <- read.csv(here::here(yr, 'data','user_input','vast_lognormal_index.csv')) %>%
  dplyr::left_join( strata_limits_tmp %>% select(area = strata, Stratum)) %>%
  dplyr::select(year = Time, area, biomass = Estimate, se = Std..Error.for.Estimate) %>%
  dplyr::mutate(type = "vast" ) %>% # units bio0mass = kg (need to convert to t later)
  dplyr::filter(year >= 1990,
                biomass > 0 )

## model runs with rema ----
# compare: db_gap_index
#          vast_index
#          db_gap_index + rema (base case)
#          vast_index + rema

# db_gap_index
db_gap_index_mod <- db_gap_index %>%
  dplyr::filter(!area == "Total") %>%
  dplyr::mutate( cv = sqrt(var)/biomass,
                 lci = (biomass- (sqrt(var)*1.96)),
                 uci =  (biomass+ (sqrt(var)*1.96)) ) %>%
  dplyr::select(year, strata= area, biomass, cv, lci, uci) %>%
  dplyr::mutate(strata = case_when(strata == 'WESTERN GOA' ~ "Western",
                                   strata == 'CENTRAL GOA' ~ 'Central',
                                   strata == 'EASTERN GOA' ~ 'Eastern' ),
                strata = factor(strata, levels = c("Western", "Central", "Eastern"))) %>%
  dplyr::ungroup()

input_db_rema <- prepare_rema_input( model_name = 'db_rema',
                                biomass_dat = db_gap_index_mod,
                                end_year = yr,
                                zeros = list(assumption = 'NA' ) )
mod_db_rema <- fit_rema(input_db_rema)
out_db_rema <- tidy_rema(mod_db_rema)

out_db_rema$proportion_biomass_by_strata

input_db_rema_1pe <- prepare_rema_input( model_name = 'db_rema',
                                     biomass_dat = db_gap_index_mod,
                                     end_year = yr,
                                     zeros = list(assumption = 'NA' ),
                                     PE_options = list(pointer_PE_biomass = c(1,1,1)))
mod_db_rema_1pe <- fit_rema(input_db_rema_1pe )
out_db_rema_1pe <- tidy_rema(mod_db_rema_1pe)

out_db_rema_1pe$proportion_biomass_by_strata

# vast_index
vast_index_mod <- vast_index %>%
  dplyr::filter(!area == "Total") %>%
  dplyr::mutate( cv = se/biomass,
                 lci = (biomass- (se*1.96))/1000, # calc lci and convert units
                 uci =  (biomass+ (se*1.96))/1000, # calc uci and convert units
                 biomass = biomass/1000 ) %>% # converting units to t (after cv calc)
  dplyr::select(year, strata= area, biomass, cv, lci, uci) %>%
  dplyr::mutate(strata = factor(strata, levels = c("Western", "Central", "Eastern"))) %>%
  dplyr::arrange(year) %>%
  dplyr::ungroup()

# Don't need to run this
input_vast_rema <- prepare_rema_input( model_name = 'vast_rema',
                                     biomass_dat = vast_index_mod,
                                     end_year = yr,
                                     zeros = list(assumption = 'NA' ) )
mod_vast_rema <- fit_rema(input_vast_rema)
out_vast_rema <- tidy_rema(mod_vast_rema)

# DON'T NEED!!!
out_vast_rema$proportion_biomass_by_strata

# organize vast data with proportions
proportion_biomass_vast <- vast_index %>%
  dplyr::select(year, area, model_name = type , biomass ) %>%
  dplyr::mutate(biomass = biomass/1000,) %>% # converting units from kg to t
  tidyr::pivot_wider( names_from = area , values_from = biomass) %>%
  dplyr::mutate( Western = Western/ Total,
                 Central = Central/ Total,
                 Eastern = Eastern/ Total)

## tables for Sept. PT ----
prop_mods_dat  <- out_db_rema$proportion_biomass_by_strata %>% # base
  dplyr::mutate(model = "db+rema (base)" ) %>%
  dplyr::select(model, year, Western, Central, Eastern) %>%
  dplyr::full_join(proportion_biomass_vast %>%
                     dplyr::select(model = model_name, year, Western, Central, Eastern)) %>%
  #dplyr::full_join(out_vast_rema$proportion_biomass_by_strata %>%
  #                   dplyr::mutate(model = "vast+rema") %>%
  #                   dplyr::select(model, year, Western, Central, Eastern)) %>%
  tidyr::pivot_longer(c(Western, Central, Eastern), names_to = "area", values_to = "proportion" ) %>%
  dplyr::mutate(area = factor(area, levels = c("Western", "Central", "Eastern")))

prop_mods_table <- prop_mods_dat %>%
  dplyr::filter(year %in% c(2021,2022, 2023)) %>%
  dplyr::mutate(percentage = round(proportion*100, 1)) %>%
  dplyr::select(Year=year, Area=area, model, percentage) %>%
  tidyr::pivot_wider(names_from = model, values_from = percentage)

write.csv( prop_mods_table ,
           here::here(yr, 'sept_pt', 'apportionment', 'table_proportions.csv' ), row.names = F)

## figures for Sept. PT ----

# biomass caught by area

biom_bts_area <- db_gap_index %>%
  dplyr::filter(!area == "Total") %>%
  dplyr::mutate(lci = (biomass- (sqrt(var)*1.96)),
                uci =  (biomass+ (sqrt(var)*1.96)),
                Area = factor(area, levels = c("Western", "Central", "Eastern")))

biom_bts_fig <- ggplot(biom_bts_area, aes( x = year, y = biomass, color = Area)) +
  geom_line(linewidth= 1.25 ) +
  geom_point(size= 2) +
  geom_errorbar(aes(ymin= lci, ymax= uci, color = Area, width = 0.5)) +
  scale_color_viridis_d(alpha = 0.75) +
  expand_limits(y=0) +
  xlab("Year") +
  ylab("Biomass (t)") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw()

ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_bts_biom_area.png'), plot= biom_bts_fig,
       width = 6, height = 4)

# biomass by region by type: db, vast, db+rema (base), vast+rema
# Do individual graphs by region
biom_mods_dat  <- vast_index_mod %>%                         # vast index
  dplyr::mutate( Model = "vast") %>%
  dplyr::select( year, area = strata, Model, biomass, lci, uci ) %>%
  dplyr::full_join( out_db_rema$biomass_by_strata %>%          # db + rema (base)
                      dplyr::mutate( Model = "db+rema (base)" ) %>%
                      dplyr::filter( !is.na(obs)) %>%
                      dplyr::select( year, area= strata, Model, biomass= pred, lci= pred_lci, uci= pred_uci ) ) %>%
  dplyr::full_join( biom_bts_area %>%          # db + rema (base)
                      dplyr::mutate( Model = "db (gap) " ) %>%
                      dplyr::select( year, area , Model, biomass) ) %>%
  #dplyr::full_join( out_vast_rema$biomass_by_strata %>%          # db + vast
  #                    dplyr::mutate( index = "vast+rema") %>%
  #                    dplyr::filter( !is.na(obs)) %>%
  #                    dplyr::select( year, area= strata, index, biomass= pred ) ) %>%
  dplyr::mutate(area = factor(area, levels = c("Western", "Central", "Eastern")))

ko_tmp <- biom_mods_dat %>%
  tidyr::pivot_wider(names_from = "Model", values_from = "biomass")

index_area_fig <- ggplot(biom_mods_dat, aes( x= year, y= biomass, color = Model)) +
  geom_line(linewidth= 1.15 ) +
  geom_point(size= 1.5) +
  geom_errorbar(aes(ymin= lci, ymax= uci, color = Model, width = 0.5)) +
  #scale_color_viridis_d(alpha = 0.75) +
  scale_color_manual(values = c("grey40","blue4", "darkgoldenrod1")) +
  expand_limits(y=0) +
  xlab("Year") +
  ylab("Estimated Biomass (t)") +
  scale_x_continuous(breaks = seq(1990,2020,5)) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~area, nrow= 3, scales= "free") +
  theme_bw()

ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_index_area.png'), plot= index_area_fig,
       width = 7, height = 6)

# proportions by regions by type: db, vast, db+rema (base), vast+rema
# Do single graph by region and models
# Individual graphs by of stacked by area

prop_mods_dat  <- out_db_rema$proportion_biomass_by_strata %>% # base
  dplyr::mutate(Model = "db+rema (base)" ) %>%
  dplyr::select(Model, year, Western, Central, Eastern) %>%
  dplyr::full_join(proportion_biomass_vast %>%
                     dplyr::select(Model = model_name, year, Western, Central, Eastern)) %>%
  #dplyr::full_join(out_vast_rema$proportion_biomass_by_strata %>%
  #                   dplyr::mutate(model = "vast+rema") %>%
  #                   dplyr::select(model, year, Western, Central, Eastern)) %>%
  tidyr::pivot_longer(c(Western, Central, Eastern), names_to = "area", values_to = "proportion" ) %>%
  dplyr::mutate(Area = factor(area, levels = c("Western", "Central", "Eastern")))

prop_mods_dat %>% group_by(Model, Area) %>%
  summarise(max_prop = max(proportion), min_prop = min(proportion), mean_prop = mean(proportion, na.rm=T))

proportion_mod_fig <- ggplot(prop_mods_dat , aes( x= year, y= proportion, fill = Area)) +
  geom_bar( position= "stack" , stat = 'identity' ) +
  scale_fill_viridis_d(alpha = 0.75) +
  #scico::scale_fill_scico_d(palette = "roma") +
  xlab("Year") +
  ylab("Proportion") +
  #scale_x_continuous(breaks = seq(1990,2020,5)) +
  facet_wrap(~Model, ncol= 3) +
  theme_bw()

ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_prop_area_mod.png'), plot= proportion_mod_fig,
       width = 6, height = 4)

proportion_mod_fig2 <- ggplot(prop_mods_dat , aes( x= year, y= proportion, color = Area, lty=Model)) +
  geom_line(linewidth= 1.15 ) +
  #geom_point(size= 1.5) +
  xlab("Year") +
  ylab("Proportion") +
  #scale_x_continuous(breaks = seq(1990,2020,5)) +
  scale_color_viridis_d(alpha = 0.75) +
  theme_bw()

ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_prop_area_mod_single.png'), plot= proportion_mod_fig2,
       width = 6, height = 4)

# Fisheries catch with ABC last 10 years (2014-2023)
# Read in fisheries catch by area; read in ABC file (use this to calc prop)
dusk_fsh_catch <- read.csv(here::here(yr, "data", "raw", "fsh_catch_data.csv")) %>%
  dplyr::filter(agency_species_code %in% c(154, 172),
                year %in% c(2014:2023)) %>%
  dplyr::mutate(area = case_when(fmp_subarea == "WG" ~ "Western",
                                 fmp_subarea == "CG" ~ "Central",
                                 fmp_subarea == "WY" ~ "Eastern",
                                 fmp_subarea == "SE" ~ "Eastern")) %>%
  dplyr::group_by(year, area) %>%
  dplyr::summarise(weight = sum(weight_posted, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::select(year, area, weight) %>%
  dplyr::mutate(Area = factor(area, levels = c("Western", "Central", "Eastern")))

dusk_abc <- read.csv(here::here(yr, "data", "user_input", "abc_tac_ko.csv")) %>%
  dplyr::filter(Year %in% c(2014:2023)) %>%
  dplyr::select(year = Year, ABC)

dusk_subABC <- prop_mods_dat %>%
  dplyr::filter(year %in% c(2014:2023)) %>%
  dplyr::full_join( dusk_abc , by = "year") %>%
  dplyr::mutate(subABC = proportion*ABC) %>%
  #dplyr::full_join(dusk_fsh_catch) %>%
  dplyr::mutate(Area = factor(area, levels = c("Western", "Central", "Eastern")))

fsh_catch_abc_fig <- ggplot() +
  geom_bar(data= dusk_fsh_catch, stat = "identity", alpha = 0.5, aes(x = year, y = weight), color = "grey") +
  geom_line(data=dusk_subABC, aes(x= year, y= subABC, color= Model), linewidth= 1) +
  scale_color_manual(values = c("blue4", "darkgoldenrod1")) +
  ylab("Catch (t)") +
  facet_wrap(~Area) +
  theme_bw()

ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_catch_abc.png'), plot= fsh_catch_abc_fig,
       width = 7, height = 4)

tmp_subABC <- prop_mods_dat %>%
  dplyr::filter(year %in% c(2014:2023)) %>%
  dplyr::full_join( dusk_abc , by = "year") %>%
  dplyr::mutate(subABC = proportion*ABC) %>%
  dplyr::full_join(dusk_fsh_catch) %>%
  dplyr::mutate(Area = factor(area, levels = c("Western", "Central", "Eastern")),
                exceeds = -(subABC-weight)) %>%
  dplyr::filter(exceeds > 0)

# map of catch biomass

map_geostat_dat <- data_geostat %>%
  dplyr::filter(Catch_KG > 0 ) %>%
  dplyr::select(Catch_KG, Year, Lat, Lon)

ak <- map_data('worldHires')
ak <- subset(ak, long<0)

surv_map_pres2 <- ggplot(map_geostat_dat) +
  geom_point(aes(x= Lon, y= Lat, size= Catch_KG), shape= 16, alpha = 0.6, color= "blue") +
  geom_polygon(data=ak, aes(long, lat, group=group), fill=8, color="grey30") +
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'N'))) +
  coord_map(xlim = c(-180, -130), ylim = c(50.5, 63)) +
  geom_vline(xintercept = c(-159, -147), color = "grey30", linetype = "dashed") +
  facet_wrap(~Year) +
  afscassess::theme_report()
  #theme_bw() +
  #theme(panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(),
  #      legend.position = "bottom")

ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_survey_map_pts.png'), plot= surv_map_pres2,
       width = 12, height = 7)

# map for general areas

ak <- map_data('worldHires')
ak <- subset(ak, long<0)

goa_map <- ggplot() +
  geom_polygon(data=ak, aes(long, lat, group=group), fill=8, color="grey30") +
  annotate("text", x = c(52, 52, 52) , y= c(-170, -155, -140), label = c("WG", "CG", "EG")) +
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'N'))) +
  coord_map(xlim = c(-180, -130), ylim = c(50.5, 63)) +
  #annotate("text", x = c(52, 52, 52) , y= c(-170, -155, -140), label = c("WG", "CG", "EG")) +
  geom_vline(xintercept = c(-159, -147), color = "grey30", linetype = "dashed")
  #annotate("text", x = c(52, 52, 52) , y= c(-170, -155, -140), label = c("WG", "CG", "EG")) +
  #afscassess::theme_report()
ggsave(here::here( yr, 'sept_pt', 'apportionment', 'figures' , 'fig_goa_map.png'), plot= goa_map,
       width = 8, height = 4)
