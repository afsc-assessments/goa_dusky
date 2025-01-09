## Figure Functions for Presentation
## November 2024
## K. Omori
main_path <- 'C:/Users/krist/Desktop/Assessments/goa_dusky'
model = 'm22.5a'

# Catch figure

# fishery catch by area
fish_catch_tmp <- read.csv(here(main_path, year, "data", "output", "fish_catch.csv")) 

catch_plot <-  ggplot2::ggplot(fish_catch_tmp, ggplot2::aes(year, catch)) +
  ggplot2::geom_line(linewidth = 2, color = "blue3") +
  ggplot2::ylab("Catch (t)") +
  ggplot2::xlab("Year") +
  ggplot2::expand_limits(y = 0) +
  ggplot2::theme(axis.text=element_text(size=12),
                 axis.title=element_text(size=16,face="bold"),
                 legend.text=element_text(size=12))

ggplot2::ggsave(here(main_path, year, model, "figs", "catch_pres2.png"),
                plot= catch_plot,
                width = 6, height = 5, units = "in", dpi = 200)

# comp data observed
# fishery size comp

fsc_dat <- read.csv(here( main_path , year, model, "processed", "fsc.csv") ) %>%
  filter(groups == "obs")

fsc_plot <- ggplot(fsc_dat, aes(x= year, y= Length)) +
  geom_point(aes(size = value) ) +
  geom_hline(yintercept = 45, color = 'red4', linetype= 2) +
  labs(x= NULL) +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) 
# survey size
ssc_dat <- read.csv(here( main_path , year, model, "processed", "ssc.csv") ) %>%
  filter(groups == "obs")

ssc_plot <- ggplot(ssc_dat, aes(x= year, y= Length)) +
  geom_point(aes(size = value) ) +
  #geom_hline(yintercept = 45, color = 'red4', linetype= 2) +
  labs(x= NULL) +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) 

#fishery age
fac_dat <- read.csv(here( main_path , year, model, "processed", "fac.csv") ) %>%
  filter(groups == "obs")

fac_plot <- ggplot(fac_dat, aes(x= year, y= Age)) +
  geom_point(aes(size = value) ) +
  #geom_hline(yintercept = 45, color = 'red4', linetype= 2, linewidth = 1.5) +
  labs(x= NULL, title = "Fishery") +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size = 16, face= 'bold')) 

# survey age
sac_dat <- read.csv(here( main_path , year, model, "processed", "sac.csv") ) %>%
  filter(groups == "obs")

sac_plot <- ggplot(sac_dat, aes(x= year, y= Age)) +
  geom_point(aes(size = value) ) +
  #geom_hline(yintercept = 45, color = 'red4', linetype= 2, linewidth = 1.5) +
  labs(x= NULL, title = "Survey") +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size = 16, face= 'bold')) 

comp_plots <- cowplot::plot_grid(fsc_plot, ssc_plot, fac_plot, sac_plot, nrow = 1)
ggplot2::ggsave(plot = comp_plots, here(main_path, year, model, "figs", "comps_presentation.png"),
                width = 12, height = 6, units = "in", dpi = 200)
fsc_plot2 <- fsc_plot +labs(title = NULL)
fac_plot2 <- fac_plot +labs(title = NULL)
fish_comp_plots <- cowplot::plot_grid(fsc_plot2, fac_plot2, nrow = 1)
ggplot2::ggsave(plot = fish_comp_plots, here(main_path, year, model, "figs", "comps_fish_presentation.png"),
                width = 8, height = 6, units = "in", dpi = 200)

ssc_plot2 <- ssc_plot +labs(title = NULL)
sac_plot2 <- sac_plot +labs(title = NULL)
bts_comp_plots <- cowplot::plot_grid(ssc_plot2, sac_plot2, nrow = 1)
ggplot2::ggsave(plot = bts_comp_plots, here(main_path, year, model, "figs", "comps_bts_presentation.png"),
                width = 8, height = 6, units = "in", dpi = 200)
