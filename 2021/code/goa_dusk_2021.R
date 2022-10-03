# load ----
library(groundfishr)
library(gfdata)
library(tidyverse)
library(vroom)
library(funcr)
theme_set(theme_report())
library(here)
library(scico)

# globals ----
year <- 2021
species <- "DUSK"
TAC <- c(3957, 3700, 3676)
afsc_user <- "WILLIAMSB"
afsc_pwd <- "H@ines2010_2021"
akfin_user <- "bwilliams"
akfin_pwd <- "jnu$6350"
model <- "base"
dat <- "goa_dr_2020"
dat2 <- "goa_dr_2021"

# create a directory for the projection
dir.create(here::here(year, "proj", "data"), recursive = TRUE)

accepted_model(2020, "m15.5a", 2021)

# update catch data
goa_dusk(year, akfin_user, akfin_pwd, afsc_user, afsc_pwd, off_yr = TRUE)
clean_catch(year = year, species = species, TAC = TAC)
ts_biomass(year, id = "db")

catch <- vroom::vroom(here::here(year, "data", "output", "fsh_catch.csv")) %>%
                      mutate(Catch = round(Catch))

# not enough time to adjust new functions - so switch to og code

# read in observer/catch landings data
vroom::vroom(here::here(year, "data", "raw", 
              paste0(fishery, "_catch_data.csv"))) -> catch_data
vroom::vroom(here::here(year, "data", "raw",  paste0(fishery, "_obs_data.csv"))) -> obs_data
endyr <- as.numeric(substr(Sys.time(),1,4))

# Estimate ratio of catch from current date to end of year (Endyr_ratio)
yr <- seq(endyr-3,endyr-1)
Endyr_C <- matrix(nrow=length(yr),ncol=2)
colnames(Endyr_C)<-c("Oct_C","Total_C")
rownames(Endyr_C)<-yr
for(y in 1:length(yr)){
  Data<-subset(obs_data,obs_data$YEAR==yr[y])
  Data_pre<-subset(Data,as.Date(Data$HAUL_DATE)<=paste(yr[y],substr(max(as.Date(catch_data$WEEK_END_DATE)),5,10),sep=""))
  Endyr_C[y,1]<-sum(Data_pre$EXTRAPOLATED_WEIGHT)
  Endyr_C[y,2]<-sum(Data$EXTRAPOLATED_WEIGHT)}
Endyr_ratio<-1+(sum(Endyr_C[,2])-sum(Endyr_C[,1]))/sum(Endyr_C[,1])
Catch_end_date<-max(as.Date(catch_data$WEEK_END_DATE))

# Compute total catch and estimate yeild ratio of TAC to catch
yr<-seq(endyr-3,endyr)
C<-matrix(nrow=length(yr),ncol=1)
rownames(C)<-yr
for(y in 1:length(yr)){
  Data<-subset(catch_data,catch_data$YEAR==yr[y])
  C[y,1]<-sum(Data$WEIGHT_POSTED)}
C[length(C[,1]),1]<-C[length(C[,1]),1]*Endyr_ratio
yld_rat<-mean(C[(length(C)-3):(length(C)-1)]/TAC)


#====================================================================================================
# Step 2: Get initial projection model data files set up
#====================================================================================================

# Read in proj.dat and wite to projection model data file
file.copy(here::here(year, "base", "proj.dat"),
          here::here(year, "proj", "data", "goa_dusky.dat"))

# Get setup.dat file setup
setup<-readLines(here::here(year, "proj","setup.dat"),warn=FALSE)
L_endyr<-grep("# beg_yr_label",setup)
# test if endyr is for a full or partial assessment
# if(endyr %% 2 == 0){ # For odd-year full assessment
if(endyr %% 2 == 1){ # For even-year full assessment
  setup[L_endyr]<-paste(endyr-1," #_Begin Year")}else{
    setup[L_endyr]<-paste(endyr," #_Begin Year")}
write.table(setup,file=here::here(year, "proj","setup.dat"),quote=F,row.names=F,col.names=F)


#=================================================================================
# Step 3: Run Max F projection scenario
#====================================================================================================

# Setup spp_catch.dat file (toggled to test for on/off year)
L_1<-"#_Number_of_years with specified catch"
# if(endyr %% 2 == 0){ # For odd-year full assessment
if(endyr %% 2 == 1){ # For even-year full assessment
  L_2<-2}else{
    L_2<-1}
L_3<-"# Number of species"
L_4<-1
L_5<-"# data files for each species"
L_6<-"data/goa_dusky.dat"
L_7<-"# ABC Multipliers"
L_8<-1
L_9<-"# Population scalars"
L_10<-1000
L_11<-"# Number of TAC model categories"
L_12<-1
L_13<-"# TAC model indices (for aggregating)"
L_14<-1
L_15<-"# Catch in each future year" # Includes toggle for on/off year
# if(endyr %% 2 == 0){ # For odd-year full assessment
if(endyr %% 2 == 1){ # For even-year full assessment
  L_16<-paste(paste(endyr-1,round(C[which(as.numeric(rownames(C))==(endyr-1)),1],digits=4),sep="\t"),"# Finalized previous year catch",sep=" ")
  L_17<-paste(paste(endyr,round(C[which(as.numeric(rownames(C))==(endyr)),1],digits=4),sep="\t"),"# Estimated from catch thru",Catch_end_date,"with expansion factor =",Endyr_ratio,sep=" ")
  spp_catch<-c(L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10,L_11,L_12,L_13,L_14,L_15,L_16,L_17)}else{
    L_16<-paste(paste(endyr,round(C[which(as.numeric(rownames(C))==(endyr)),1],digits=4),sep="\t"),"# Estimated from catch thru",Catch_end_date,"with expansion factor =",Endyr_ratio,sep=" ")
    spp_catch<-c(L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10,L_11,L_12,L_13,L_14,L_15,L_16)}
write.table(spp_catch,file=here::here(year, "proj", "spp_catch.dat"),quote=F,row.names=F,col.names=F)
write.table(spp_catch,file=here::here(year, "proj", "data","goa_dusky_max_spcat.dat"),quote=F,row.names=F,col.names=F)

# Run model
#	NOTE: there's a warning message that comes up with this, you can ignore it,
#		it has to do with the ADMB version the projection model was originally
#		compiled with but has not effect on results
setwd(here::here(year, "proj"))
shell("main.exe")

# Read results from Max scenario
bigfile_write<-readLines(here::here(year, "proj", "bigfile.out"),warn=FALSE)
bigfile<-read.delim(here::here(year, "proj", "bigfile.out"),sep="", header=T)
F_profile<-readLines(here::here(year, "proj", "F_profile.out"),warn=FALSE)
means<-readLines(here::here(year, "proj", "means.out"),warn=FALSE)
percentdb<-readLines(here::here(year, "proj", "percentdb.out"),warn=FALSE)
percentiles<-readLines(here::here(year, "proj", "percentiles.out"),warn=FALSE)

# Make bigsum file
bigsum<-matrix(nrow=98,ncol=9)
colnames(bigsum)<-c("Alt","Stock","Year","ABC","OFL","Catch","SSB","F","Total_Biom")
Alt1<-matrix(nrow=14,ncol=9);Alt1[,1]<-1;Alt1[,2]<-as.character(bigfile$Spp[1]);yrs<-sort(unique(bigfile$Yr));Alt1[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==1);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt1[y,4:9]<-colMeans(bigsum_data)}
Alt2<-matrix(nrow=14,ncol=9);Alt2[,1]<-2;Alt2[,2]<-as.character(bigfile$Spp[1]);Alt2[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==2);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt2[y,4:9]<-colMeans(bigsum_data)}
Alt3<-matrix(nrow=14,ncol=9);Alt3[,1]<-3;Alt3[,2]<-as.character(bigfile$Spp[1]);Alt3[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==3);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt3[y,4:9]<-colMeans(bigsum_data)}
Alt4<-matrix(nrow=14,ncol=9);Alt4[,1]<-4;Alt4[,2]<-as.character(bigfile$Spp[1]);Alt4[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==4);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt4[y,4:9]<-colMeans(bigsum_data)}
Alt5<-matrix(nrow=14,ncol=9);Alt5[,1]<-5;Alt5[,2]<-as.character(bigfile$Spp[1]);Alt5[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==5);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt5[y,4:9]<-colMeans(bigsum_data)}
Alt6<-matrix(nrow=14,ncol=9);Alt6[,1]<-6;Alt6[,2]<-as.character(bigfile$Spp[1]);Alt6[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==6);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt6[y,4:9]<-colMeans(bigsum_data)}
Alt7<-matrix(nrow=14,ncol=9);Alt7[,1]<-7;Alt7[,2]<-as.character(bigfile$Spp[1]);Alt7[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==7);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt7[y,4:9]<-colMeans(bigsum_data)}
bigsum<-rbind(Alt1,Alt2,Alt3,Alt4,Alt5,Alt6,Alt7)

# Write results from max scenario
write.table(bigfile_write,file= here::here(year, "proj", "maxout", "bigfile.out"),quote=F,row.names=F,col.names=F)
write.table(F_profile,file= here::here(year, "proj", "maxout"," F_profile.out"),quote=F,row.names=F,col.names=F)
write.table(means,file= here::here(year, "proj", "maxout", "means.out"),quote=F,row.names=F,col.names=F)
write.table(percentdb,file= here::here(year, "proj", "maxout", "percentdb.out"),quote=F,row.names=F,col.names=F)
write.table(percentiles,file= here::here(year, "proj", "maxout", "percentiles.out"),quote=F,row.names=F,col.names=F)
write.table(bigsum,file= here::here(year, "proj", "maxout", "bigsum.dat"),quote=F,row.names=F,col.names=c("Alt","Stock","Year","ABC","OFL","Catch","SSB","F","Total_Biom"))


#=================================================================================
# Step 4: Run Author's F scenario
#====================================================================================================

# Setup spp_catch.dat file
L_1<-"#_Number_of_years with specified catch"
# if(endyr %% 2 == 0){ # For odd-year full assessment
if(endyr %% 2 == 1){ # For even-year full assessment
  L_2<-4}else{
    L_2<-3}
L_3<-"# Number of species"
L_4<-1
L_5<-"# data files for each species"
L_6<-paste("data/goa_dusky.dat",sep="")
L_7<-"# ABC Multipliers"
L_8<-1
L_9<-"# Population scalars"
L_10<-1000
L_11<-"# Number of TAC model categories"
L_12<-1
L_13<-"# TAC model indices (for aggregating)"
L_14<-1
L_15<-"# Catch in each future year" # Includes toggle for on/off year
# if(endyr %% 2 == 0){ # For odd-year full assessment
if(endyr %% 2 == 1){ # For even-year full assessment
  L_16<-paste(paste(endyr-1,round(C[which(as.numeric(rownames(C))==(endyr-1)),1],digits=4),sep="\t"),"# Finalized previous year catch",sep=" ")
  L_17<-paste(paste(endyr,round(C[which(as.numeric(rownames(C))==(endyr)),1],digits=4),sep="\t"),"# Estimated from catch thru",Catch_end_date,"with expansion factor =",Endyr_ratio,sep=" ")
  p1<-percentiles[grep("Catch",percentiles)[1]:grep("Spawning_Biomass",percentiles)[1]]
  L_18<-paste(paste(endyr+1,round(as.numeric(strsplit(p1[5], split=" ")[[1]][8])*1000*yld_rat,digits=4),sep="\t"),"# Estimated as Max F scenario catch*yieldratio =",yld_rat,sep=" ")
  L_19<-paste(paste(endyr+2,round(as.numeric(strsplit(p1[6], split=" ")[[1]][8])*1000*yld_rat,digits=4),sep="\t"),"# Estimated as Max F scenario catch*yieldratio",yld_rat,sep=" ")
  spp_catch<-c(L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10,L_11,L_12,L_13,L_14,L_15,L_16,L_17,L_18,L_19)
  }else{
    L_16<-paste(paste(endyr,round(C[which(as.numeric(rownames(C))==(endyr)),1],digits=4),sep="\t"),"# Estimated from catch thru",Catch_end_date,"with expansion factor =",Endyr_ratio,sep=" ")
    p1<-percentiles[grep("Catch",percentiles)[1]:grep("Spawning_Biomass",percentiles)[1]]
    L_17<-paste(paste(endyr+1,round(as.numeric(strsplit(p1[4],split=" ")[[1]][8])*1000*yld_rat,digits=4),sep="\t"),"# Estimated as Max F scenario catch*yieldratio",yld_rat,sep=" ")
    L_18<-paste(paste(endyr+2,round(as.numeric(strsplit(p1[5],split=" ")[[1]][8])*1000*yld_rat,digits=4),sep="\t"),"# Estimated as Max F scenario catch*yieldratio",yld_rat,sep=" ")
    spp_catch<-c(L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10,L_11,L_12,L_13,L_14,L_15,L_16,L_17,L_18)
    }
write.table(spp_catch,file=here::here(year, "proj", "spp_catch.dat"),quote=F,row.names=F,col.names=F)
write.table(spp_catch,file=here::here(year, "proj", "data", "goa_dusky_spcat.dat"),quote=F,row.names=F,col.names=F)

# Run model
#	NOTE: there's a warning message that comes up with this, you can ignore it,
#		it has to do with the ADMB version the projection model was originally
#		compiled with but has not effect on results
# setwd(paste(path,"/model",sep=""))
shell("main.exe")
setwd(here::here())

# Read results from Author's F scenario
bigfile_write<-readLines(here::here(year, "proj", "bigfile.out"),warn=FALSE)
bigfile<-read.delim(here::here(year, "proj", "bigfile.out"),sep="", header=T)
F_profile<-readLines(here::here(year, "proj", "F_profile.out"),warn=FALSE)
means<-readLines(here::here(year, "proj", "means.out"),warn=FALSE)
percentdb<-readLines(here::here(year, "proj", "percentdb.out"),warn=FALSE)
percentiles<-readLines(here::here(year, "proj", "percentiles.out"),warn=FALSE)

# Make bigsum file

bigfile %>% 
  as.data.frame() %>% 
  group_by(Alternative, Yr) %>% 
  summarise_all(list(mean)) %>% 
  filter(Yr %in% c(2022, 2023))  %>% 
  mutate(ABC = round(ABC * 1000),
         OFL = round(OFL * 1000),
         Catch = round(Catch * 1000),
         SSB = round(SSB * 1000),
         Tot_biom = round(Tot_biom * 1000))

rename(year = yr) %>% 
  vroom::vroom_write(here::here(year, "proj", "author_f", "bigsum.csv"))
bigsum<-matrix(nrow=98,ncol=9)
colnames(bigsum)<-c("Alt","Stock","Year","ABC","OFL","Catch","SSB","F","Total_Biom")
Alt1<-matrix(nrow=14,ncol=9);Alt1[,1]<-1;Alt1[,2]<-as.character(bigfile$Spp[1]);yrs<-sort(unique(bigfile$Yr));Alt1[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==1);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt1[y,4:9]<-colMeans(bigsum_data)}
Alt2<-matrix(nrow=14,ncol=9);Alt2[,1]<-2;Alt2[,2]<-as.character(bigfile$Spp[1]);Alt2[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==2);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt2[y,4:9]<-colMeans(bigsum_data)}
Alt3<-matrix(nrow=14,ncol=9);Alt3[,1]<-3;Alt3[,2]<-as.character(bigfile$Spp[1]);Alt3[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==3);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt3[y,4:9]<-colMeans(bigsum_data)}
Alt4<-matrix(nrow=14,ncol=9);Alt4[,1]<-4;Alt4[,2]<-as.character(bigfile$Spp[1]);Alt4[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==4);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt4[y,4:9]<-colMeans(bigsum_data)}
Alt5<-matrix(nrow=14,ncol=9);Alt5[,1]<-5;Alt5[,2]<-as.character(bigfile$Spp[1]);Alt5[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==5);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt5[y,4:9]<-colMeans(bigsum_data)}
Alt6<-matrix(nrow=14,ncol=9);Alt6[,1]<-6;Alt6[,2]<-as.character(bigfile$Spp[1]);Alt6[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==6);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt6[y,4:9]<-colMeans(bigsum_data)}
Alt7<-matrix(nrow=14,ncol=9);Alt7[,1]<-7;Alt7[,2]<-as.character(bigfile$Spp[1]);Alt7[,3]<-yrs
for(y in 1:length(yrs)){bigsum_data<-subset(bigfile,bigfile$Alternative==7);bigsum_data<-subset(bigsum_data[,4:9],bigsum_data$Yr==yrs[y]);Alt7[y,4:9]<-colMeans(bigsum_data)}
bigsum<-rbind(Alt1,Alt2,Alt3,Alt4,Alt5,Alt6,Alt7)

# Write results from max scenario
write.table(bigfile_write,file=here::here(year, "proj", "author_f", "bigfile.out"),quote=F,row.names=F,col.names=F)
write.table(F_profile,file=here::here(year, "proj", "author_f", "F_profile.out"),quote=F,row.names=F,col.names=F)
write.table(means,file=here::here(year, "proj", "author_f", "means.out"),quote=F,row.names=F,col.names=F)
write.table(percentdb,file=here::here(year, "proj", "author_f", "percentdb.out"),quote=F,row.names=F,col.names=F)
write.table(percentiles,file=here::here(year, "proj", "author_f", "percentiles.out"),quote=F,row.names=F,col.names=F)
write.table(bigsum,file=here::here(year, "proj", "author_f", "bigsum.dat"),quote=F,row.names=F,col.names=c("Alt","Stock","Year","ABC","OFL","Catch","SSB","F","Total_Biom"))


# data ----
setup(year)

# plot survey results ----
vast <- vroom(here("2021", "data", "user_input", "2021_vast_survey.csv"))
sb <- vroom(here("2021", "data", "output",  "goa_ts_biomass_db.csv"))

vast |> 
  mutate(lci = t - sd * 1.96,
         uci = t + sd * 1.96) |> 
  dplyr::select(-sd) |> 
  mutate(Model = "VAST") |> 
  bind_rows(sb |> 
              rename(t = biomass) |> 
              mutate(Model = "Design-based")) |> 
  group_by(Model) %>% 
  mutate(mean = mean(t)) %>% 
  ungroup() %>% 
  ggplot(aes(year, t, fill = Model, color = Model)) + 
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = tickr(vast, year)$breaks,
                     labels = tickr(vast, year)$labels) +
  scale_y_continuous(labels = scales::comma) +
  ylab("Metric tons\n") +
  xlab("\nYear") +
  expand_limits(y = 0) +
  scale_fill_scico_d(palette = "roma", begin = 0.2) +
  scale_color_scico_d(palette = "roma", begin = 0.2) +
  geom_line(aes(y = mean), lty = 3) +
  # theme_present() +
  theme(legend.position = c(0.2, 0.8))

ggsave(here::here(year, "figs", "ts_biomass.png"), width = 6.5, height = 5.5, units = "in", dpi = 200)
# ggsave(here::here(year, "figs", "ts_biomass_pres.png"), width = 6.5, height = 5.5, units = "in", dpi = 200)
# plot catch/biomass
rep <-  readLines(here::here(here::here(year, "base", "base.rep")))

std <- read.delim(here::here(here::here(year, "base", "base.std")), sep="", header = TRUE)

filter(catch, Year == year) %>%
  left_join(read.delim(here::here(here::here(year, "proj", "author_f", "bigsum.dat")), sep="", header = TRUE) %>%
  filter(Year == year, Alt == 2) %>%
  dplyr::select(Year, value = Total_Biom))

std %>% 
  filter(name=="tot_biom") %>%
  bind_cols(filter(catch, Year < year)) %>% 
  filter(Year >= 1991) %>%
  dplyr::select(Year, Catch, value, std.dev) %>%
  bind_rows(filter(catch, Year == year) %>%
              left_join(read.delim(here::here(here::here(year, "proj", "author_f", "bigsum.dat")), sep="", header = TRUE) %>%
                          filter(Year == year, Alt == 2) %>%
                          mutate(value = Total_Biom * 1000) %>%
                          dplyr::select(Year, value))) %>%
  mutate(std.dev = ifelse(is.na(std.dev), std.dev[Year==year-1], std.dev)) %>% 
  mutate(lci = value - std.dev * 1.96,
         uci = value + std.dev * 1.96) %>% 
  mutate(perc = Catch / value,
         lci = Catch / lci,
         uci = Catch / uci,
         mean = mean(perc)) %>% 
  dplyr::select(Year, value, mean, perc, lci, uci) -> df 

df %>% 
  summarise(min = min(perc),
            max = max(perc))
(0.03065468-0.02214721)/0.02214721

df %>% 
  ggplot(aes(Year, perc)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_hline(yintercept = df$mean, lty = 3) +
  expand_limits(y = c(0, 0.08)) +
  scale_x_continuous(breaks = funcr::tickr(df, Year, start = 1990)$breaks,
                     labels = funcr::tickr(df, Year, start = 1990)$labels) +
  # theme_present() +
  xlab("\nYear") +
  ylab("Catch/Biomass\n") 


ggsave(here::here(year, "figs", "catch_bio.png"), width = 6.5, height = 5.5, units = "in", dpi = 200)
# ggsave(here::here(year, "figs", "catch_bio_pres.png"), width = 6.5, height = 5.5, units = "in", dpi = 200)
