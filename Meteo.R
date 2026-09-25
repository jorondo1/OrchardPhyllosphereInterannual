#Script to calculate and plot meteo data



# Load packages -----------------------------------------------------------

library(tidyverse)
library(scales)



# Define paths ------------------------------------------------------------

path2 <- paste0("S:/LaforestLapointeI/ANJA_WERZ/PhD/A_Apple_microbiome/Data/Meteo/Horaire/")



# Load data ---------------------------------------------------------------
#Horaire, Sherbrooke
meteo_0522 <- read.csv(paste0(path2, "Meteo_052022.csv"), header = T)

meteo_0722 <- read.csv(paste0(path2, "Meteo_072022.csv"), header = T)

meteo_0523 <- read.csv(paste0(path2, "Meteo_052023.csv"), header = T)

meteo_0723 <- read.csv(paste0(path2, "Meteo_072023.csv"), header = T)

meteo_0524 <- read.csv(paste0(path2, "Meteo_052024.csv"), header = T)

meteo_0724 <- read.csv(paste0(path2, "Meteo_072024.csv"), header = T)


#Horaire, individual stations
meteo_ASB_0522 <- read.csv(paste0(path2, "Meteo_ASB_052022.csv"), header = T)

meteo_ASB_0523 <- read.csv(paste0(path2, "Meteo_ASB_052023.csv"), header = T)

meteo_ASB_0524 <- read.csv(paste0(path2, "Meteo_ASB_052024.csv"), header = T)

meteo_ASB_0722 <- read.csv(paste0(path2, "Meteo_ASB_072022.csv"), header = T)

meteo_ASB_0723 <- read.csv(paste0(path2, "Meteo_ASB_072023.csv"), header = T)

meteo_ASB_0724 <- read.csv(paste0(path2, "Meteo_ASB_072024.csv"), header = T)


meteo_COMP_0522 <- read.csv(paste0(path2, "Meteo_COMP_052022.csv"), header = T)

meteo_COMP_0523 <- read.csv(paste0(path2, "Meteo_COMP_052023.csv"), header = T)

meteo_COMP_0524 <- read.csv(paste0(path2, "Meteo_COMP_052024.csv"), header = T)

meteo_COMP_0722 <- read.csv(paste0(path2, "Meteo_COMP_072022.csv"), header = T)

meteo_COMP_0723 <- read.csv(paste0(path2, "Meteo_COMP_072023.csv"), header = T)

meteo_COMP_0724 <- read.csv(paste0(path2, "Meteo_COMP_072024.csv"), header = T)


meteo_MILT_0522 <- read.csv(paste0(path2, "Meteo_MILT_052022.csv"), header = T)

meteo_MILT_0523 <- read.csv(paste0(path2, "Meteo_MILT_052023.csv"), header = T)

meteo_MILT_0524 <- read.csv(paste0(path2, "Meteo_MILT_052024.csv"), header = T)

meteo_MILT_0722 <- read.csv(paste0(path2, "Meteo_MILT_072022.csv"), header = T)

meteo_MILT_0723 <- read.csv(paste0(path2, "Meteo_MILT_072023.csv"), header = T)

meteo_MILT_0724 <- read.csv(paste0(path2, "Meteo_MILT_072024.csv"), header = T)



# Combine and filter dataset (Horaire, closest station) ----------------------------------------------

meteo_hi <- rbind(meteo_ASB_0522, meteo_ASB_0523) %>%
  rbind(meteo_ASB_0524) %>%
  rbind(meteo_ASB_0722) %>%
  rbind(meteo_ASB_0723) %>%
  rbind(meteo_ASB_0724) %>%
  rbind(meteo_COMP_0522) %>%
  rbind(meteo_COMP_0523) %>%
  rbind(meteo_COMP_0524) %>%
  rbind(meteo_COMP_0722) %>%
  rbind(meteo_COMP_0723) %>%
  rbind(meteo_COMP_0724) %>%
  rbind(meteo_MILT_0522) %>%
  rbind(meteo_MILT_0523) %>%
  rbind(meteo_MILT_0524) %>%
  rbind(meteo_MILT_0722) %>%
  rbind(meteo_MILT_0723) %>%
  rbind(meteo_MILT_0724) %>%
  rbind(meteo_0522) %>%
  rbind(meteo_0722) %>%
  rbind(meteo_0523) %>%
  rbind(meteo_0723) %>%
  rbind(meteo_0524) %>%
  rbind(meteo_0724) %>%
  select("Nom.de.la.Station", Année, Mois, Jour, "Heure..HNL.", Temp,  Hum_rel, Hauteur_precip) %>%
  mutate(Temp = gsub(",", ".", Temp), Hauteur_precip = gsub(",", ".", Hauteur_precip)) %>%
  mutate(Temp = as.numeric(Temp), Hauteur_precip = as.numeric(Hauteur_precip)) %>%
  rename(Year = Année, Heure = "Heure..HNL.", Station = "Nom.de.la.Station") %>%
  mutate(Date = paste0(Jour, "-", Mois, "-", Year)) %>%
  mutate(Date = as.Date(Date, format = "%d-%m-%Y"), Year = factor(Year, levels = c(2022, 2023, 2024)), Date2 = paste0(Jour, "-", Mois),
         Mois2 = ifelse(Mois == 05, paste0("May"), paste0("July")), Mois2 = factor(Mois2, levels = c("May", "July")),
         Temp = ifelse(Temp == "", NA, Temp))



# Calculate meteo data (Horaire, closest meteo station) ---------------------------------

PMB_EIP3 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(13, 14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(7, 8, 9, 10) |
           Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(13, 14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), deg_h = sum(Temp[Temp > 18.3], na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  ungroup() %>%
  select(-Year, -Mois)

PMB_EIP4 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(8, 9, 10) |
           Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("PMB")) %>%
  ungroup() %>%
  select(-Year, -Mois)

PMB_EIP_i <- merge(PMB_EIP3, PMB_EIP4, by = "Time")


ASB_EIP3 <- meteo_hi %>%
  filter(Station == "LAC MEMPHREMAGOG") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(13, 14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(3, 4, 5, 6, 7) |
           Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(15, 16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), deg_h = sum(Temp[Temp > 18.3], na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  ungroup() %>%
  select(-Year, -Mois)

ASB_EIP4 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(4, 5, 6) |
           Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("ASB")) %>%
  ungroup() %>%
  select(-Year, -Mois)

ASB_EIP_i <- merge(ASB_EIP3, ASB_EIP4, by = "Time")


VBS_EIP3 <- meteo_hi %>%
  filter(Station == "SHERBROOKE") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(16, 17, 18, 19) |
           Year == 2022 & Mois == 7 & Jour %in% c(3, 4, 5, 6) |
           Year == 2023 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2023 & Mois == 7 & Jour %in% c(15, 16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), deg_h = sum(Temp[Temp > 18.3], na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  ungroup() %>%
  select(-Year, -Mois)

VBS_EIP4 <- meteo_hi %>%
  filter(Station == "SHERBROOKE") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(17, 18, 19) |
           Year == 2022 & Mois == 7 & Jour %in% c(4, 5, 6) |
           Year == 2023 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2023 & Mois == 7 & Jour %in% c(16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("VBS")) %>%
  ungroup() %>%
  select(-Year, -Mois)

VBS_EIP_i <- merge(VBS_EIP3, VBS_EIP4, by = "Time")


COM_EIP3 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(13, 14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), deg_h = sum(Temp[Temp > 18.3], na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  ungroup() %>%
  select(-Year, -Mois)

COM_EIP4 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("COM")) %>%
  ungroup() %>%
  select(-Year, -Mois)

COM_EIP_i <- merge(COM_EIP3, COM_EIP4, by = "Time")


MIB_EIP3 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(21, 22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(12, 13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), deg_h = sum(Temp[Temp > 18.3], na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  ungroup() %>%
  select(-Year, -Mois)

MIB_EIP4 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("MIB")) %>%
  ungroup() %>%
  select(-Year, -Mois)

MIB_EIP_i <- merge(MIB_EIP3, MIB_EIP4, by = "Time")


MIC_EIP3 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(21, 22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(12, 13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), deg_h = sum(Temp[Temp > 18.3], na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  ungroup() %>%
  select(-Year, -Mois)

MIC_EIP4 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T)) %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("MIC")) %>%
  ungroup() %>%
  select(-Year, -Mois)

MIC_EIP_i <- merge(MIC_EIP3, MIC_EIP4, by = "Time")


meteo_ind_i <- rbind(PMB_EIP_i, ASB_EIP_i, VBS_EIP_i, COM_EIP_i, MIC_EIP_i, MIB_EIP_i) %>%
  mutate(Location = case_when(Site %in% c("PMB", "COM") ~ "B", Site %in% c("MIC", "MIB") ~ "D", Site == "ASB" ~ "A", Site == "VBS" ~ "C"),
         Time = case_when(Time == "5_2022" ~ "May 2022", Time == "5_2023" ~ "May 2023", Time == "5_2024" ~ "May 2024",
                          Time == "7_2022" ~ "July 2022", Time == "7_2023" ~ "July 2023", Time == "7_2024" ~ "July 2024"),
         Time = factor(Time, levels = c("May 2022", "May 2023", "May 2024", "July 2022", "July 2023", "July 2024")))

#write_csv(meteo_ind_i, "Meteo_indexes.csv")


###Supp. figure###
ggplot(meteo_ind_i, aes(x = Time, y = mean_temp, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  labs(x = "", y = "Mean temperature [°C] (last 96h)") +
  scale_fill_manual(values = fill_loc, breaks = labels_loc) +
  theme_barplot

ggsave("Meteo_mean_temp_ind_96h.png", dpi = 300)


###Supp. figure###
ggplot(meteo_ind_i, aes(x = Time, y = deg_h, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  labs(x = "", y = "Degree hours >18.3 °C (last 96h)") +
  scale_fill_manual(values = fill_loc, breaks = labels_loc) +
  theme_barplot

ggsave("Meteo_deg_h_ind_96h.png", dpi = 300)


###Supp. figure###
ggplot(meteo_ind_i, aes(x = Time, y = precip_72h, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")  +
  labs(x = "", y = "Precipitation [mm] (last 72h)") +
  scale_fill_manual(values = fill_loc, breaks = labels_loc) +
  theme_barplot

ggsave("Meteo_precip_i_72h.png", dpi = 300)

