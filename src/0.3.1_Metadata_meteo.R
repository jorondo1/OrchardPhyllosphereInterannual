# Script to plot temperature and precipitation of May and July 2022, 2023, 2024
# Author: Anja Werz


# Load packages -----------------------------------------------------------

library(tidyverse)
library(scales)


# Define paths ------------------------------------------------------------

path1 <- paste0("data/Meteo/Quotidien")

path2 <- paste0("data/Meteo/Horaire")


# Load data ---------------------------------------------------------------

#Quotidien
meteo_22 <- read.csv(file.path(path1, "Meteo_2022.csv"), header = T)

meteo_23 <- read.csv(file.path(path1, "Meteo_2023.csv"), header = T)

meteo_24 <- read.csv(file.path(path1, "Meteo_2024.csv"), header = T)


#Horaire, Sherbrooke
meteo_0522 <- read.csv(file.path(path2, "Meteo_052022.csv"), header = T)

meteo_0722 <- read.csv(file.path(path2, "Meteo_072022.csv"), header = T)

meteo_0523 <- read.csv(file.path(path2, "Meteo_052023.csv"), header = T)

meteo_0723 <- read.csv(file.path(path2, "Meteo_072023.csv"), header = T)

meteo_0524 <- read.csv(file.path(path2, "Meteo_052024.csv"), header = T)

meteo_0724 <- read.csv(file.path(path2, "Meteo_072024.csv"), header = T)


#Horaire, individual stations
meteo_ASB_0522 <- read.csv(file.path(path2, "Meteo_ASB_052022.csv"), header = T)

meteo_ASB_0523 <- read.csv(file.path(path2, "Meteo_ASB_052023.csv"), header = T)

meteo_ASB_0524 <- read.csv(file.path(path2, "Meteo_ASB_052024.csv"), header = T)

meteo_ASB_0722 <- read.csv(file.path(path2, "Meteo_ASB_072022.csv"), header = T)

meteo_ASB_0723 <- read.csv(file.path(path2, "Meteo_ASB_072023.csv"), header = T)

meteo_ASB_0724 <- read.csv(file.path(path2, "Meteo_ASB_072024.csv"), header = T)


meteo_COMP_0522 <- read.csv(file.path(path2, "Meteo_COMP_052022.csv"), header = T)

meteo_COMP_0523 <- read.csv(file.path(path2, "Meteo_COMP_052023.csv"), header = T)

meteo_COMP_0524 <- read.csv(file.path(path2, "Meteo_COMP_052024.csv"), header = T)

meteo_COMP_0722 <- read.csv(file.path(path2, "Meteo_COMP_072022.csv"), header = T)

meteo_COMP_0723 <- read.csv(file.path(path2, "Meteo_COMP_072023.csv"), header = T)

meteo_COMP_0724 <- read.csv(file.path(path2, "Meteo_COMP_072024.csv"), header = T)


meteo_MILT_0522 <- read.csv(file.path(path2, "Meteo_MILT_052022.csv"), header = T)

meteo_MILT_0523 <- read.csv(file.path(path2, "Meteo_MILT_052023.csv"), header = T)

meteo_MILT_0524 <- read.csv(file.path(path2, "Meteo_MILT_052024.csv"), header = T)

meteo_MILT_0722 <- read.csv(file.path(path2, "Meteo_MILT_072022.csv"), header = T)

meteo_MILT_0723 <- read.csv(file.path(path2, "Meteo_MILT_072023.csv"), header = T)

meteo_MILT_0724 <- read.csv(file.path(path2, "Meteo_MILT_072024.csv"), header = T)



# Combine and filter dataset (Quotidien) ----------------------------------------------

meteo_q <- rbind(meteo_22, meteo_23) %>%
  rbind(meteo_24) %>%
  filter(Mois %in% c(5, 7)) %>%
  select(Date.Heure, Année, Mois, Jour, Temp_max, Temp_min, Temp_moy, Precip_tot) %>%
  mutate(Temp_max = gsub(",", ".", Temp_max), Temp_min = gsub(",", ".", Temp_min), 
         Temp_moy = gsub(",", ".", Temp_moy), Precip_tot = gsub(",", ".", Precip_tot)) %>%
  mutate(Temp_max = as.numeric(Temp_max), Temp_min = as.numeric(Temp_min), Temp_moy = as.numeric(Temp_moy), Precip_tot = as.numeric(Precip_tot)) %>%
  rename(Date = Date.Heure, Year = Année) %>%
  mutate(Date = as.Date(Date, format = "%d.%m.%Y"), Year = factor(Year, levels = c(2022, 2023, 2024)), Date2 = paste0(Jour, "-", Mois),
         Mois2 = ifelse(Mois == 5, paste0("May"), paste0("July")), Mois2 = factor(Mois2, levels = c("May", "July")))



# Plot data (Quotidien) ---------------------------------------------------------------

meteo_q %>%
  filter(Mois2 == "May") %>%
ggplot(aes(x = reorder(Date2, Date), y = Temp_moy, colour = Year, group = Year)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = c("red3", "royalblue", "forestgreen"), breaks = c("2022", "2023", "2024")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("rect", xmin = "17-5", xmax = "20-5", ymin = 0, ymax = 25, alpha = 0.2, fill = "red3") +
  annotate("rect", xmin = "15-5", xmax = "25-5", ymin = 0, ymax = 25, alpha = 0.2, fill = "royalblue") +
  annotate("rect", xmin = "21-5", xmax = "23-5", ymin = 0, ymax = 25, alpha = 0.2, fill = "forestgreen") +
  annotate("text", x = 18.5, y = -0.5, label = "Sampling 2022", size = 4) +
  annotate("text", x = 18.5, y = 0.5, label = "Sampling 2023", size = 4) +
  annotate("text", x = 18.5, y = 1.5, label = "Sampling 2024", size = 4) +
  annotate("segment", x = 17, xend = 20, y = 0, yend = 0, color = "black") +
  annotate("segment", x = 15, xend = 23, y = 1, yend = 1, color = "black") +
  annotate("segment", x = 17, xend = 21, y = 2, yend = 2, color = "black")

ggsave("out/meteo/Meteo_temp_may.png", dpi = 300)



meteo_q %>%
  filter(Mois2 == "July") %>%
  ggplot(aes(x = reorder(Date2, Date), y = Temp_moy, colour = Year, group = Year)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = c("red3", "royalblue", "forestgreen"), breaks = c("2022", "2023", "2024")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("rect", xmin = "7-7", xmax = "11-7", ymin = 0, ymax = 25, alpha = 0.2, fill = "red3") +
  annotate("rect", xmin = "17-7", xmax = "21-7", ymin = 0, ymax = 25, alpha = 0.2, fill = "royalblue") +
  annotate("text", x = 9, y = -0.5, label = "Sampling 2022", size = 4) +
  annotate("text", x = 19, y = -0.5, label = "Sampling 2023", size = 4) +
  annotate("segment", x = 7, xend = 11, y = 0, yend = 0, color = "black") +
  annotate("segment", x = 17, xend = 21, y = 0, yend = 0, color = "black")

ggsave("out/meteo/Meteo_temp_july.png", dpi = 300)



meteo_q %>%
  filter(Mois2 == "May") %>%
ggplot(aes(x = reorder(Date2, Date), y = Precip_tot, colour = Year, group = Year)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = c("red3", "royalblue", "forestgreen"), breaks = c("2022", "2023", "2024")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("rect", xmin = "17-5", xmax = "20-5", ymin = 0, ymax = 51, alpha = 0.2, fill = "red3") +
  annotate("rect", xmin = "15-5", xmax = "25-5", ymin = 0, ymax = 51, alpha = 0.2, fill = "royalblue") +
  annotate("rect", xmin = "21-5", xmax = "23-5", ymin = 0, ymax = 51, alpha = 0.2, fill = "forestgreen") +
  annotate("text", x = 18.5, y = -0.5, label = "Sampling 2022", size = 4) +
  annotate("text", x = 18.5, y = -1.5, label = "Sampling 2023", size = 4) +
  annotate("text", x = 18.5, y = -2.5, label = "Sampling 2024", size = 4) +
  annotate("segment", x = 17, xend = 20, y = 0, yend = 0, color = "black") +
  annotate("segment", x = 15, xend = 23, y = -1, yend = -1, color = "black") +
  annotate("segment", x = 17, xend = 21, y = -2, yend = -2, color = "black")

ggsave("out/meteo/Meteo_prec_may.png", dpi = 300)



meteo_q %>%
  filter(Mois2 == "July") %>%
  ggplot(aes(x = reorder(Date2, Date), y = Precip_tot, colour = Year, group = Year)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = c("red3", "royalblue", "forestgreen"), breaks = c("2022", "2023", "2024")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("rect", xmin = "7-7", xmax = "11-7", ymin = 0, ymax = 81, alpha = 0.2, fill = "red3") +
  annotate("rect", xmin = "17-7", xmax = "21-7", ymin = 0, ymax = 81, alpha = 0.2, fill = "royalblue") +
  annotate("text", x = 9, y = -0.5, label = "Sampling 2022", size = 4) +
  annotate("text", x = 19, y = -0.5, label = "Sampling 2023", size = 4) +
  annotate("segment", x = 7, xend = 11, y = 0, yend = 0, color = "black") +
  annotate("segment", x = 17, xend = 21, y = 0, yend = 0, color = "black")

ggsave("out/meteo/Meteo_prec_july.png", dpi = 300)



# Combine and filter dataset (Horaire, Sherbrooke) ----------------------------------------------

meteo_h <- rbind(meteo_0522, meteo_0722) %>%
  rbind(meteo_0523) %>%
  rbind(meteo_0723) %>%
  rbind(meteo_0524) %>%
  rbind(meteo_0724) %>%
  select(Année, Mois, Jour, 
         #Heure, 
         Temp,  Hum_rel, Hauteur_precip) %>%
  mutate(Temp = gsub(",", ".", Temp), Hauteur_precip = gsub(",", ".", Hauteur_precip)) %>%
  mutate(Temp = as.numeric(Temp), Hauteur_precip = as.numeric(Hauteur_precip)) %>%
  rename(Year = Année) %>%
  mutate(Date = paste0(Jour, "-", Mois, "-", Year)) %>%
  mutate(Date = as.Date(Date, format = "%d-%m-%Y"), Year = factor(Year, levels = c(2022, 2023, 2024)), Date2 = paste0(Jour, "-", Mois),
         Mois2 = ifelse(Mois == 05, paste0("May"), paste0("July")), Mois2 = factor(Mois2, levels = c("May", "July")),
         Temp = ifelse(Temp == "", NA, Temp))



# Calculate Epiphytic Infection Potential (Sherbrooke meteo station) ---------------------------------

PMB_EIP1 <- meteo_h %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(13, 14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(7, 8, 9, 10) |
           Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(13, 14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  #ungroup() %>%
  select(-Year, -Mois)

PMB_EIP2 <- meteo_h %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(8, 9, 10) |
           Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("PMB")) %>%
  select(-Year, -Mois)

PMB_EIP <- merge(PMB_EIP1, PMB_EIP2, by = "Time")


ASB_EIP1 <- meteo_h %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(13, 14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(3, 4, 5, 6, 7) |
           Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(15, 16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T),
            , .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

ASB_EIP2 <- meteo_h %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(4, 5, 6) |
           Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("ASB")) %>%
  select(-Year, -Mois)

ASB_EIP <- merge(ASB_EIP1, ASB_EIP2, by = "Time")


VBS_EIP1 <- meteo_h %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(16, 17, 18, 19) |
           Year == 2022 & Mois == 7 & Jour %in% c(3, 4, 5, 6) |
           Year == 2023 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2023 & Mois == 7 & Jour %in% c(15, 16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

VBS_EIP2 <- meteo_h %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(17, 18, 19) |
           Year == 2022 & Mois == 7 & Jour %in% c(4, 5, 6) |
           Year == 2023 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2023 & Mois == 7 & Jour %in% c(16, 17, 18) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("VBS")) %>%
  select(-Year, -Mois)

VBS_EIP <- merge(VBS_EIP1, VBS_EIP2, by = "Time")


COM_EIP1 <- meteo_h %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(13, 14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

COM_EIP2 <- meteo_h %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("COM")) %>%
  select(-Year, -Mois)

COM_EIP <- merge(COM_EIP1, COM_EIP2, by = "Time")


MIB_EIP1 <- meteo_h %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(21, 22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(12, 13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T),
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

MIB_EIP2 <- meteo_h %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("MIB")) %>%
  select(-Year, -Mois)

MIB_EIP <- merge(MIB_EIP1, MIB_EIP2, by = "Time")


MIC_EIP1 <- meteo_h %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(21, 22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(12, 13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T),
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

MIC_EIP2 <- meteo_h %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("MIC")) %>%
  select(-Year, -Mois)

MIC_EIP <- merge(MIC_EIP1, MIC_EIP2, by = "Time")


meteo_ind <- rbind(PMB_EIP, ASB_EIP, VBS_EIP, COM_EIP, MIC_EIP, MIB_EIP) %>%
  mutate(Location = case_when(
    Site %in% c("PMB", "COM") ~ "Compton", 
    Site %in% c("MIC", "MIB") ~ "Milton", 
    Site == "ASB" ~ "Saint-Benoît", 
    Site == "VBS" ~ "Windsor"))


ggplot(meteo_ind, aes(x = Time, y = mean_temp, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")

ggsave("out/meteo/Meteo_mean_temp_96h.png", dpi = 300)


ggplot(meteo_ind, aes(x = Time, y = deg_h, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")

ggsave("out/meteo/Meteo_deg_h_96h.png", dpi = 300)


ggplot(meteo_ind, aes(x = Time, y = precip_72h, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")

ggsave("out/meteo/Meteo_precip_72h.png", dpi = 300)


saveRDS(meteo_ind, "data/meteo/processed/Meteo_indicators.rds")


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



# Calculate Epiphytic Infection Potential (Horaire, closest meteo station) ---------------------------------

PMB_EIP3 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2022 & Mois == 5 & Jour %in% c(13, 14, 15, 16) |
           Year == 2022 & Mois == 7 & Jour %in% c(7, 8, 9, 10) |
           Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(13, 14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
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
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("PMB")) %>%
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
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
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
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("ASB")) %>%
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
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
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
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("VBS")) %>%
  select(-Year, -Mois)

VBS_EIP_i <- merge(VBS_EIP3, VBS_EIP4, by = "Time")


COM_EIP3 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(11, 12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(13, 14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(18, 19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

COM_EIP4 <- meteo_hi %>%
  filter(Station == "LENNOXVILLE") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(12, 13, 14) |
           Year == 2023 & Mois == 7 & Jour %in% c(14, 15, 16) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21) |
           Year == 2024 & Mois == 7 & Jour %in% c(14, 15, 16)) %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("COM")) %>%
  select(-Year, -Mois)

COM_EIP_i <- merge(COM_EIP3, COM_EIP4, by = "Time")


MIB_EIP3 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(21, 22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(12, 13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T), 
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

MIB_EIP4 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("MIB")) %>%
  select(-Year, -Mois)

MIB_EIP_i <- merge(MIB_EIP3, MIB_EIP4, by = "Time")


MIC_EIP3 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(21, 22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(17, 18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(19, 20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(12, 13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(mean_temp = mean(Temp, na.rm = T), 
            deg_h = sum(Temp[Temp > 18.3], na.rm = T),
            .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year)) %>%
  select(-Year, -Mois)

MIC_EIP4 <- meteo_hi %>%
  filter(Station == "GRANBY") %>%
  filter(Year == 2023 & Mois == 5 & Jour %in% c(22, 23, 24) |
           Year == 2023 & Mois == 7 & Jour %in% c(18, 19, 20) |
           Year == 2024 & Mois == 5 & Jour %in% c(20, 21, 22) |
           Year == 2024 & Mois == 7 & Jour %in% c(13, 14, 15))  %>%
  group_by(Year, Mois) %>%
  summarise(precip_72h = sum(Hauteur_precip, na.rm = T), .groups = 'drop') %>%
  mutate(Time = paste0(Mois, "_", Year), Site = paste0("MIC")) %>%
  select(-Year, -Mois)

MIC_EIP_i <- merge(MIC_EIP3, MIC_EIP4, by = "Time")


meteo_ind_i <- rbind(PMB_EIP_i, ASB_EIP_i, VBS_EIP_i, COM_EIP_i, MIC_EIP_i, MIB_EIP_i) %>%
  mutate(Location = case_when(
    Site %in% c("PMB", "COM") ~ "Compton", 
    Site %in% c("MIC", "MIB") ~ "Milton",
    Site == "ASB" ~ "Saint-Benoît", 
    Site == "VBS" ~ "Windsor"))

saveRDS(meteo_ind_i, "data/meteo/processed/Meteo_indexes.rds")


ggplot(meteo_ind_i, aes(x = Time, y = mean_temp, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")

ggsave("out/meteo/Meteo_mean_temp_ind_96h.png", dpi = 300)


ggplot(meteo_ind_i, aes(x = Time, y = deg_h, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")

ggsave("out/meteo/Meteo_deg_h_ind_96h.png", dpi = 300)


ggplot(meteo_ind_i, aes(x = Time, y = precip_72h, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")

ggsave("out/meteo/Meteo_precip_i_72h.png", dpi = 300)




