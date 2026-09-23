#color schemes and themes


library(ggplot2)


#year
labels_year <- c("2022", "2023", "2024")

colors_year <- c("black", "black", "black")

shapes_year <- c(21,21,21)

fill_year <- c("#3B839E","#BB4D3C", "#9AB88A") #first two colors adopted from Sophie


#time
labels_time <- c("May", "July", "other")

colors_time <- c("black", "black", "grey80")

shapes_time <- c(21, 21, 21)

shapes_time2 <- c(21, 22)

fill_time <- c("#CC8FBB", "#7AAB32", "grey90") #colors adopted from Sophie


#management
labels_man <- c("Conventional", "Organic", "other")

colors_man <- c("black", "black", "grey80")

shapes_man <- c(21, 21, 21)

fill_man <- c("#E4A804", "#527798", "grey90") #colors adopted from Sophie


#time x management
labels_timman <- c("May_Conventional", "May_Organic", "July_Conventional", "July_Organic", "other")

colors_timman <- c("black", "black", "black", "black", "grey80")

shapes_timman <- c(21, 21, 21, 21)

fill_timman <- c("#D89C60", "#8F83AA", "#AFAA1B", "#669165", "grey90") #https://colorkit.co/color-mixer/?mix=7aab32-e4a804&ratios=1-1&space=rgb


#site
labels_site <- c("A", "B1", "B2", "C", "D1", "D2")

colors_site <- c("black", "black", "black", "black", "black", "black")

shapes_site <- c(22,21,21,22,23,23)

fill_site <- c("#D17913", "#F3A44A","#096EA4","#1C9EE4","#FFC787","#89CEF3") #colors adopted from Sophie


#location
labels_loc <- c("A", "B", "C", "D")

colors_loc <- c("black", "black", "black", "black")

shapes_loc <- c(21, 21, 21, 21)

fill_loc <- c("#5DB63B", "#2C9EE3", "gold", "#F8A11C")


#cultivar
labels_cult <- c("Cortland", "Liberty", "Paulared", "Honeycrisp", "Spartan")

colors_cult <- c("black", "black", "black", "black", "black")

shapes_cult <- c(21, 21, 21, 21, 21)

fill_cult <- c("#7DB16B","#45818E","#AB4F84", "#AE9FCB", "#99CFE1") #first three colors adopted from Sophie


#gradient
grad_low <- c("red")

grad_mid <- c("white")

grad_hi <- c("blue")

grad_na <- c("white")


#families --> Amy
labels_fam <- c("Acetobacteraceae", "Bacillaceae", "Beijerinckiaceae", "Comamonadaceae",
                "Deinococcaceae", "Enterobacteriaceae", "Erwiniaceae", "Geodermatophilaceae",
                "Hymenobacteraceae", "Lactobacillaceae", "Microbacteriaceae", "Micrococcaceae",
                "Nocardioidaceae", "Oxalobacteraceae", "Peptostreptococcaceae", "Pseudomonadaceae",
                "Pseudonocardiaceae", "Roseiflexaceae", "Sphingomonadaceae", "Spirosomataceae", "Weeksellaceae",
                "Botryosphaeriaceae", "Buckleyzymaceae", "Bulleraceae", "Bulleribasidiaceae",
                "Cladosporiaceae", "Cryptococcaceae", "Didymellaceae", "Didymosphaeriaceae",
                "Erysiphaceae", "Erythrobasidiaceae", "Filobasidiaceae", "Intrasporangiaceae",
                "Kineosporiaceae", "Mycosphaerellaceae", "Nectriaceae", "Phaeosphaeriaceae",
                "Polyporaceae", "Pseudeurotiaceae", "Saccharimonadaceae", "Saccotheciaceae",
                "Sclerotiniaceae", "Sporidiobolaceae", "Sporocadaceae", "Taphrinaceae", "Venturiaceae",
                "Others", "Unclassified")

colors_fam <- c(rep("black"), length(labels_fam))

fill_fam <- c("#F1F42C", "#F4D801", "#F8B11B", "#E17919",
              "#E13E11", "#711939", "#C85A78", "#DA1C91",
              "#AA9ADD", "#AA5ADD", "#A16591", "#A1C991",  
              "#A3B78C", "#A3E18C", "#37A346", "#43BCCF", 
              "#118998FF", "#1F677A", "#1F549A", "#4A559A", "#244154",
              "#F6F199", "#F7E42C", "#E8D101", "#E8B11B",
              "#F9B36D", "#D17919", "#D13E11", "#611939", 
              "#B85A78", "#E35959", "#CA1C91", "#D5AEE6",   
              "#A99ADD", "#A95ADD", "#A19591", "#A1B991", 
              "#A3A78C", "#A3D18C", "#379346", "#43ACCF",
              "#117998FF", "#1E677A", "#1E549A", "#3A559A", "#144154",
              "grey40", "grey90")




#themes
theme_pcoa <- theme(plot.title = element_text(size = 18), 
                    legend.position = "right",
                    legend.title = element_text(size=18),
                    legend.text = element_text(size=16),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    strip.text = element_text(size = 16),
                    strip.background = element_rect(color = "black", fill = "white"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(color = "white", fill = "white"),
                    panel.border = element_rect(color = "black"),
                    plot.background = element_rect(color = "white", fill = "white")
                    )


theme_barplot <- theme(legend.title = element_text(size=18),
                       legend.text = element_text(size=16, face = "italic"),
                       plot.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       axis.text = element_text(size = 14),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_rect(color = "white", fill = "white"),
                       panel.border = element_rect(color = "black"),
                       plot.background = element_rect(color = "white", fill = "white"),
                       legend.key.height = grid::unit(0.45, "cm"))












