library(DBI)
library(tidyverse)
library(viridis)


library(ape); #The Swiss army knife of comparative phylogenetics R packages
library(geiger); #A package for "macroevolutionary simulation and estimating parameters related to diversification"
library(phytools); #Very useful for visualization particularly, great blog support
library(RColorBrewer); # Accessory package with better color
library(plotrix); #Contains a useful way of producing a color legend
library(ggplot2)
library(Biostrings)
library(dplyr)
library(ggtree)
library(treeio)
library(phyloseq)
library(ggjoy)
library(tidytree)
library(RCurl)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ggtree")

if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")




mltree <- read.tree("final_iqMF.merge.contree")

to_drop<-c("QCAZ27653", "AES2439", "MZUA2392", "MZUA785",
           "MZUA2366", "MZUA3012", "MZUA3012", "QCAZ61154", "MZUA2514",
           "MZUA1763", "VLU397", "MECN13062")

ml_reduced<- drop.tip(mltree, to_drop)

ggtree(ml_reduced, , ladderize = FALSE)  + geom_text(aes(label=node)) + geom_tiplab() 

tree_orestes <- extract.clade(ml_reduced, node = 423)

ggtree(tree_orestes, ladderize = FALSE, size = 1) +
  geom_text(aes(label=node))
  geom_tiplab(size = 6) + 
  geom_text(aes(label=node))
  geom_cladelabel(node=69, label="P. sp1", color="gold",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=67, label="P. sp2", color="black",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=64, label="P. colodactylus", color = "gold",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=62, label="P. sp3", color = "black",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=60, label="P. sp4", color = "gold",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=59, label="P. sp5", color = "black",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=58, label="P. sp6", color = "gold",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
    geom_cladelabel(node=54, label="P. muranunka", color = "gold",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=42, label="P. sp7", color = "black",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3) +
  geom_cladelabel(node=41, label="P. sp8", color = "gold",
                  offset= 0.012, , align=TRUE, fontsize = 5, barsize = 3)
  
  
  
  habitat <- read.csv("4_habitat.csv", fileEncoding="UTF-8-BOM", 
                      stringsAsFactors = FALSE)
  
  records <- read.csv("2_records_new.csv", fileEncoding="UTF-8-BOM", 
                      stringsAsFactors = FALSE) 
  
  pristimantis <- read.csv("1_pristimantis.csv", fileEncoding="UTF-8-BOM", 
                           stringsAsFactors = FALSE) 
  morphometrics <- read.csv("5_morphometrics.csv", fileEncoding="UTF-8-BOM", 
                            stringsAsFactors = FALSE) 
  
  data2 <- habitat %>% 
    left_join(pristimantis, by = "museum_id") %>% 
    left_join(records, by = "museum_id") %>% 
    group_by(ecosystem) %>% 
    summarise(species_number = n_distinct(species), elevation =) 
    
   
  names(data2)
  
  graph4.1 <- ggplot(data = data2, aes(x = ecosystem,
                                       y = species_number, fill = ecosystem, na.rm = TRUE)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    labs(x = "Ecosystem" , y = "number of species") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_viridis_d()
  
  graph4.1
  
  
  data3 <- pristimantis %>% 
    left_join(morphometrics, by = "museum_id") %>% 
    left_join(records, by = "museum_id") %>% 
    group_by(species)   %>% 
    summarise( altitude = mean(elevation)) 
  
  %>% 
    left_join(records, by = "museum_id") %>% 
    mutate(hand = Finger.width_mm/HANDL_mm) %>% 
    group_by(species) %>% 
    summarise(hand_value = mean(hand), altitude = mean(elevation)) 
    

  
  data5 <- pristimantis %>% 
    left_join(records, by = "museum_id") %>% 
    left_join(morphometrics, by = "museum_id") %>% 
    filter(elevation >1) %>% 
    mutate(hand = Finger.width_mm/HANDL_mm) %>% 
    group_by(species) %>% 
    summarise(hand_value = mean(hand), altitude = mean(elevation)) 
  
  data6 <- pristimantis %>% 
    left_join(morphometrics, by = "museum_id") %>% 
    left_join(records, by = "museum_id") %>% 
    filter(HANDL_mm > 0 & Finger.width_mm >0 & elevation >1) %>% 
    mutate(hand = Finger.width_mm/HANDL_mm) %>% 
    group_by(species) %>% 
    summarise(hand_value = mean(hand), altitude = mean(elevation)) 
  
  
  
  m1 <- lm(formula = altitude ~ hand_value, data = data6)
  
  pred <- predict(m1, se.fit = TRUE)
  
  preds <- data.frame(mean = pred$fit,
                      upr = pred$fit + 1.96 * pred$se.fit,
                      lwr = pred$fit - 1.96 * pred$se.fit)
  
  ggplot(data6, aes(x = altitude, y = hand_value, color = altitude)) +
    geom_point() +
    stat_smooth(aes(color= altitude, fill = altitude),
                method = "lm", formula = y ~ x, size = 1, 
                color = "darkmagenta", se=TRUE) +
    labs(x = "Altitude" , y = "FDW/HAL") +
    theme_classic() 
    
  
  
  p2 <- data6 %>% 
    ggplot(mapping = aes(x = hand_value, y = altitude)) +
    geom_ribbon(aes(ymin = preds$lwr, ymax = preds$upr), 
                fill = "gray90") +
    geom_line(aes(y = preds$mean), color = "darkmagenta", size = 0.8) + 
    geom_point(color = "darkmagenta", alpha = 0.5) +
    labs(x = "FDW/HANDL", y = "Altitude") +
    theme_minimal() 
  
p2