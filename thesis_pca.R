
setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")
# Load required library for PCA
library(stats)
library(dplyr)
library(vegan)
library(PERMANOVA)
library(iNEXT)
library(tidyverse)
library(tidyr)
library(ggplot2)


# Assuming your data is stored in a matrix called "bird_data"
# Replace this with your actual data
bird_data <- read.csv("Namdapha_Flock_Data.csv", header = T)
flock_csv = read.csv("flocktypes.csv")

FlockType = flock_csv$flock_type

# Perform PCA
pca_result <- prcomp(bird_data, scale. = TRUE)

# Summary of PCA
summary(pca_result)

# Biplot of the PCA
biplot(pca_result, scale = .1)

# Scree plot to visualize the proportion of variance explained by each principal component
plot(pca_result)

# Eigenvalues
pca_result$sdev^2

# Proportion of variance explained by each principal component
prop_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)
prop_var

#GGPLOT CODE FROM CLASS!!!!!!!
nmds = metaMDS(bird_data, distance = "bray", k = 2, maxit = 999)
nmds_dat = as_tibble(nmds$points)
nmds_dat = nmds_dat %>% mutate(FlockType = as.factor(FlockType))

ordiellipse(nmds, groups = FlockType, col = c("#0387ad", "#02a890", "#83c273"), 
            kind = "se", conf = 0.95, draw = "lines")

#ggplot code for NMDS
ggplot(nmds_dat, aes(x = MDS1, y = MDS2, color = FlockType)) +
  coord_cartesian(xlim = c(-2,2.5), ylim = c(-2,2)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = FlockType), type = "t", linewidth = 1, level = 0.95) +
  theme_classic() +
  scale_color_manual(values = c("#027d68", "#4682B4", "#CD5B45")) +
  scale_fill_manual(values = c("#027d68", "#4682B4", "#CD5B45")) +
  geom_point(size = 1.5, alpha = .9)
# Add + geom_text(size = 3, hjust = 0.5, vjust = 0.5, label = flock_csv$flock_ID) for labels

nmds$stress
adonis2()

