.libPaths("C:\\Users\\ronit\\OneDrive\\Documents\\R\\win-library\\4.1")
setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")

library(vegan)
library(PERMANOVA)
library(iNEXT)

paradiseus <- bird_data %>% filter(dicrurus_paradiseus == 1)
remifer <- bird_data %>% filter(dicrurus_remifer == 1)

lrtd = bird_data %>% filter(dicrurus_remifer == 1)
grtd = bird_data %>% filter(dicrurus_paradiseus == 1)


drongos <- bird_data %>%
  filter((dicrurus_paradiseus == 1 & dicrurus_remifer == 0) |
           (dicrurus_paradiseus == 0 & dicrurus_remifer == 1))


species_presence <- drongos %>%
  mutate(species = case_when(dicrurus_paradiseus == 1 ~ "Dicrurus paradiseus", 
                             dicrurus_remifer == 1 ~ "Dicrurus remifer"))



drongo_pca = prcomp(drongos)
drongo_nmds = metaMDS(drongos, distance = "jaccard", k = 2, maxit = 999)
drongo_dat = as_tibble(drongo_nmds$points)
drongo_dat = drongo_dat %>% mutate(levels = as.factor(species_presence$species))

# NMDS for the 2 drongo flocks
ggplot(drongo_dat, aes(x = MDS1, y = MDS2, color = levels)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = levels), type = "t", linewidth = 1.2, level = 0.95) +
  theme_classic() + geom_point(size = 1.5, alpha = .9) +
  scale_color_manual(values = c("#4682B4", "#CD5B45")) +
  scale_fill_manual(values = c("#4682B4", "#CD5B45"))

#======================================================================================================================


bodymass_table2 = merge(weighted_vocal_table, bodymasses, all.x = T)


flock_bodymass <- bodymass_table2 %>%
  group_by(flock_ID) %>%
  mutate(species_present = case_when(
    any(species_ID == "dicrurus_remifer") ~ "Dicrurus remifer",
    any(species_ID == "dicrurus_paradiseus") ~ "Dicrurus paradiseus",
    TRUE ~ "None"
  )) %>%
  filter(species_present != "None") %>%
  group_by(flock_ID, species_present, flocktype) %>%
  summarize(avg_bodymass = sum(count * bodymass) / sum(count), .groups = 'drop')


ggplot(flock_bodymass, aes(x = species_present, y = avg_bodymass, fill = species_present)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  labs(x = "", y = "Mean Body Mass (g)") +  
  stat_summary(fun = mean, geom = "point", shape = 24, size = 2, color = "black", fill = NA,
    position = position_dodge(0.9)) + theme_classic() + 
  theme(axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 12),  # Increase the size of x-axis text
        legend.position = "none") +
  scale_fill_manual(values = c("#4682B4", "#CD5B45"))

