setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")

library(ggplot2)
library(mgcv)
library(dplyr)
library(patchwork)

# Replace 'your_data.csv' with the path to your CSV file
data <- read.csv("weighted_vocal_table.csv")

# Ensure that the data is loaded correctly
head(data)

# Split the data by flock type
flock_types <- unique(data$flock_type)
models <- list()
plots <- list()

for (flock in flock_types) {
  # Filter data for the current flock type
  flock_data <- data %>% filter(flock_type == flock)
  
  # Fit a GLM
  glm_model <- glm(weighted_vocal_participation ~ half_ratio_index, data = flock_data)
  models[[flock]] <- glm_model
  
  # Create the plot
  p <- ggplot(flock_data, aes(x = half_ratio_index, y = weighted_vocal_participation)) +
    geom_point() +
    geom_smooth(method = "glm", method.args = list(family = "gaussian"), color = "blue", fill = "lightblue") +
    labs(x = "Half Ratio Index", y = "Weighted Vocal Participation", title = paste("Flock Type:", flock)) +
    theme_classic() +
    theme(axis.title = element_text(size = 20))
  
  plots[[flock]] <- p
}

# Combine the plots using patchwork
combined_plot <- Reduce(`+`, plots) + plot_layout(ncol = 1)
print(combined_plot)
