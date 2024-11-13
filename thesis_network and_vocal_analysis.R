
.libPaths("C:\\Users\\ronit\\OneDrive\\Documents\\R\\win-library\\4.1")

library(igraph)
library(igraphdata)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(glmmTMB)
install.packages("Rtools")


setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")

bird_data = read.csv("Namdapha_Flock_Data.csv", header = T)
flocktypes = read.csv("flocktypes.csv")
species_names <- colnames(bird_data)


S_flocks = bird_data %>% mutate(flocktype = flocktypes$flock_type) %>% filter(flocktype == "S") %>% 
  select(-flocktype) %>% select(where(~ sum(.) != 0)) %>% as.matrix()
LC_flocks = bird_data %>% mutate(flocktype = flocktypes$flock_type) %>% filter(flocktype == "LC") %>% 
  select(-flocktype) %>% select(where(~ sum(.) != 0)) %>% as.matrix()
LU_flocks = bird_data %>% mutate(flocktype = flocktypes$flock_type) %>% filter(flocktype == "LU") %>% 
  select(-flocktype) %>% select(where(~ sum(.) != 0)) %>% as.matrix()

all_flocks = ifelse(bird_data >= 1, 1, 0)
number_of_flocks = as.data.frame(colSums(all_flocks))
number_of_flocks = rownames_to_column(number_of_flocks, var = "species_ID")
colnames(number_of_flocks) <- c("species_ID", "number_of_flocks")

S_flocks = ifelse(S_flocks > 1, 1, S_flocks)
LC_flocks = ifelse(LC_flocks > 1, 1, LC_flocks)
LU_flocks = ifelse(LU_flocks > 1, 1, LU_flocks)


#Cutoff for number of flocks a species is present in
S_flocks2 = S_flocks[, colSums(S_flocks) > 7] #80 flocks 
LC_flocks2 = LC_flocks[, colSums(LC_flocks) > 4] # 44 flocks
LU_flocks2 = LU_flocks[, colSums(LU_flocks) > 2] # 24 flocks


# Calculate co-occurrence matrix
s_adj_matrix = t(S_flocks) %*% S_flocks
lc_adj_matrix = t(LC_flocks) %*% LC_flocks
lu_adj_matrix = t(LU_flocks) %*% LU_flocks

#For filtered networks
s_adj_matrix2 = t(S_flocks2) %*% S_flocks2
lc_adj_matrix2 = t(LC_flocks2) %*% LC_flocks2
lu_adj_matrix2 = t(LU_flocks2) %*% LU_flocks2

# Set diagonal to zero
diag(s_adj_matrix) = 0
diag(lc_adj_matrix) = 0
diag(lu_adj_matrix) = 0

#Filtered
diag(s_adj_matrix2) = 0
diag(lc_adj_matrix2) = 0
diag(lu_adj_matrix2) = 0


# Compute graph from adjacency matrix
s <- graph_from_adjacency_matrix(s_adj_matrix, weighted = TRUE, mode = "undirected")
lc <- graph_from_adjacency_matrix(lc_adj_matrix, weighted = TRUE, mode = "undirected")
lu <- graph_from_adjacency_matrix(lu_adj_matrix, weighted = TRUE, mode = "undirected")

#Filtered networks
s2 <- graph_from_adjacency_matrix(s_adj_matrix2, weighted = TRUE, mode = "undirected")
lc2 <- graph_from_adjacency_matrix(lc_adj_matrix2, weighted = TRUE, mode = "undirected")
lu2 <- graph_from_adjacency_matrix(lu_adj_matrix2, weighted = TRUE, mode = "undirected")

#-------------------------------------------------------------------------------------------------------------------------
# Calculate the degree centrality
s_degree <- degree(s, mode = "all", loops = FALSE, normalized = T)
lc_degree <- degree(lc, mode = "all", loops = FALSE, normalized = T)
lu_degree <- degree(lu, mode = "all", loops = FALSE, normalized = T)

#Weighted degree
s_strength = strength(s, mode = "all", loops = FALSE)
lc_strength = strength(lc)
lu_strength = strength(lu)
#-------------------------------------------------------------------------------------------------------------------------
# Set vertex attributes for S flocks
V(s2)$color <- "lightblue"                      # Vertex color
V(s2)$frame.color <- "black"
V(s2)$label <- V(s2)$name                        # Vertex label (if available)
V(s2)$label.cex <- 1.5                          # Vertex label font size

# Set edge attributes
E(s2)$color <- "orange"                           # Edge color
E(s2)$width <- E(s2)$weight / 5                      # Edge width based on weight

#-------------------------------------------------------------------------------------------------------------------------
# Set vertex attributes  for LC flocks
V(lc2)$color <- "lightblue"                      # Vertex color
V(lc2)$frame.color <- "black"
V(lc2)$label <- V(lc2)$name                        # Vertex label (if available)
V(lc2)$label.cex <- 1.5                          # Vertex label font size

# Set edge attributes
E(lc2)$color <- "orange"                           # Edge color
E(lc2)$width <- E(lc2)$weight / 2                      # Edge width based on weight

#-------------------------------------------------------------------------------------------------------------------------

# Set vertex attributes for LU flocks
V(lu2)$color <- "lightblue"                      # Vertex color
V(lu2)$frame.color <- "black"
V(lu2)$label <- V(lu2)$name                        # Vertex label (if available)
V(lu2)$label.cex <- 1.5                          # Vertex label font size

# Set edge attributes
E(lu2)$color <- "orange"                           # Edge color
E(lu2)$width <- E(lu2)$weight / 2                      # Edge width based on weight

#-------------------------------------------------------------------------------------------------------------------------

# Choose a layout for the graph
layout_s <- layout_nicely(s2)
layout_lc <- layout_nicely(lc2)
layout_lu <- layout_nicely(lu2)


#Set layout boundaries to max
par(mar=c(0,0,0,0))

# Plot the network with customized attributes
plot(s2, layout=layout_s, edge.curved=FALSE, vertex.label.color="black", vertex.label.cex=0.8)
plot(lc2, layout=layout_lc, edge.curved=FALSE, vertex.label.color="black", vertex.label.cex=0.8)
plot(lu2, layout=layout_lu, edge.curved=FALSE, vertex.label.color="black", vertex.label.cex=0.8)

#==========================================================================================================================


## Vocal participation part!!!!!!!

vocaldata = read.csv("vocal_participation.csv")
binary_data <- vocaldata %>% select(!(flock_ID:species_ID))


vocal_participation = rowMeans(binary_data == 1, na.rm = TRUE)

#Make a consolidated table with species IDs, flock number, counts, percentages & number of 1s and 0s
vocal_table = vocaldata %>% 
  select(flock_ID:species_ID) %>% mutate(vocal_participation, 
  "num_ones" = rowSums(binary_data == 1, na.rm = TRUE), "num_zeros" = rowSums(binary_data == 0, na.rm = TRUE))

#==========================================================================================================================


# Determine log curve for negative weights

x = vocal_table$count
y = vocal_table$vocal_participation

model <- lm(y ~ log(x))

# Get coefficients
a <- exp(coef(model)[1])
b <- coef(model)[2]

par(mar=c(4,4,4,4))


# Step 3: Add the fitted curve (linear regression)
#Fitting a logarithmic curve to the data:
x_seq <- seq(min(x), max(x), length.out = 100)
y_pred <- predict(model, newdata = data.frame(x = x_seq))

plot(x, y, ylab = "Unweighted Vocal Participation", xlab = "Number of Individuals")
lines(x_seq, y_pred, col = "red", lwd = 3)

#==========================================================================================================================

# Create a data frame
logcurve_df <- data.frame(count = vocal_table$count, 
                          vocal_participation = vocal_table$vocal_participation)

# Fit the model
logcurve_model <- lm(vocal_participation ~ log(count), data = logcurve_df)

# Get predictions for the fitted curve
logcurve_df$predicted <- predict(logcurve_model)

# Generate the ggplot
ggplot(logcurve_df, aes(x = count, y = vocal_participation)) +
  geom_point(alpha = 0.3, size = 3, color = "#0387ad") +
  geom_line(aes(y = predicted), color = "#d9534f", linewidth = 2) +
  labs(x = "Number of Individuals", y = "Vocal Participation") +
  theme_classic() + theme(axis.title = element_text(size = 20))

#Weighting vocal participation: t / 1.28logn + 1
#weighted_vocal_table = vocal_table %>% mutate(weight = log(count, base = 100), 
                                              #weighted_vocal_participation = vocal_participation / (weight+1))




#WEIGHTING EQUATION

weighted_vocal_table = vocal_table %>% 
  mutate(weight = y_pred[x], weighted_vocal_participation = vocal_participation / (weight+(1-b)))


#==========================================================================================================================

frequencies = as.data.frame(table(weighted_vocal_table$flock_ID))

#Add flocktypes column with number of reps based species count in each flock
weighted_vocal_table = weighted_vocal_table %>% mutate(flocktype = 
  unlist(lapply(seq_len(nrow(frequencies)), function(i) rep(flocktypes$flock_type[i], frequencies$Freq[i]))))

# Tables for each flock type - REMEMBER TO REMOVE BIRDS WITH LESS THAN 5-6 FLOCKS IN DATA!!!!!
s_table = weighted_vocal_table %>% filter(flocktype == "S") %>% 
  group_by(species_ID) %>% summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation),
                                     percent_flocks = (table(species_ID)/80))
lc_table = weighted_vocal_table %>% filter(flocktype == "LC") %>% 
  group_by(species_ID) %>% summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation),
                                     percent_flocks = (table(species_ID)/44))
lu_table = weighted_vocal_table %>% filter(flocktype == "LU") %>% 
  group_by(species_ID) %>% summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation),
                                     percent_flocks = (table(species_ID)/24))

#Degree centrality vs. Vocal Participation graph for S flocks--------------------------------------------------------------


s_degree = as.data.frame(s_degree)
s_degree$species_ID = rownames(s_degree)  
s_plot = merge(s_degree, s_table, by = "species_ID")
s_plot = s_plot %>% filter(species_ID %in% colnames(S_flocks2))
#s_plot$s_degree = s_plot$s_degree/max(s_plot$s_degree)

s_fit <- lm(mean_vp ~ s_degree, data = s_plot)
s_slope <- coef(s_fit)[2]
s_rsquared <- summary(s_fit)$r.squared
summary(s_fit) # From p-values, not statistically significant

#Code for regression plot
ggplot(s_plot, aes(x = s_degree, y = mean_vp)) + geom_point() + 
  geom_smooth(method = "lm", se = T, color = "blue", linewidth = 1.5) +  
  geom_text(aes(label = species_ID), size = 2.5, hjust = 1, vjust = 1.2) +  
  labs(x = "Degree Centrality", y = "Mean Vocal Participation") + theme_classic() +
  annotate("text", x = min(s_plot$s_degree), y = max(s_plot$mean_vp),
             label = paste("Slope =", round(s_slope, 3), "\nR-squared =", round(s_rsquared, 3)),
             hjust = 0, vjust = 1, size = 4) + theme(axis.title = element_text(size = 20))

hist(s_plot$mean_vp) # Normal

#Degree centrality vs. Vocal Participation graph for LC flocks-------------------------------------------------------------


lc_degree = as.data.frame(lc_degree)
lc_degree$species_ID = rownames(lc_degree)  
lc_plot = merge(lc_degree, lc_table, by = "species_ID")
lc_plot = lc_plot %>% filter(species_ID %in% colnames(LC_flocks2))
#lc_plot$lc_degree = lc_plot$lc_degree/max(lc_plot$lc_degree)

lc_fit <- lm(mean_vp ~ lc_degree, data = lc_plot)
lc_slope <- coef(lc_fit)[2]
lc_rsquared <- summary(lc_fit)$r.squared
summary(lc_fit) # From p-values, not statistically significant

#Code for regression plot
ggplot(lc_plot, aes(x = lc_degree, y = mean_vp)) + geom_point() + 
  geom_smooth(method = "lm", se = T, color = "blue", linewidth = 1.5) +  
  geom_text(aes(label = species_ID), size = 2.5, hjust = 1, vjust = 1.2) +  
  labs(x = "Degree Centrality", y = "Mean Vocal Participation") + theme_classic() +
  annotate("text", x = min(lc_plot$lc_degree), y = max(lc_plot$mean_vp),
           label = paste("Slope =", round(lc_slope, 3), "\nR-squared =", round(lc_rsquared, 3)),
           hjust = 0, vjust = 1, size = 4) + theme(axis.title = element_text(size = 20))
  
hist(lc_plot$mean_vp) # Normal

#Degree centrality vs. Vocal Participation graph for LU flocks-------------------------------------------------------------


lu_degree = as.data.frame(lu_degree)
lu_degree$species_ID = rownames(lu_degree)  
lu_plot = merge(lu_degree, lu_table, by = "species_ID")
lu_plot = lu_plot %>% filter(species_ID %in% colnames(LU_flocks2))
#lu_plot$lu_degree = lu_plot$lu_degree/max(lu_plot$lu_degree)


lu_fit <- lm(mean_vp ~ lu_degree, data = lu_plot)
lu_slope <- coef(lu_fit)[2]
lu_rsquared <- summary(lu_fit)$r.squared
summary(lu_fit) # Statistically significant

#Code for regression plot
ggplot(lu_plot, aes(x = lu_degree, y = mean_vp)) + geom_point() + 
  geom_smooth(method = "lm", se = T, color = "blue", linewidth = 1.5) +  
  geom_text(aes(label = species_ID), size = 2.5, hjust = 1, vjust = 1.2) +  
  labs(x = "Degree Centrality", y = "Mean Vocal Participation") + theme_classic() +
  annotate("text", x = min(lu_plot$lu_degree), y = max(lu_plot$mean_vp),
           label = paste("Slope =", round(lu_slope, 3), "\nR-squared =", round(lu_rsquared, 3)),
           hjust = 0, vjust = 1, size = 4) + theme(axis.title = element_text(size = 20))

hist(lu_plot$mean_vp) # Normal

#==========================================================================================================================

## Using Outside flock data to calculate flocking propensity and related metrics

#Flocking Propensity vs. Vocalness graphs
outsideflockdata = read.csv("OutsideFlockData.csv")
outsideflockcounts = outsideflockdata %>% group_by(species_ID) %>% summarise(outside_sum = sum(count))
insideflockcounts = vocaldata %>% group_by(species_ID) %>% summarise(inside_sum = sum(count))
all_table = weighted_vocal_table %>% group_by(species_ID) %>% 
  summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation))


propensity.table = merge(insideflockcounts, outsideflockcounts, by = "species_ID", all.x = T)
propensity.table[is.na(propensity.table)] <- 0
propensity.table = propensity.table %>% group_by(species_ID) %>% 
  mutate(total = sum(inside_sum, outside_sum)) %>% mutate(propensity = inside_sum/total)
propensity.table = merge(propensity.table, all_table, by = "species_ID")
propensity.table = propensity.table %>% mutate(percentage_abundance = (inside_sum/sum(inside_sum))*100)


#Remove species in 2 or less (>5%) flocks
bird_data2 = bird_data %>% as.matrix() %>% ifelse(. > 1, 1, .)
propensity.table2 = propensity.table %>% filter(species_ID %in% colnames(bird_data2[, colSums(bird_data2) > 2]))

#Vocal participation vs. Flocking Propensity Graph
ggplot(propensity.table2, aes(x = propensity, y = mean_vp)) + geom_point() + 
  geom_text(aes(label = species_ID), size = 2.5, hjust = 1, vjust = 1.2) +  
  labs(x = "Flocking Propensity", y = "Weighted Mean Vocal Participation") + theme_classic() +
  theme(axis.title = element_text(size = 20))

#Vocal participation vs. Abundance Graph for birds in >5% flocks -- FIT A GLM INSTEAD!!!!!!!!!!!!!
ggplot(propensity.table2, aes(x = inside_sum, y = mean_vp)) + geom_point() + 
  labs(x = "Abundance in Flocks", y = "Weighted Mean Vocal Participation") + theme_classic() +
  geom_smooth(method = "lm", se = T, color = "blue") + theme(axis.title = element_text(size = 20))


#Abundance in flocks vs. Overall Abundance Graph -- Not required
ggplot(propensity.table, aes(x = total, y = inside_sum)) + geom_point() + 
  geom_text(aes(label = species_ID), size = 2.5, hjust = 1, vjust = 1.2) +  
  labs(x = "Absolute Abundance", y = "Abundance in Flocks") + theme_classic() +
  geom_smooth(method = "lm", se = FALSE, color = "blue")

#==========================================================================================================================

library(moments)

# Histogram for Species-wise Flocking Propensity
ggplot(propensity.table, aes(x = propensity)) +
  geom_histogram(binwidth = .1, fill = "#0387ad", color = "white", boundary = 0) +
  labs(x = "Flocking Propensity", y = "Number of Species") +
  theme_classic() + theme(axis.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

#Histogram for Species-wise Abundances
ggplot(propensity.table, aes(x = percentage_abundance)) +
  geom_histogram(binwidth = 1, fill = "#0387ad", color = "white", boundary = 0) +
  labs(x = "Percentage (%) in Flocks", y = "Number of Species") +
  theme_classic() + theme(axis.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(0:8))

skewness(propensity.table$propensity)
skewness(propensity.table$percentage_abundance)

#==========================================================================================================================



propensity.table %>% filter(propensity >= .9)


#==========================================================================================================================


# Leadership vs. Vocal Participation Analysis

weighted_vocal_table2 = weighted_vocal_table %>% group_by(species_ID) %>% 
  summarise(leading_flocks = sum(leader), total_flocks = table(species_ID), mean_vp = mean(vocal_participation), 
            mean_weighted_vp = mean(weighted_vocal_participation)) %>%
  mutate(lead_proportion = leading_flocks/total_flocks)

weighted_vocal_table2 = weighted_vocal_table2 %>% filter(total_flocks > 2) %>% 
  mutate(leader_yn = ifelse(lead_proportion > 0, "Leaders", "Non-leaders"))

# Regression for all flocks -- Not required --- FIT A GLM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ggplot(weighted_vocal_table2, aes(x = lead_proportion, y = mean_vp)) + geom_point() + 
  labs(x = "Proportion of Flocks Lead", y = "Mean Vocal Participation") + theme_classic() +
  geom_smooth(method = "lm", se = T, color = "blue")

# Boxplot for Leaders vs. Non-Leaders
ggplot(weighted_vocal_table2, aes(x = factor(leader_yn), y = mean_vp, fill = factor(leader_yn))) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  geom_violin(trim = FALSE, alpha = 0.3, width = 0.4) +    # Decrease the width of the violin plot
  labs(x = "", y = "Mean Vocal Participation") + 
  theme_classic() + 
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = c("#CD5B45", "#0387ad")) +  # Purple and green colors
  scale_y_continuous(limits = c(-0.15, 1.08))

# T-test
t.test(mean_vp ~ leader_yn, data = weighted_vocal_table2)
?t.test

hist(weighted_vocal_table2$mean_vp)


#==========================================================================================================================


# SILENT PARTICIPANTS ANALYSIS 

weighted_vocal_table <- weighted_vocal_table %>%
  mutate(silent_participant = ifelse(weighted_vocal_participation > 0.0, 0, 1))

flock_total = weighted_vocal_table %>% group_by(flock_ID) %>% summarise(total_count = sum(count))
silent_total = weighted_vocal_table %>% filter(silent_participant == 1) %>% 
  group_by(flock_ID) %>% summarise(silent_count = sum(count))

silent_analysis_table = merge(flock_total, silent_total, by = "flock_ID")
silent_analysis_table = silent_analysis_table %>% 
  mutate(silent_ratio = silent_count/total_count, vocal_ratio = (total_count - silent_count)/total_count)

mean(silent_analysis_table$silent_ratio)
mean(silent_analysis_table$vocal_ratio)
#==========================================================================================================================

# Histogram for vocalness

vocal_distribution = weighted_vocal_table %>% group_by(species_ID) %>% 
  summarise(mean_vp = mean(weighted_vocal_participation))
hist(vocal_distribution$mean_vp)


# Histogram for Species-wise Vocalness Distribution

ggplot(vocal_distribution, aes(x = mean_vp)) +
  geom_histogram(binwidth = 0.1, fill = "#0387ad", color = "white", boundary = 0) +
  labs(x = "Mean Vocal Participation", y = "Number of Species") +
  theme_classic() + theme(axis.title = element_text(size = 18)) +
  xlim(0, 1) + ylim(0, 25) +
  scale_x_continuous(breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

skewness(vocal_distribution$mean_vp)
#==========================================================================================================================

# BEAUTIFYING GRAPHS

# Define y-axis limits
y_limits <- c(-.12, .85)  # Adjust these values as needed
x_limits <- c(0.4, 1)  # Adjust these values as needed

# Small-bodied Flocks plot
p1 <- ggplot(s_plot, aes(x = s_degree, y = mean_vp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", linewidth = 1.5, level = 0.95) +  
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  theme(axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ggtitle("a) Small-bodied Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + xlim(x_limits)


# Large-bodied Canopy Flocks plot
p2 <- ggplot(lc_plot, aes(x = lc_degree, y = mean_vp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", linewidth = 1.5, level = 0.95) +  
  labs(x = NULL, y = "Mean Vocal Participation") + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title = element_text(size = 15)) +
  ggtitle("b) Large-bodied Canopy Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + xlim(x_limits)

# Large-bodied Undergrowth Flocks plot
p3 <- ggplot(lu_plot, aes(x = lu_degree, y = mean_vp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", linewidth = 1.5, level = 0.95) +  
  labs(x = "Degree Centrality") + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank()) +
  ggtitle("c) Large-bodied Undergrowth Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + xlim(x_limits)

# Arrange the plots vertically
library(patchwork)

# Arrange the plots vertically
p1 + p2 + p3 + plot_layout(ncol = 1)


#==========================================================================================================================


# Vocal Participation vs. Average Group Size Graph -- GLM


propensity.table2 = merge(propensity.table2, number_of_flocks, by = "species_ID")
propensity.table2 = propensity.table2 %>% mutate(avg_group_size = inside_sum/number_of_flocks)

#Remove outlier
filtered_table <- propensity.table2 %>% filter(species_ID != "psittiparus_gularis")

# Harish code for the same plot, binomial instead of gamma
ggplot(data = filtered_table, aes(x = avg_group_size, y = mean_vp)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, 
              linewidth = 1.8, colour = "#02607A", fill = "#87CEEB") +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) + 
  theme_classic() +
  labs(x = "Mean Intraspecific Group Size", y = "Mean Vocal Participation") +
  theme(axis.title = element_text(size = 18)) + 
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = c(0:9))

#==========================================================================================================================

# Find slope, intercept, p values
# Fit the binomial regression model
binomial_groupsize <- glm(mean_vp ~ avg_group_size, data = filtered_table, family = binomial)

# Summarize the model
summary(binomial_groupsize)

# Extract coefficients (slope and intercept)
groupsize_coefficients <- summary(binomial_groupsize)$coefficients

# Intercept (linear/odds ratio)
exp(groupsize_coefficients[1, 1])

# Slope (linear/odds ratio)
exp(groupsize_coefficients[2, 1])

# P-values
groupsize_coefficients[, 4]

#==========================================================================================================================
#==========================================================================================================================
#==========================================================================================================================


#Beta regression instead of Binomial
filtered_table$mean_vp = filtered_table$mean_vp + 0.0000001
beta_groupsize <- glmmTMB(mean_vp ~ avg_group_size, data = filtered_table, family = beta_family(link = "logit"))
predict_groupsize <- predict(beta_groupsize, type = "response", se.fit = TRUE)

ggplot(filtered_table, aes(x = avg_group_size, y = mean_vp)) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  geom_line(aes(y = fitted(beta_groupsize)), color = '#02607A', linewidth = 1.5) +
  geom_ribbon(aes(ymin = predict_groupsize$fit - 1.96 * predict_groupsize$se.fit,
                  ymax = predict_groupsize$fit + 1.96 * predict_groupsize$se.fit),
                  alpha = 0.3, fill = "#0387AD") +
  labs(x = "Mean Intraspecific Group Size", y = "Mean Vocal Participation") +
  theme_classic() +
  theme(plot.title = element_text(size=10), axis.title = element_text(size = 20)) +
  scale_x_continuous(breaks = c(0:9))

# Coefficients
coefficients_groupsize <- summary(beta_groupsize)$coefficients$cond
intercept_groupsize_logit <- coefficients_groupsize["(Intercept)", "Estimate"]
slope_groupsize_logit <- coefficients_groupsize["avg_group_size", "Estimate"]

plogis(0.33888-0.06817)

# Transform the intercept and slope for groupsize
intercept_groupsize <- inv_logit(intercept_groupsize_logit)
slope_groupsize <- inv_logit(intercept_groupsize_logit + slope_groupsize_logit) - inv_logit(intercept_groupsize_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_groupsize, "\n")
cat("Slope on the linear scale:", slope_groupsize, "\n")


# Fit a GLM with a binomial distribution -- VP vs. Avg. Group Size
glm_fit_groupsize <- glm(mean_vp ~ avg_group_size, family = "binomial", data = filtered_table)

summary(glm_fit_groupsize) # Statistically significant
hist(resid(glm_fit_groupsize)) # Good fit
#==========================================================================================================================
#==========================================================================================================================
#==========================================================================================================================

# Fit a GLM with a binomial distribution -- VP vs. Leadership p/n graph


ggplot(data = weighted_vocal_table2, aes(x = lead_proportion, y = mean_vp)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), 
              se = TRUE, colour = "#02607A", fill = "#87CEEB", linewidth = 1.8) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  theme_classic() +
  labs(x = "Proportion of Flocks Led", y = "Mean Vocal Participation") +
  theme(axis.title = element_text(size = 18))

# Fit a GLM with a beta distribution -- VP vs. Leadership p/n graph

weighted_vocal_table2$mean_vp = weighted_vocal_table2$mean_vp + 0.0000001
beta_leadership = glmmTMB(mean_vp ~ lead_proportion, data = weighted_vocal_table2, family = beta_family(link = "logit"))
predict_leadership <- predict(beta_leadership, type = "response", se.fit = TRUE)


ggplot(weighted_vocal_table2, aes(x = lead_proportion, y = mean_vp)) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  geom_line(aes(y = fitted(beta_leadership)), color = '#02607A', linewidth = 1.5) +
  geom_ribbon(aes(ymin = predict_leadership$fit - 1.96 * predict_leadership$se.fit,
                  ymax = predict_leadership$fit + 1.96 * predict_leadership$se.fit),
              alpha = 0.3, fill = "#0387AD") +
  labs(x = "Proportion of Flocks Led", y = "Mean Vocal Participation") +
  theme_classic() +
  theme(plot.title = element_text(size=10), axis.title = element_text(size = 20))

#Coefficients
coefficients_lead <- summary(beta_leadership)$coefficients$cond
intercept_lead_logit <- coefficients_lead["(Intercept)", "Estimate"]
slope_lead_logit <- coefficients_lead["lead_proportion", "Estimate"]

plogis(6.0508-1.1508)


# Transform the intercept and slope for s flocks
intercept_lead <- inv_logit(intercept_lead_logit)
slope_lead <- inv_logit(intercept_lead_logit + slope_lead_logit) - inv_logit(intercept_lead_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_lead, "\n")
cat("Slope on the linear scale:", slope_lead, "\n")

#==========================================================================================================================
# Find slope, intercept, p values
# Fit the binomial regression model
binomial_leadership <- glm(mean_vp ~ lead_proportion, data = weighted_vocal_table2, family = binomial)

# Summarize the model
summary(binomial_leadership)

# Extract coefficients (slope and intercept)
leadership_coefficients <- summary(binomial_leadership)$coefficients

# Intercept (linear/odds ratio)
exp(leadership_coefficients[1, 1])

# Slope (linear/odds ratio)
exp(scale(leadership_coefficients[2, 1]))

# P-values
leadership_coefficients[, 4]

#==========================================================================================================================

lc_plot$mean_vp = lc_plot$mean_vp + 0.0000001
# Use beta GLM, logistic


library(betareg)
library(glmmTMB)

# Fit the beta regression model
beta_s <- glmmTMB(mean_vp ~ s_degree, data = s_plot, family = beta_family(link = "logit"))
beta_lc <- glmmTMB(mean_vp ~ lc_degree, data = lc_plot, family = beta_family(link = "logit"))
beta_lu <- glmmTMB(mean_vp ~ lu_degree, data = lu_plot, family = beta_family(link = "logit"))


summary(beta_s) # *
summary(beta_lc) #
summary(beta_lu) # *


hist(resid(beta_s))
hist(resid(beta_lc))
hist(resid(beta_lu))



# Get predicted values and standard errors
predict_s <- predict(beta_s, type = "response", se.fit = TRUE)
predict_lc <- predict(beta_lc, type = "response", se.fit = TRUE)
predict_lu <- predict(beta_lu, type = "response", se.fit = TRUE)


# Plot with ggplot2 --- FITTED ERRORBAR SHADING TO ALL!!!!!!!!!!!!!!!!!!!!!


aa = ggplot(s_plot, aes(x = s_degree, y = mean_vp)) +
  geom_point(size = 3, shape = 21, fill = "#FF7F50", colour = "#CD5B45", stroke = 1.2, alpha = .7) +
  geom_line(aes(y = fitted(beta_s)), colour = "#CD5B45", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_s$fit - 1.96 * predict_s$se.fit,
                  ymax = predict_s$fit + 1.96 * predict_s$se.fit),
              alpha = 0.3, fill = "#FF7F50") +
  labs(x = NULL, y = "Mean Vocal Participation") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("a) Small-bodied Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) +
  xlim(x_limits)

bb = ggplot(lc_plot, aes(x = lc_degree, y = mean_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#02a890", colour = "#027d68", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lc)), colour = "#027d68", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_lc$fit - 1.96 * predict_lc$se.fit,
                  ymax = predict_lc$fit + 1.96 * predict_lc$se.fit),
              alpha = 0.3, fill = "#02a890") +
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("b) Large-bodied Canopy Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + 
  xlim(x_limits)

cc = ggplot(lu_plot, aes(x = lu_degree, y = mean_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#87CEEB", colour = "#4682B4", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lu)), colour = "#4682B4", linewidth = 1.8) + 
  geom_ribbon(aes(ymin = predict_lu$fit - 1.96 * predict_lu$se.fit,
                  ymax = predict_lu$fit + 1.96 * predict_lu$se.fit),
              alpha = 0.3, fill = "#87CEEB") +
  labs(x = 'Degree Centrality', y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("c) Large-bodied Undergrowth Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + 
  xlim(x_limits)


 

aa + bb + cc + plot_layout(ncol = 1)
aa + bb + plot_layout(ncol = 2) + plot_annotation(caption = "Degree Centrality", 
                                                  theme = theme(plot.caption = element_text(hjust = 0.5, size=15)))


summary(beta_s)

#==========================================================================================================================

inv_logit <- function(x) exp(x) / (1 + exp(x))

# Find coefficients on the linear scale:
coefficients_s <- summary(beta_s)$coefficients$cond
intercept_s_logit <- coefficients_s["(Intercept)", "Estimate"]
slope_s_logit <- coefficients_s["s_degree", "Estimate"]



# Transform the intercept and slope for s flocks
intercept_s <- inv_logit(intercept_s_logit)
slope_s <- inv_logit(intercept_s_logit + slope_s_logit) - inv_logit(intercept_s_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_s, "\n")
cat("Slope on the linear scale:", slope_s, "\n")

#==========================================================================================================================

# LC flocks
coefficients_lc <- summary(beta_lc)$coefficients$cond
intercept_lc_logit <- coefficients_lc["(Intercept)", "Estimate"]
slope_lc_logit <- coefficients_lc["lc_degree", "Estimate"]

# Transform the intercept and slope for lc flocks
intercept_lc <- inv_logit(intercept_lc_logit)
slope_lc <- inv_logit(intercept_lc_logit + slope_lc_logit) - inv_logit(intercept_lc_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_lc, "\n")
cat("Slope on the linear scale:", slope_lc, "\n")

#==========================================================================================================================

# LU flocks
coefficients_lu <- summary(beta_lu)$coefficients$cond
intercept_lu_logit <- coefficients_lu["(Intercept)", "Estimate"]
slope_lu_logit <- coefficients_lu["lu_degree", "Estimate"]

# Transform the intercept and slope for lu flocks
intercept_lu <- inv_logit(intercept_lu_logit)
slope_lu <- inv_logit(intercept_lu_logit + slope_lu_logit) - inv_logit(intercept_lu_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_lu, "\n")
cat("Slope on the linear scale:", slope_lu, "\n")


#==========================================================================================================================

# Diversity across flock types boxplots


species_per_flock <- weighted_vocal_table %>%
  group_by(flock_ID, flocktype) %>%
  summarise(number_of_species = n_distinct(species_ID))

# Convert flocktype to a factor and specify the order of levels
species_per_flock$flocktype <- factor(species_per_flock$flocktype, 
                                      levels = c("S", "LC", "LU"),
                                      labels = c("Small-bodied\nFlocks", 
                                                 "Large-bodied\nCanopy Flocks", "Large-bodied\nUndergrowth Flocks"))

# Generate the box-violin plot
ggplot(species_per_flock, aes(x = flocktype, y = number_of_species, fill = flocktype)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 24,  # Inverted hollow triangle shape
    size = 3,    # Size of the triangle
    color = "black",
    fill = NA,   # No fill to make it hollow
    position = position_dodge(0.9)  # Align triangles with the boxes
  ) +
  labs(x = "", y = "Number of Species") +  
  theme_classic() + 
  theme(axis.title = element_text(size = 18), 
        axis.text.x = element_text(size = 12),  # Increase the size of x-axis text
        legend.position = "none") +
  scale_fill_manual(values = c("Small-bodied\nFlocks" = "#CD5B45", 
                               "Large-bodied\nCanopy Flocks" = "#027d68", 
                               "Large-bodied\nUndergrowth Flocks" = "#4682B4"))

species_per_flock %>% group_by(flocktype) %>% summarize(mean(number_of_species))

#==========================================================================================================================

# Body mass graphs

bodymasses = read.csv("body_mass_table.csv", header = T)

bodymass_table = merge(weighted_vocal_table, bodymasses, all.x = T)

bodymass_table = bodymass_table %>% group_by(flock_ID, flocktype) %>% summarise(mean_bodymass = mean(bodymass))
bodymass_table$flocktype <- factor(bodymass_table$flocktype, 
                                      levels = c("S", "LC", "LU"),
                                      labels = c("Small-bodied\nFlocks", 
                                                 "Large-bodied\nCanopy Flocks", "Large-bodied\nUndergrowth Flocks"))

ggplot(bodymass_table, aes(x = flocktype, y = mean_bodymass, fill = flocktype)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  geom_violin(trim = FALSE, alpha = 0.3, width = 0.6) +    # Decrease the width of the violin plot
  labs(x = "", y = "Mean Body Mass (g)") +  
  theme_classic() + 
  theme(axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 12),  # Increase the size of x-axis text
        legend.position = "none") +
  scale_fill_manual(values = c("Small-bodied\nFlocks" = "#CD5B45", 
                               "Large-bodied\nCanopy Flocks" = "#027d68", 
                               "Large-bodied\nUndergrowth Flocks" = "#4682B4"))

bodymass_table %>% group_by(flocktype) %>% summarise(mean_bm = mean(mean_bodymass))

#==========================================================================================================================


# Fit the beta regression model
beta_sp <- glmmTMB(mean_vp ~ percent_flocks, data = s_plot, family = beta_family(link = "logit"))
beta_lcp <- glmmTMB(mean_vp ~ percent_flocks, data = lc_plot, family = beta_family(link = "logit"))
beta_lup <- glmmTMB(mean_vp ~ percent_flocks, data = lu_plot, family = beta_family(link = "logit"))


summary(beta_sp) # *
summary(beta_lcp) #
summary(beta_lup) # *

plogis(1.813)

hist(resid(beta_sp))
hist(resid(beta_lcp))
hist(resid(beta_lup))



# Get predicted values and standard errors
predict_sp <- predict(beta_sp, type = "response", se.fit = TRUE)
predict_lcp <- predict(beta_lcp, type = "response", se.fit = TRUE)
predict_lup <- predict(beta_lup, type = "response", se.fit = TRUE)


# Plot with ggplot2 --- FITTED ERRORBAR SHADING TO ALL!!!!!!!!!!!!!!!!!!!!!

y_limits2 <- c(-.12, 1.05)


ap = ggplot(s_plot, aes(x = percent_flocks, y = mean_vp)) +
  geom_point(size = 3, shape = 21, fill = "#FF7F50", colour = "#CD5B45", stroke = 1.2, alpha = .7) +
  geom_line(aes(y = fitted(beta_sp)), colour = "#CD5B45", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_sp$fit - 1.96 * predict_sp$se.fit,
                  ymax = predict_sp$fit + 1.96 * predict_sp$se.fit),
              alpha = 0.3, fill = "#FF7F50") +
  labs(x = NULL, y = "Mean Vocal Participation") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("a) Small-bodied Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits2) +
  scale_x_continuous(limits = c(0,1))

# Plot with bb colors
bp = ggplot(lc_plot, aes(x = percent_flocks, y = mean_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#02a890", colour = "#027d68", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lcp)), colour = "#027d68", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_lcp$fit - 1.96 * predict_lcp$se.fit,
                  ymax = predict_lcp$fit + 1.96 * predict_lcp$se.fit),
              alpha = 0.3, fill = "#02a890") +
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("b) Large-bodied Canopy Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits2) +
  scale_x_continuous(limits = c(0,1))

# Plot with cc colors
cp = ggplot(lu_plot, aes(x = percent_flocks, y = mean_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#87CEEB", colour = "#4682B4", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lup)), colour = "#4682B4", linewidth = 1.8) + 
  geom_ribbon(aes(ymin = predict_lup$fit - 1.96 * predict_lup$se.fit,
                  ymax = predict_lup$fit + 1.96 * predict_lup$se.fit),
              alpha = 0.3, fill = "#87CEEB") +
  labs(x = 'Proportion of Flocks', y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("c) Large-bodied Undergrowth Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits2) +
  scale_x_continuous(limits = c(0,1))


ap + bp + cp + plot_layout(ncol = 1)
ap + bp + plot_layout(ncol = 2) + plot_annotation(caption = "Proportion in Flocks", 
                                                  theme = theme(plot.caption = element_text(hjust = 0.5, size=15)))

summary(beta_lcp)
plogis(3.5779+1.9739)
#==========================================================================================================================


inv_logit <- function(x) exp(x) / (1 + exp(x))

# Find coefficients on the linear scale:
coefficients_sp <- summary(beta_sp)$coefficients$cond
intercept_sp_logit <- coefficients_sp["(Intercept)", "Estimate"]
slope_sp_logit <- coefficients_sp["percent_flocks", "Estimate"]

# Transform the intercept and slope for sp flocks
intercept_sp <- inv_logit(intercept_sp_logit)
slope_sp <- inv_logit(intercept_sp_logit + slope_sp_logit) - inv_logit(intercept_sp_logit)

# Print the results
summary(beta_sp)
cat("Intercept on the linear scale:", intercept_sp, "\n")
cat("Slope on the linear scale:", slope_sp, "\n")

#==========================================================================================================================

# LCP flocks
coefficients_lcp <- summary(beta_lcp)$coefficients$cond
intercept_lcp_logit <- coefficients_lcp["(Intercept)", "Estimate"]
slope_lcp_logit <- coefficients_lcp["percent_flocks", "Estimate"]

# Transform the intercept and slope for lcp flocks
intercept_lcp <- inv_logit(intercept_lcp_logit)
slope_lcp <- inv_logit(intercept_lcp_logit + slope_lcp_logit) - inv_logit(intercept_lcp_logit)

# Print the results
summary(beta_lcp)
cat("Intercept on the linear scale:", intercept_lcp, "\n")
cat("Slope on the linear scale:", slope_lcp, "\n")

#==========================================================================================================================

# LUP flocks
coefficients_lup <- summary(beta_lup)$coefficients$cond
intercept_lup_logit <- coefficients_lup["(Intercept)", "Estimate"]
slope_lup_logit <- coefficients_lup["percent_flocks", "Estimate"]

# Transform the intercept and slope for lup flocks
intercept_lup <- inv_logit(intercept_lup_logit)
slope_lup <- inv_logit(intercept_lup_logit + slope_lup_logit) - inv_logit(intercept_lup_logit)

# Print the results
summary(beta_lup)
cat("Intercept on the linear scale:", intercept_lup, "\n")
cat("Slope on the linear scale:", slope_lup, "\n")

#==========================================================================================================================

# Tables for supplementary figures

write.csv(s_plot, file= "s_plot_table.csv")
write.csv(lc_plot, file= "lc_plot_table.csv")
write.csv(lu_plot, file= "lu_plot_table.csv")

#==========================================================================================================================

# Degree centrality vs. proportion of flocks

# ggplot(s_plot, aes(x = percent_flocks, y = s_degree)) +
#   geom_point(size = 3, shape = 21, fill = "#FF7F50", colour = "#CD5B45", stroke = 1.2, alpha = .7) +
#   geom_smooth(method = "lm", colour = "#CD5B45", fill = "#FF7F50", alpha = 0.3, se = TRUE, size = 1.8) +
#   labs(x = NULL, y = NULL) +
#   theme_classic() +
#   theme(plot.title = element_text(size=10)) +
#   ggtitle("a) Small-bodied Flocks") +
#   scale_y_continuous(limits = c(0.4,1)) +
#   scale_x_continuous(limits = c(0, 1))
# 
# 
# ggplot(lc_plot, aes(x = percent_flocks, y = lc_degree)) + 
#   geom_point(size = 3, shape = 21, fill = "#02a890", colour = "#027d68", stroke = 1.2, alpha = .7) +
#   geom_smooth(method = "lm", colour = "#027d68", fill = "#02a890", alpha = 0.3, se = TRUE, size = 1.8) +
#   labs(x = NULL, y = "Degree Centrality") + 
#   theme_classic() +
#   theme(axis.title = element_text(size = 15),
#         axis.title.x = element_blank(), plot.title = element_text(size=10)) +
#   ggtitle("b) Large-bodied Canopy Flocks") +
#   scale_y_continuous(limits = c(0.4,1)) +
#   scale_x_continuous(limits = c(0,1))
# 
# 
# ggplot(lu_plot, aes(x = percent_flocks, y = lu_degree)) + 
#   geom_point(size = 3, shape = 21, fill = "#87CEEB", colour = "#4682B4", stroke = 1.2, alpha = .7) + 
#   geom_smooth(method = "lm", colour = "#4682B4", fill = "#87CEEB", alpha = 0.3, se = TRUE, size = 1.8) +
#   labs(x = 'Proportion of Flocks', y = NULL) + 
#   theme_classic() +
#   theme(axis.title = element_text(size = 15),
#         axis.title.y = element_blank(), plot.title = element_text(size=10)) +
#   ggtitle("c) Large-bodied Undergrowth Flocks") +
#   scale_y_continuous(limits = c(0.4,1)) +
#   scale_x_continuous(limits = c(0,1))


#==========================================================================================================================
# 
# library(car)
# s_plot <- s_plot %>%
#   left_join(select(weighted_vocal_table2, species_ID, lead_proportion), by = "species_ID")
# 
# s_plot2 = s_plot %>% select(s_degree, percent_flocks, lead_proportion)
# s_plot2 = scale(s_plot2)
# 
# smodel = lm(mean_vp ~ s_degree + percent_flocks + lead_proportion, data = s_plot)
# 


library(car)
vif(smodel)
# s_pca = prcomp(s_plot2, center = TRUE, scale. = TRUE)
# summary(s_pca)

# glmmTMB(mean_vp ~ s_pca$x[,1] + s_pca$x[,2], data = s_plot, family = beta_family(link = "logit"))
# glmmTMB(mean_vp ~ s_pca$x[,1] + s_pca$x[,3], data = s_plot, family = beta_family(link = "logit"))
# glmmTMB(mean_vp ~ s_pca$x[,2] + s_pca$x[,3], data = s_plot, family = beta_family(link = "logit"))
# 
# glmmTMB(mean_vp ~ s_pca$x[,1], data = s_plot, family = beta_family(link = "logit"))
# glmmTMB(mean_vp ~ s_pca$x[,2], data = s_plot, family = beta_family(link = "logit"))
# glmmTMB(mean_vp ~ s_pca$x[,3], data = s_plot, family = beta_family(link = "logit"))
# 
# glmmTMB(mean_vp ~ s_pca$x[,1] * s_pca$x[,2] * s_pca$x[,3], data = s_plot, family = beta_family(link = "logit"))


# S flocks
s_plot <- s_plot %>%
  left_join(select(weighted_vocal_table2, species_ID, lead_proportion), by = "species_ID")

s_plot <- s_plot %>%
  left_join(select(filtered_table, species_ID, avg_group_size), by = "species_ID")

glmmTMB(mean_vp ~ s_degree + avg_group_size + lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--22.24856
glmmTMB(mean_vp ~ s_degree * avg_group_size * lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--21.59023
glmmTMB(mean_vp ~ s_degree + avg_group_size, data = s_plot, family = beta_family(link = "logit")) #--20.72585
glmmTMB(mean_vp ~ avg_group_size + lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--15.18746
glmmTMB(mean_vp ~ s_degree + lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--23.48310

glmmTMB(mean_vp ~ s_degree * avg_group_size, data = s_plot, family = beta_family(link = "logit")) #--18.77546
glmmTMB(mean_vp ~ avg_group_size * lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--16.40909
glmmTMB(mean_vp ~ s_degree * lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--21.49121

glmmTMB(mean_vp ~ s_degree, data = s_plot, family = beta_family(link = "logit")) #--11.5654
glmmTMB(mean_vp ~ avg_group_size, data = s_plot, family = beta_family(link = "logit")) #--11.66446
glmmTMB(mean_vp ~ lead_proportion, data = s_plot, family = beta_family(link = "logit")) #--17.13762


# LC flocks
lc_plot <- lc_plot %>%
  left_join(select(weighted_vocal_table2, species_ID, lead_proportion), by = "species_ID")

lc_plot <- lc_plot %>%
  left_join(select(filtered_table, species_ID, avg_group_size), by = "species_ID")

glmmTMB(mean_vp ~ lc_degree + avg_group_size + lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--17.97094
glmmTMB(mean_vp ~ lc_degree * avg_group_size * lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--16.78376
glmmTMB(mean_vp ~ lc_degree + avg_group_size, data = lc_plot, family = beta_family(link = "logit")) #--19.77086
glmmTMB(mean_vp ~ avg_group_size + lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--12.86252
glmmTMB(mean_vp ~ lc_degree + lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--15.83052

glmmTMB(mean_vp ~ lc_degree * avg_group_size, data = lc_plot, family = beta_family(link = "logit")) #--18.78845
glmmTMB(mean_vp ~ avg_group_size * lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--17.10235
glmmTMB(mean_vp ~ lc_degree * lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--14.23428

glmmTMB(mean_vp ~ lc_degree, data = lc_plot, family = beta_family(link = "logit")) #--8.10534
glmmTMB(mean_vp ~ avg_group_size, data = lc_plot, family = beta_family(link = "logit")) #--13.85661
glmmTMB(mean_vp ~ lead_proportion, data = lc_plot, family = beta_family(link = "logit")) #--13.90715


# LU Flocks
lu_plot <- lu_plot %>%
  left_join(select(weighted_vocal_table2, species_ID, lead_proportion), by = "species_ID")

glmmTMB(mean_vp ~ lu_degree + percent_flocks + lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-5.91411
glmmTMB(mean_vp ~ lu_degree * percent_flocks * lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-14.36937
glmmTMB(mean_vp ~ lu_degree + percent_flocks, data = lu_plot, family = beta_family(link = "logit")) #-5.81655
glmmTMB(mean_vp ~ percent_flocks + lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-7.76460
glmmTMB(mean_vp ~ lu_degree + lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-7.65616

glmmTMB(mean_vp ~ lu_degree * percent_flocks, data = lu_plot, family = beta_family(link = "logit")) #-20.72776
glmmTMB(mean_vp ~ percent_flocks * lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-7.74139
glmmTMB(mean_vp ~ lu_degree * lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-7.86385

glmmTMB(mean_vp ~ lu_degree, data = lu_plot, family = beta_family(link = "logit")) #-3.34297
glmmTMB(mean_vp ~ percent_flocks, data = lu_plot, family = beta_family(link = "logit")) #-7.09244
glmmTMB(mean_vp ~ lead_proportion, data = lu_plot, family = beta_family(link = "logit")) #-9.45848


library(mgcv)
#==========================================================================================================================

# Model Selection

library(AICcmodavg)
library(glmulti)

s_full_model = glmmTMB(mean_vp ~ avg_group_size + lead_proportion + s_degree + percent_flocks + 
                         avg_group_size:lead_proportion + avg_group_size:s_degree + 
                         avg_group_size:percent_flocks + lead_proportion:s_degree + 
                         lead_proportion:percent_flocks + s_degree:percent_flocks, 
                       family = beta_family(link = "logit"), data = s_plot)

s_model_selection <- aictab(s_full_model)

# Get AIC-based model selection table
s_model_selection_table <- s_model_selection$tab

s_selection_result <- glmulti(y = "mean_vp", xr = c("avg_group_size", "lead_proportion", "s_degree", "percent_flocks"), 
                              data = s_plot, intercept = F, level = 2, maxsize = 2, method = "h", crit = "aic", plotty = F, 
                              family = "binomial")
                              
get.models(s_selection_result, subset = "bestaic")



weighted_vocal_table %>%
  filter(leader == 1) %>%
  distinct(flock_ID) %>%
  nrow()


paper_table = propensity.table2 <- propensity.table2 %>%
  left_join(select(weighted_vocal_table2, species_ID, lead_proportion), by = "species_ID")


write.csv(paper_table, "paper_table.csv")
