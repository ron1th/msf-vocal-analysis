#Initial Data processing using the flock matrix (Namdapha_Flock_Data.csv) and flock types (flocktypes.csv)

bird_data = read.csv("./data/Namdapha_Flock_Data.csv", header = T)
flocktypes = read.csv("./data/flocktypes.csv")
species_names <- colnames(bird_data)

select <- dplyr::select

#Split matrices by flocktype
S_flocks = bird_data %>% mutate(flocktype = flocktypes$flock_type) %>% filter(flocktype == "S") %>% 
  select(-flocktype) %>% select(where(~ sum(.) != 0)) %>% as.matrix()
LC_flocks = bird_data %>% mutate(flocktype = flocktypes$flock_type) %>% filter(flocktype == "LC") %>% 
  select(-flocktype) %>% select(where(~ sum(.) != 0)) %>% as.matrix()
LU_flocks = bird_data %>% mutate(flocktype = flocktypes$flock_type) %>% filter(flocktype == "LU") %>% 
  select(-flocktype) %>% select(where(~ sum(.) != 0)) %>% as.matrix()

# Convert counts to presence/absence (1/0)
all_flocks = ifelse(bird_data >= 1, 1, 0)

S_flocks = ifelse(S_flocks > 1, 1, S_flocks)
LC_flocks = ifelse(LC_flocks > 1, 1, LC_flocks)
LU_flocks = ifelse(LU_flocks > 1, 1, LU_flocks)

# No. of Flocks each species was found in
number_of_flocks <- enframe(colSums(all_flocks), name = "species_ID", value = "number_of_flocks")

#Cutoffs for number of flocks a species is present in
S_flocks2 = S_flocks[, colSums(S_flocks) > 7] #80 flocks 
LC_flocks2 = LC_flocks[, colSums(LC_flocks) > 4] # 43 flocks
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

# Calculate the degree centrality & Weighted Degree ============================================

s_degree <- degree(s, mode = "all", loops = FALSE, normalized = T)
lc_degree <- degree(lc, mode = "all", loops = FALSE, normalized = T)
lu_degree <- degree(lu, mode = "all", loops = FALSE, normalized = T)

#Weighted degree
s_strength = strength(s, mode = "all", loops = FALSE)
lc_strength = strength(lc, mode = "all", loops = FALSE)
lu_strength = strength(lu, mode = "all", loops = FALSE)

# Read Vocal Data and create vocal table =====================================================================

vocaldata = read.csv("./data/vocal_participation.csv")
binary_data <- vocaldata %>% select(!(flock_ID:species_ID))


vocal_participation = rowMeans(binary_data == 1, na.rm = TRUE)

#Make a consolidated table with species IDs, flock number, counts, percentages & number of 1s and 0s
vocal_table = vocaldata %>% 
  select(flock_ID:species_ID) %>% mutate(vocal_participation, 
                                         "num_ones" = rowSums(binary_data == 1, na.rm = TRUE), "num_zeros" = rowSums(binary_data == 0, na.rm = TRUE))

#Weighting the vocal table

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

weighted_vocal_table = vocal_table %>% 
  mutate(weight = y_pred[x], weighted_vocal_participation = vocal_participation / (weight+(1-b)))
# Tables for each flock type ==========================================================================================================================

frequencies = as.data.frame(table(weighted_vocal_table$flock_ID))

#Add flocktypes column with number of reps based species count in each flock
weighted_vocal_table = weighted_vocal_table %>% mutate(flocktype = 
                                                         unlist(lapply(seq_len(nrow(frequencies)), function(i) rep(flocktypes$flock_type[i], frequencies$Freq[i]))))

# Tables for each flock type
s_table = weighted_vocal_table %>% filter(flocktype == "S") %>% 
  group_by(species_ID) %>% summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation),
                                     percent_flocks = (table(species_ID)/80))
lc_table = weighted_vocal_table %>% filter(flocktype == "LC") %>% 
  group_by(species_ID) %>% summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation),
                                     percent_flocks = (table(species_ID)/44))
lu_table = weighted_vocal_table %>% filter(flocktype == "LU") %>% 
  group_by(species_ID) %>% summarise(mean_vp = mean(vocal_participation), mean_weighted_vp = mean(weighted_vocal_participation),
                                     percent_flocks = (table(species_ID)/24))

# Degree centrality vs. Vocal Participation graph for S flocks ==========================================================================================

s_degree = as.data.frame(s_degree)
s_degree$species_ID = rownames(s_degree)  
s_plot = merge(s_degree, s_table, by = "species_ID")
s_plot = s_plot %>% filter(species_ID %in% colnames(S_flocks2))
#s_plot$s_degree = s_plot$s_degree/max(s_plot$s_degree)

s_fit <- lm(mean_vp ~ s_degree, data = s_plot)
s_slope <- coef(s_fit)[2]
s_rsquared <- summary(s_fit)$r.squared

#Degree centrality vs. Vocal Participation graph for LC flocks =============================================================================================

lc_degree = as.data.frame(lc_degree)
lc_degree$species_ID = rownames(lc_degree)  
lc_plot = merge(lc_degree, lc_table, by = "species_ID")
lc_plot = lc_plot %>% filter(species_ID %in% colnames(LC_flocks2))
#lc_plot$lc_degree = lc_plot$lc_degree/max(lc_plot$lc_degree)

lc_fit <- lm(mean_vp ~ lc_degree, data = lc_plot)
lc_slope <- coef(lc_fit)[2]
lc_rsquared <- summary(lc_fit)$r.squared

#Degree centrality vs. Vocal Participation graph for LU flocks =====================================================================================================================

lu_degree = as.data.frame(lu_degree)
lu_degree$species_ID = rownames(lu_degree)  
lu_plot = merge(lu_degree, lu_table, by = "species_ID")
lu_plot = lu_plot %>% filter(species_ID %in% colnames(LU_flocks2))
#lu_plot$lu_degree = lu_plot$lu_degree/max(lu_plot$lu_degree)


lu_fit <- lm(mean_vp ~ lu_degree, data = lu_plot)
lu_slope <- coef(lu_fit)[2]
lu_rsquared <- summary(lu_fit)$r.squared

# Using Outside flock data to calculate flocking propensity and related metrics =====================================================================================================================

#Flocking Propensity vs. Vocal activity graphs
outsideflockdata = read.csv("./data/OutsideFlockData.csv")
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

# Split by Leaders and Non-leaders ============================================================================================================================================

weighted_vocal_table2 = weighted_vocal_table %>% group_by(species_ID) %>% 
  summarise(leading_flocks = sum(leader), total_flocks = table(species_ID), mean_vp = mean(vocal_participation), 
            mean_weighted_vp = mean(weighted_vocal_participation)) %>%
  mutate(lead_proportion = leading_flocks/total_flocks)

weighted_vocal_table2 = weighted_vocal_table2 %>% filter(total_flocks > 2) %>% 
  mutate(leader_yn = ifelse(lead_proportion > 0, "Leaders", "Non-leaders"))

# Vocal Distribution Table ============================================================================================================================================

vocal_distribution = weighted_vocal_table %>% group_by(species_ID) %>% 
  summarise(mean_vp = mean(weighted_vocal_participation))

# Modify & Filter Propensity Table  ==========================================================================================================================

propensity.table2 = merge(propensity.table2, number_of_flocks, by = "species_ID")
propensity.table2 = propensity.table2 %>% mutate(avg_group_size = inside_sum/number_of_flocks)

#Remove outlier
filtered_table <- propensity.table2 %>% filter(species_ID != "psittiparus_gularis")


#For Groupsize vs. VP Beta regression ==========================================================================================================================

filtered_table$mean_vp = filtered_table$mean_vp + 0.0000001
beta_groupsize <- glmmTMB(mean_vp ~ avg_group_size, data = filtered_table, family = beta_family(link = "logit"))
predict_groupsize <- predict(beta_groupsize, type = "response", se.fit = TRUE)

#For Leadership Proportion vs. VP Beta regression ==========================================================================================================================

weighted_vocal_table2$mean_vp = weighted_vocal_table2$mean_vp + 0.0000001
beta_leadership = glmmTMB(mean_vp ~ lead_proportion, data = weighted_vocal_table2, family = beta_family(link = "logit"))
predict_leadership <- predict(beta_leadership, type = "response", se.fit = TRUE)

# For Centrality vs. VP Beta GLMs ==========================================================================================================================

lc_plot$mean_vp = lc_plot$mean_vp + 0.0000001 #jitter

library(betareg)
library(glmmTMB)

# Fit the beta regression model
beta_s <- glmmTMB(mean_vp ~ s_degree, data = s_plot, family = beta_family(link = "logit"))
beta_lc <- glmmTMB(mean_vp ~ lc_degree, data = lc_plot, family = beta_family(link = "logit"))
beta_lu <- glmmTMB(mean_vp ~ lu_degree, data = lu_plot, family = beta_family(link = "logit"))

# Get predicted values and standard errors
predict_s <- predict(beta_s, type = "response", se.fit = TRUE)
predict_lc <- predict(beta_lc, type = "response", se.fit = TRUE)
predict_lu <- predict(beta_lu, type = "response", se.fit = TRUE)

# For Body mass box-violin plots ==========================================================================================================================

bodymasses = read.csv("./data/body_mass_table.csv", header = T)

bodymass_table = merge(weighted_vocal_table, bodymasses, all.x = T)

bodymass_table = bodymass_table %>% group_by(flock_ID, flocktype) %>% summarise(mean_bodymass = mean(bodymass))
bodymass_table$flocktype <- factor(bodymass_table$flocktype, 
                                   levels = c("S", "LC", "LU"),
                                   labels = c("Small-bodied\nFlocks", 
                                              "Large-bodied\nCanopy Flocks", "Large-bodied\nUndergrowth Flocks"))

# For Richness box-violin plots ==========================================================================================================================

species_per_flock <- weighted_vocal_table %>%
  group_by(flock_ID, flocktype) %>%
  summarise(number_of_species = n_distinct(species_ID))

# Convert flocktype to a factor and specify the order of levels
species_per_flock$flocktype <- factor(species_per_flock$flocktype, 
                                      levels = c("S", "LC", "LU"),
                                      labels = c("Small-bodied\nFlocks", 
                                                 "Large-bodied\nCanopy Flocks", "Large-bodied\nUndergrowth Flocks"))

# For Proportion in Flocks vs. VP GLM ==========================================================================================================================

# Fit the beta regression model
beta_sp <- glmmTMB(mean_vp ~ percent_flocks, data = s_plot, family = beta_family(link = "logit"))
beta_lcp <- glmmTMB(mean_vp ~ percent_flocks, data = lc_plot, family = beta_family(link = "logit"))
beta_lup <- glmmTMB(mean_vp ~ percent_flocks, data = lu_plot, family = beta_family(link = "logit"))

# Get predicted values and standard errors
predict_sp <- predict(beta_sp, type = "response", se.fit = TRUE)
predict_lcp <- predict(beta_lcp, type = "response", se.fit = TRUE)
predict_lup <- predict(beta_lup, type = "response", se.fit = TRUE)

# For Flock NMDS ==========================================================================================================================

bird_data <- read.csv("./data/Namdapha_Flock_Data.csv", header = T)

nmds = metaMDS(bird_data, distance = "bray", k = 2, maxit = 999)
nmds_dat = as_tibble(nmds$points)
nmds_dat = nmds_dat %>% mutate(FlockType = as.factor(flocktypes$flock_type))
