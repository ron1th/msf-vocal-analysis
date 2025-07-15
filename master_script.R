library(igraph)
library(igraphdata)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(glmmTMB)
library(moments)
library(PERMANOVA)
library(iNEXT)
library(vegan)
library(patchwork)
library(DHARMa)
library(performance)

source("./scripts/data_processing.R")


# Histogram for Species-wise Flocking Propensity======================================================================

flocking_propensity_hist = 
  ggplot(propensity.table, aes(x = propensity)) +
  geom_histogram(binwidth = .1, fill = "#0387ad", color = "white", boundary = 0) +
  labs(x = "Flocking Propensity", y = "Number of Species") +
  theme_classic() + theme(axis.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

ggsave(filename = "./outputs/plots/flocking_propensity_hist.jpg", plot = flocking_propensity_hist, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")


#Histogram for Species-wise Abundances================================================================================

abundance_hist = 
  ggplot(propensity.table, aes(x = percentage_abundance)) +
  geom_histogram(binwidth = 1, fill = "#0387ad", color = "white", boundary = 0) +
  labs(x = "Percentage (%) in Flocks", y = "Number of Species") +
  theme_classic() + theme(axis.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(0:8))

ggsave(filename = "./outputs/plots/abundance_hist.jpg", plot = abundance_hist, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")

# Histogram for Species-wise Vocal Distribution=======================================================================

vocal_dist_hist = 
  ggplot(vocal_distribution, aes(x = mean_vp)) +
  geom_histogram(binwidth = 0.1, fill = "#0387ad", color = "white", boundary = 0) +
  labs(x = "Mean Vocal Participation", y = "Number of Species") +
  theme_classic() + theme(axis.title = element_text(size = 18)) +
  xlim(0, 1) + ylim(0, 25) +
  scale_x_continuous(breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

ggsave(filename = "./outputs/plots/vocal_dist_hist.jpg", plot = vocal_dist_hist, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")


# Boxplot for Leaders vs. Non-Leaders=================================================================================

leaders_nonleaders_boxplot = 
  ggplot(weighted_vocal_table2, aes(x = factor(leader_yn), y = mean_vp, fill = factor(leader_yn))) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  labs(x = "", y = "Mean Vocal Participation") + 
  theme_classic() + 
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c("#CD5B45", "#0387ad")) +  # Purple and green colors
  scale_y_continuous(limits = c(0, 1.00), breaks = c(0, 0.25,0.50, 0.75, 1.00))

ggsave(filename = "./outputs/plots/leaders_nonleaders_boxplot.jpg", plot = leaders_nonleaders_boxplot, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")

# Richness across flock types boxplots ===================================================================================

# Generate the Diversity across flock types box-violin plot

richness_boxplot = 
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

ggsave(filename = "./outputs/plots/richness_boxplot.jpg", plot = richness_boxplot, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")

# Body mass box-violin plots ===================================================================================================================

bodymass_boxplot = 
  ggplot(bodymass_table, aes(x = flocktype, y = mean_bodymass, fill = flocktype)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  stat_summary(fun = mean, geom = "point", shape = 24, size = 2, color = "black", fill = NA, position = position_dodge(0.9)) +
  labs(x = "", y = "Mean Body Mass (g)") +  
  theme_classic() + 
  theme(axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 12),  # Increase the size of x-axis text
        legend.position = "none") +
  scale_fill_manual(values = c("Small-bodied\nFlocks" = "#CD5B45", 
                               "Large-bodied\nCanopy Flocks" = "#027d68", 
                               "Large-bodied\nUndergrowth Flocks" = "#4682B4"))

ggsave(filename = "./outputs/plots/bodymass_boxplot.jpg", plot = bodymass_boxplot, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")


# For flock NMDS ======================================================================================================

flock_nmds = 
  ggplot(nmds_dat, aes(x = MDS1, y = MDS2, color = FlockType)) +
  coord_cartesian(xlim = c(-2,2.5), ylim = c(-2,2)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = FlockType), type = "t", linewidth = 1, level = 0.95) +
  theme_classic() +
  scale_color_manual(values = c("#027d68", "#4682B4", "#CD5B45")) +
  scale_fill_manual(values = c("#027d68", "#4682B4", "#CD5B45")) +
  geom_point(size = 1.5, alpha = .9)

ggsave(filename = "./outputs/plots/flock_nmds.jpg", plot = flock_nmds, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")


# Groupsize vs. VP Beta Regression ===================================================================================

groupsize_beta = 
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

ggsave(filename = "./outputs/plots/groupsize_beta.jpg", plot = groupsize_beta, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")

# Leadership proportion vs. VP Beta GLM ==============================================================================

leadership_beta = 
  ggplot(weighted_vocal_table2, aes(x = lead_proportion, y = mean_vp)) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  geom_line(aes(y = fitted(beta_leadership)), color = '#02607A', linewidth = 1.5) +
  geom_ribbon(aes(ymin = predict_leadership$fit - 1.96 * predict_leadership$se.fit,
                  ymax = predict_leadership$fit + 1.96 * predict_leadership$se.fit),
              alpha = 0.3, fill = "#0387AD") +
  labs(x = "Proportion of Flocks Led", y = "Mean Vocal Participation") +
  theme_classic() +
  theme(plot.title = element_text(size=10), axis.title = element_text(size = 20))

ggsave(filename = "./outputs/plots/leadership_beta.jpg", plot = leadership_beta, device = "jpeg", dpi = 300, width = 9, height = 6, units = "in")

# Degree Centrality vs. VP Beta GLM =================================================================================

y_limits <- c(-.12, 1.05)  # Adjust these values as needed
x_limits <- c(0.25, 1)  # Adjust these values as needed

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
  scale_x_continuous(limits = c(0.3,1))

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
  scale_x_continuous(limits = c(0.3,1))

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
  scale_x_continuous(limits = c(0.3,1))

# aa + bb + cc + plot_layout(ncol = 1)
centrality_beta = 
  aa + bb + plot_layout(ncol = 2) + plot_annotation(caption = "Degree Centrality", 
                                                  theme = theme(plot.caption = element_text(hjust = 0.5, size=15)))


ggsave(filename = "./outputs/plots/centrality_beta.jpg", plot = centrality_beta, device = "jpeg", dpi = 300, width = 8, height = 4, units = "in")


# Proportion in Flocks vs. VP GLM ==========================================================================================================================

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
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) +
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
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) +
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
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) +
  scale_x_continuous(limits = c(0,1))


# ap + bp + cp + plot_layout(ncol = 1)
proportion_beta = 
  ap + bp + plot_layout(ncol = 2) + plot_annotation(caption = "Proportion in Flocks", 
                                                  theme = theme(plot.caption = element_text(hjust = 0.5, size=15)))

ggsave(filename = "./outputs/plots/proportion_beta.jpg", plot = proportion_beta, device = "jpeg", dpi = 300, width = 8, height = 4, units = "in")


# Species PCA ==========================================================================================================================

pca = ggplot(pca_scores, aes(x = PC1, y = PC2, label = species_ID, color = leader_yn, fill = leader_yn)) +
  geom_point(size = 3, shape = 21, stroke = 1.2, alpha = 0.7) +
  #geom_text(vjust = -0.5, size = 3) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  labs(title = "", x = "PC1", y = "PC2", color = "Levels", fill = "Levels") +
  theme_classic()

ggsave(filename = "./outputs/plots/species_pca.jpg", plot = pca, device = "jpeg", dpi = 300, width = 8, height = 4, units = "in")




####################################################################################################################
####################################################################################################################
##################################### Creating Coefficient Tables ##################################################
####################################################################################################################
####################################################################################################################



#S Centrality flocks=========================================================================================================

coefficients_s <- summary(beta_s)$coefficients$cond
intercept_s_logit <- coefficients_s["(Intercept)", "Estimate"]
slope_s_logit <- coefficients_s["s_degree", "Estimate"]
se_s <- coefficients_s["s_degree", 'Std. Error']

s_slope_lower <- slope_s_logit - 1.96 * se_s
s_slope_upper <- slope_s_logit + 1.96 * se_s

marginal_effect_s <- plogis(intercept_s_logit + slope_s_logit) - plogis(intercept_s_logit)
marginal_effect_s_lower <- plogis(intercept_s_logit + s_slope_lower) - plogis(intercept_s_logit)
marginal_effect_s_upper <- plogis(intercept_s_logit + s_slope_upper) - plogis(intercept_s_logit)
z_s <- coefficients_s["s_degree", 'z value']

s_coefficients <- data.frame(
  `Slope` = marginal_effect_s,
  `Lower.95.CI` = marginal_effect_s_lower,
  `Upper.95.CI` = marginal_effect_s_upper,
  `z.value` = z_s
)

write.csv(s_coefficients, "./outputs/tables/s_coefficients.csv", row.names = FALSE)

#LC Centrality flocks=========================================================================================================

coefficients_lc <- summary(beta_lc)$coefficients$cond
intercept_lc_logit <- coefficients_lc["(Intercept)", "Estimate"]
slope_lc_logit <- coefficients_lc["lc_degree", "Estimate"]
se_lc <- coefficients_lc["lc_degree", 'Std. Error']

lc_slope_lower <- slope_lc_logit - 1.96 * se_lc
lc_slope_upper <- slope_lc_logit + 1.96 * se_lc

marginal_effect_lc <- plogis(intercept_lc_logit + slope_lc_logit) - plogis(intercept_lc_logit)
marginal_effect_lc_lower <- plogis(intercept_lc_logit + lc_slope_lower) - plogis(intercept_lc_logit)
marginal_effect_lc_upper <- plogis(intercept_lc_logit + lc_slope_upper) - plogis(intercept_lc_logit)
z_lc <- coefficients_lc["lc_degree", 'z value']

lc_coefficients <- data.frame(
  `Slope` = marginal_effect_lc,
  `Lower.95.CI` = marginal_effect_lc_lower,
  `Upper.95.CI` = marginal_effect_lc_upper,
  `z.value` = z_lc
)

write.csv(lc_coefficients, "./outputs/tables/lc_coefficients.csv", row.names = FALSE)

#S Proportion flocks=========================================================================================================

coefficients_sp <- summary(beta_sp)$coefficients$cond
intercept_sp_logit <- coefficients_sp["(Intercept)", "Estimate"]
slope_sp_logit <- coefficients_sp["percent_flocks", "Estimate"]
se_sp <- coefficients_sp["percent_flocks", 'Std. Error']

sp_slope_lower <- slope_sp_logit - 1.96 * se_sp
sp_slope_upper <- slope_sp_logit + 1.96 * se_sp

marginal_effect_sp <- plogis(intercept_sp_logit + slope_sp_logit) - plogis(intercept_sp_logit)
marginal_effect_sp_lower <- plogis(intercept_sp_logit + sp_slope_lower) - plogis(intercept_sp_logit)
marginal_effect_sp_upper <- plogis(intercept_sp_logit + sp_slope_upper) - plogis(intercept_sp_logit)
z_sp <- coefficients_sp["percent_flocks", 'z value']

sp_coefficients <- data.frame(
  `Slope` = marginal_effect_sp,
  `Lower.95.CI` = marginal_effect_sp_lower,
  `Upper.95.CI` = marginal_effect_sp_upper,
  `z.value` = z_sp
)

write.csv(sp_coefficients, "./outputs/tables/sp_coefficients.csv", row.names = FALSE)

#LC Proportion flocks=========================================================================================================

coefficients_lcp <- summary(beta_lcp)$coefficients$cond
intercept_lcp_logit <- coefficients_lcp["(Intercept)", "Estimate"]
slope_lcp_logit <- coefficients_lcp["percent_flocks", "Estimate"]
se_lcp <- coefficients_lcp["percent_flocks", 'Std. Error']

lcp_slope_lower <- slope_lcp_logit - 1.96 * se_lcp
lcp_slope_upper <- slope_lcp_logit + 1.96 * se_lcp

marginal_effect_lcp <- plogis(intercept_lcp_logit + slope_lcp_logit) - plogis(intercept_lcp_logit)
marginal_effect_lcp_lower <- plogis(intercept_lcp_logit + lcp_slope_lower) - plogis(intercept_lcp_logit)
marginal_effect_lcp_upper <- plogis(intercept_lcp_logit + lcp_slope_upper) - plogis(intercept_lcp_logit)
z_lcp <- coefficients_lcp["percent_flocks", 'z value']

lcp_coefficients <- data.frame(
  `Slope` = marginal_effect_lcp,
  `Lower.95.CI` = marginal_effect_lcp_lower,
  `Upper.95.CI` = marginal_effect_lcp_upper,
  `z.value` = z_lcp
)

write.csv(lcp_coefficients, "./outputs/tables/lcp_coefficients.csv", row.names = FALSE)

#Leadership vs. VP flocks=========================================================================================================

coefficients_lead <- summary(beta_leadership)$coefficients$cond
intercept_lead_logit <- coefficients_lead["(Intercept)", "Estimate"]
slope_lead_logit <- coefficients_lead["lead_proportion", "Estimate"]
se_lead <- coefficients_lead["lead_proportion", 'Std. Error']

lead_slope_lower <- slope_lead_logit - 1.96 * se_lead
lead_slope_upper <- slope_lead_logit + 1.96 * se_lead

marginal_effect_lead <- plogis(intercept_lead_logit + slope_lead_logit) - plogis(intercept_lead_logit)
marginal_effect_lead_lower <- plogis(intercept_lead_logit + lead_slope_lower) - plogis(intercept_lead_logit)
marginal_effect_lead_upper <- plogis(intercept_lead_logit + lead_slope_upper) - plogis(intercept_lead_logit)
z_lead <- coefficients_lead["lead_proportion", 'z value']

lead_coefficients <- data.frame(
  `Slope` = marginal_effect_lead,
  `Lower.95.CI` = marginal_effect_lead_lower,
  `Upper.95.CI` = marginal_effect_lead_upper,
  `z.value` = z_lead
)

write.csv(lead_coefficients, "./outputs/tables/lead_coefficients.csv", row.names = FALSE)

#Groupsize vs. VP flocks=========================================================================================================

coefficients_groupsize <- summary(beta_groupsize)$coefficients$cond
intercept_groupsize_logit <- coefficients_groupsize["(Intercept)", "Estimate"]
slope_groupsize_logit <- coefficients_groupsize["avg_group_size", "Estimate"]
se_groupsize <- coefficients_groupsize["avg_group_size", 'Std. Error']

groupsize_slope_lower <- slope_groupsize_logit - 1.96 * se_groupsize
groupsize_slope_upper <- slope_groupsize_logit + 1.96 * se_groupsize

marginal_effect_groupsize <- plogis(intercept_groupsize_logit + slope_groupsize_logit) - plogis(intercept_groupsize_logit)
marginal_effect_groupsize_lower <- plogis(intercept_groupsize_logit + groupsize_slope_lower) - plogis(intercept_groupsize_logit)
marginal_effect_groupsize_upper <- plogis(intercept_groupsize_logit + groupsize_slope_upper) - plogis(intercept_groupsize_logit)
z_groupsize <- coefficients_groupsize["avg_group_size", 'z value']

groupsize_coefficients <- data.frame(
  `Slope` = marginal_effect_groupsize,
  `Lower.95.CI` = marginal_effect_groupsize_lower,
  `Upper.95.CI` = marginal_effect_groupsize_upper,
  `z.value` = z_groupsize
)

write.csv(groupsize_coefficients, "./outputs/tables/groupsize_coefficients.csv", row.names = FALSE)

#Stress for NMDS Plot=========================================================================================================

stress <- nmds$stress
write.csv(stress, "./outputs/tables/nmds_stress.csv", row.names = FALSE)

#Bodymass table=========================================================================================================

write.csv(bodymass_table, "./outputs/tables/bodymass_table.csv", row.names = FALSE)
