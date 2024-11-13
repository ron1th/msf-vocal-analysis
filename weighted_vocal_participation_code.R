setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")

#==========================================================================================================================
library(igraph)
library(igraphdata)
library(ggplot2)
library(dplyr)
library(tidyverse)

lc_plot$mean_weighted_vp = lc_plot$mean_weighted_vp + 0.0000001

# Fit the beta regression model
beta_sw <- glmmTMB(mean_weighted_vp ~ s_degree, data = s_plot, family = beta_family(link = "logit"))
beta_lcw <- glmmTMB(mean_weighted_vp ~ lc_degree, data = lc_plot, family = beta_family(link = "logit"))
beta_luw <- glmmTMB(mean_weighted_vp ~ lu_degree, data = lu_plot, family = beta_family(link = "logit"))


summary(beta_sw) # **
summary(beta_lcw) #
summary(beta_luw) # *


hist(resid(beta_s))
hist(resid(beta_lc))
hist(resid(beta_lu))



# Get predicted values and standard errors
predict_sw <- predict(beta_sw, type = "response", se.fit = TRUE)
predict_lcw <- predict(beta_lcw, type = "response", se.fit = TRUE)
predict_luw <- predict(beta_luw, type = "response", se.fit = TRUE)


# Plot with ggplot2 --- FITTED ERRORBAR SHADING TO ALL!!!!!!!!!!!!!!!!!!!!!


aa2 = ggplot(s_plot, aes(x = s_degree, y = mean_weighted_vp)) +
  geom_point(size = 3, shape = 21, fill = "#FF7F50", colour = "#CD5B45", stroke = 1.2, alpha = .7) +
  geom_line(aes(y = fitted(beta_sw)), colour = "#CD5B45", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_sw$fit - 1.96 * predict_sw$se.fit,
                  ymax = predict_sw$fit + 1.96 * predict_sw$se.fit),
              alpha = 0.3, fill = "#FF7F50") +
  labs(x = NULL, y = "Weighted Mean Vocal Participation") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("a) Small-bodied Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) +
  xlim(x_limits)

bb2 = ggplot(lc_plot, aes(x = lc_degree, y = mean_weighted_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#02a890", colour = "#027d68", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lcw)), colour = "#027d68", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_lcw$fit - 1.96 * predict_lcw$se.fit,
                  ymax = predict_lcw$fit + 1.96 * predict_lcw$se.fit),
              alpha = 0.3, fill = "#02a890") +
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("b) Large-bodied Canopy Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + 
  xlim(x_limits)

cc2 = ggplot(lu_plot, aes(x = lu_degree, y = mean_weighted_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#87CEEB", colour = "#4682B4", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_luw)), colour = "#4682B4", linewidth = 1.8) + 
  geom_ribbon(aes(ymin = predict_luw$fit - 1.96 * predict_luw$se.fit,
                  ymax = predict_luw$fit + 1.96 * predict_luw$se.fit),
              alpha = 0.3, fill = "#87CEEB") +
  labs(x = 'Degree Centrality', y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("c) Large-bodied Undergrowth Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits) + 
  xlim(x_limits)



aa2 + bb2 + cc2 + plot_layout(ncol = 1)
aa2 + bb2 + plot_layout(ncol = 2) + plot_annotation(caption = "Degree Centrality", 
                                                  theme = theme(plot.caption = element_text(hjust = 0.5, size=15)))


#==========================================================================================================================

inv_logit <- function(x) exp(x) / (1 + exp(x))

# Find coefficients on the linear scale:
coefficients_sw <- summary(beta_sw)$coefficients$cond
intercept_sw_logit <- coefficients_sw["(Intercept)", "Estimate"]
slope_sw_logit <- coefficients_sw["s_degree", "Estimate"]

# Transform the intercept and slope for sw flocks
intercept_sw <- inv_logit(intercept_sw_logit)
slope_sw <- inv_logit(intercept_sw_logit + slope_sw_logit) - inv_logit(intercept_sw_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_sw, "\n")
cat("Slope on the linear scale:", slope_sw, "\n")

#==========================================================================================================================

# LCW flocks
coefficients_lcw <- summary(beta_lcw)$coefficients$cond
intercept_lcw_logit <- coefficients_lcw["(Intercept)", "Estimate"]
slope_lcw_logit <- coefficients_lcw["lc_degree", "Estimate"]

# Transform the intercept and slope for lcw flocks
intercept_lcw <- inv_logit(intercept_lcw_logit)
slope_lcw <- inv_logit(intercept_lcw_logit + slope_lcw_logit) - inv_logit(intercept_lcw_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_lcw, "\n")
cat("Slope on the linear scale:", slope_lcw, "\n")

#==========================================================================================================================

# LUW flocks
coefficients_luw <- summary(beta_luw)$coefficients$cond
intercept_luw_logit <- coefficients_luw["(Intercept)", "Estimate"]
slope_luw_logit <- coefficients_luw["lu_degree", "Estimate"]

# Transform the intercept and slope for luw flocks
intercept_luw <- inv_logit(intercept_luw_logit)
slope_luw <- inv_logit(intercept_luw_logit + slope_luw_logit) - inv_logit(intercept_luw_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_luw, "\n")
cat("Slope on the linear scale:", slope_luw, "\n")


#==========================================================================================================================

# Fitting GLM to vocal participation data

vocal_table$vocal_participation = vocal_table$vocal_participation * 0.999999
vocal_table$vocal_participation = vocal_table$vocal_participation + 0.000000001

beta_vp = glmmTMB(vocal_participation ~ count, data = vocal_table, family = beta_family(link = "logit"))
predict_vp <- predict(beta_vp, type = "response", se.fit = TRUE)

ggplot(vocal_table, aes(x = count, y = vocal_participation)) +
  geom_point() +
  geom_line(aes(y = fitted(beta_vp)), color = 'blue', size = 1)

#==========================================================================================================================


# Boxplot for Leaders vs. Non-Leaders

ggplot(weighted_vocal_table2, aes(x = factor(leader_yn), y = mean_weighted_vp, fill = factor(leader_yn))) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  geom_violin(trim = FALSE, alpha = 0.3, width = 0.4) +    # Decrease the width of the violin plot
  labs(x = "", y = "Weighted Mean Vocal Participation") + 
  theme_classic() + 
  theme(legend.position = "none", axis.title = element_text(size = 16), axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = c("#0387ad", "#02a890")) + scale_y_continuous(limits = c(-.15,1.08))

weighted_vocal_table2 %>% filter(leader_yn == "Leaders") %>% summarise(mean(mean_weighted_vp))
weighted_vocal_table2 %>% filter(leader_yn == "Non-leaders") %>% summarise(mean(mean_weighted_vp))

#==========================================================================================================================

# Binomial GLM for group size vs. VP
ggplot(data = filtered_table, aes(x = avg_group_size, y = mean_weighted_vp)) +
  geom_point(size = 2.2) + theme_classic() +
  labs(x = "Mean Intraspecific Group Size", y = "Weighted Mean Vocal Participation") +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, linewidth = 1.5) +
  theme(axis.title = element_text(size = 20)) + scale_y_continuous(limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = c(0:9))

ggplot(data = filtered_table, aes(x = avg_group_size, y = mean_weighted_vp)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, 
              linewidth = 1.8, colour = "#02607A", fill = "#87CEEB") +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) + 
  theme_classic() +
  labs(x = "Mean Intraspecific Group Size", y = "Weighted Mean Vocal Participation") +
  theme(axis.title = element_text(size = 18)) + 
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = c(0:9))

# Fit the binomial regression model
binomial_weighted_groupsize <- glm(mean_weighted_vp ~ avg_group_size, data = filtered_table, family = binomial)

# Summarize the model
summary(binomial_weighted_groupsize)

# Extract coefficients (slope and intercept)
weighted_groupsize_coefficients <- summary(binomial_weighted_groupsize)$coefficients

# Intercept (linear/odds ratio)
exp(weighted_groupsize_coefficients[1, 1])

# Slope (linear/odds ratio)
exp(weighted_groupsize_coefficients[2, 1])

# P-values
weighted_groupsize_coefficients[, 4]


#Beta regression instead of Binomial
filtered_table$mean_weighted_vp = filtered_table$mean_weighted_vp + 0.0000001
beta_wgroupsize <- glmmTMB(mean_weighted_vp ~ avg_group_size, data = filtered_table, family = beta_family(link = "logit"))
predict_wgroupsize <- predict(beta_wgroupsize, type = "response", se.fit = TRUE)

ggplot(filtered_table, aes(x = avg_group_size, y = mean_weighted_vp)) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  geom_line(aes(y = fitted(beta_wgroupsize)), color = '#02607A', linewidth = 1.5) +
  geom_ribbon(aes(ymin = predict_wgroupsize$fit - 1.96 * predict_wgroupsize$se.fit,
                  ymax = predict_wgroupsize$fit + 1.96 * predict_wgroupsize$se.fit),
              alpha = 0.3, fill = "#0387AD") +
  labs(x = "Mean Intraspecific Group Size", y = "Weighted Mean Vocal Participation") +
  theme_classic() +
  theme(plot.title = element_text(size=10), axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = c(0:9))

# Coefficients
coefficients_wgroupsize <- summary(beta_wgroupsize)$coefficients$cond
intercept_wgroupsize_logit <- coefficients_wgroupsize["(Intercept)", "Estimate"]
slope_wgroupsize_logit <- coefficients_wgroupsize["avg_group_size", "Estimate"]

plogis(0.2838-0.0620)

# Transform the intercept and slope for groupsize
intercept_wgroupsize <- inv_logit(intercept_wgroupsize_logit)
slope_wgroupsize <- inv_logit(intercept_wgroupsize_logit + slope_wgroupsize_logit) - inv_logit(intercept_wgroupsize_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_wgroupsize, "\n")
cat("Slope on the linear scale:", slope_wgroupsize, "\n")

summary(beta_wgroupsize)

#==========================================================================================================================
#==========================================================================================================================

# Fit the beta regression model
beta_spw <- glmmTMB(mean_weighted_vp ~ percent_flocks, data = s_plot, family = beta_family(link = "logit"))
beta_lcpw <- glmmTMB(mean_weighted_vp ~ percent_flocks, data = lc_plot, family = beta_family(link = "logit"))
beta_lupw <- glmmTMB(mean_weighted_vp ~ percent_flocks, data = lu_plot, family = beta_family(link = "logit"))

summary(beta_spw) # *
summary(beta_lcpw) #
summary(beta_lupw) # *

plogis(3.6973-1.8647)

hist(resid(beta_spw))
hist(resid(beta_lcpw))
hist(resid(beta_lupw))

# Get predicted values and standard errors
predict_spw <- predict(beta_spw, type = "response", se.fit = TRUE)
predict_lcpw <- predict(beta_lcpw, type = "response", se.fit = TRUE)
predict_lupw <- predict(beta_lupw, type = "response", se.fit = TRUE)

# Plot with ggplot2 --- FITTED ERRORBAR SHADING TO ALL!!!!!!!!!!!!!!!!!!!!!

y_limits2 <- c(-.12, 1.05)

apw = ggplot(s_plot, aes(x = percent_flocks, y = mean_weighted_vp)) +
  geom_point(size = 3, shape = 21, fill = "#FF7F50", colour = "#CD5B45", stroke = 1.2, alpha = .7) +
  geom_line(aes(y = fitted(beta_spw)), colour = "#CD5B45", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_spw$fit - 1.96 * predict_spw$se.fit,
                  ymax = predict_spw$fit + 1.96 * predict_spw$se.fit),
              alpha = 0.3, fill = "#FF7F50") +
  labs(x = NULL, y = "Weighted Mean Vocal Participation") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("a) Small-bodied Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits2) +
  scale_x_continuous(limits = c(0,1))

# Plot with bb colors
bpw = ggplot(lc_plot, aes(x = percent_flocks, y = mean_weighted_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#02a890", colour = "#027d68", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lcpw)), colour = "#027d68", linewidth = 1.8) +
  geom_ribbon(aes(ymin = predict_lcpw$fit - 1.96 * predict_lcpw$se.fit,
                  ymax = predict_lcpw$fit + 1.96 * predict_lcpw$se.fit),
              alpha = 0.3, fill = "#02a890") +
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("b) Large-bodied Canopy Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits2) +
  scale_x_continuous(limits = c(0,1))

# Plot with cc colors
cpw = ggplot(lu_plot, aes(x = percent_flocks, y = mean_weighted_vp)) + 
  geom_point(size = 3, shape = 21, fill = "#87CEEB", colour = "#4682B4", stroke = 1.2, alpha = .7) + 
  geom_line(aes(y = fitted(beta_lupw)), colour = "#4682B4", linewidth = 1.8) + 
  geom_ribbon(aes(ymin = predict_lupw$fit - 1.96 * predict_lupw$se.fit,
                  ymax = predict_lupw$fit + 1.96 * predict_lupw$se.fit),
              alpha = 0.3, fill = "#87CEEB") +
  labs(x = 'Proportion in Flocks', y = NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(), plot.title = element_text(size=10)) +
  ggtitle("c) Large-bodied Undergrowth Flocks") +
  scale_y_continuous(expand = c(0.05, 0.05), limits = y_limits2) +
  scale_x_continuous(limits = c(0,1))

apw + bpw + cpw + plot_layout(ncol = 1)
apw + bpw + plot_layout(ncol = 2) + plot_annotation(caption = "Proportion in Flocks", 
                                                    theme = theme(plot.caption = element_text(hjust = 0.5, size=15)))


#==========================================================================================================================

inv_logit <- function(x) exp(x) / (1 + exp(x))

# Find coefficients on the linear scale:
coefficients_spw <- summary(beta_spw)$coefficients$cond
intercept_spw_logit <- coefficients_spw["(Intercept)", "Estimate"]
slope_spw_logit <- coefficients_spw["percent_flocks", "Estimate"]

# Transform the intercept and slope for sp flocks
intercept_spw <- inv_logit(intercept_spw_logit)
slope_spw <- inv_logit(intercept_spw_logit + slope_spw_logit) - inv_logit(intercept_spw_logit)

# Print the results
summary(beta_spw)
cat("Intercept on the linear scale:", intercept_spw, "\n")
cat("Slope on the linear scale:", slope_spw, "\n")

#==========================================================================================================================

# LCP flocks
coefficients_lcpw <- summary(beta_lcpw)$coefficients$cond
intercept_lcpw_logit <- coefficients_lcpw["(Intercept)", "Estimate"]
slope_lcpw_logit <- coefficients_lcpw["percent_flocks", "Estimate"]

# Transform the intercept and slope for lcp flocks
intercept_lcpw <- inv_logit(intercept_lcpw_logit)
slope_lcpw <- inv_logit(intercept_lcpw_logit + slope_lcpw_logit) - inv_logit(intercept_lcpw_logit)

# Print the results
summary(beta_lcpw)
cat("Intercept on the linear scale:", intercept_lcpw, "\n")
cat("Slope on the linear scale:", slope_lcpw, "\n")

#==========================================================================================================================

# LUP flocks
coefficients_lupw <- summary(beta_lupw)$coefficients$cond
intercept_lupw_logit <- coefficients_lupw["(Intercept)", "Estimate"]
slope_lupw_logit <- coefficients_lupw["percent_flocks", "Estimate"]

# Transform the intercept and slope for lup flocks
intercept_lupw <- inv_logit(intercept_lupw_logit)
slope_lupw <- inv_logit(intercept_lupw_logit + slope_lupw_logit) - inv_logit(intercept_lupw_logit)

# Print the results
summary(beta_lupw)
cat("Intercept on the linear scale:", intercept_lupw, "\n")
cat("Slope on the linear scale:", slope_lupw, "\n")

#==========================================================================================================================
#==========================================================================================================================


# Plot with ggplot2
ggplot(data = weighted_vocal_table2, aes(x = lead_proportion, y = mean_weighted_vp)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), 
              se = TRUE, colour = "#02607A", fill = "#87CEEB", linewidth = 1.8) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  theme_classic() +
  labs(x = "Proportion of Flocks Led", y = "Weighted Mean Vocal Participation") +
  theme(axis.title = element_text(size = 18))


# Fit a GLM with a beta distribution -- weighted VP vs. Leadership p/n graph

weighted_vocal_table2$mean_weighted_vp = weighted_vocal_table2$mean_weighted_vp + 0.0000001
beta_wleadership = glmmTMB(mean_weighted_vp ~ lead_proportion, data = weighted_vocal_table2, family = beta_family(link = "logit"))
predict_wleadership <- predict(beta_wleadership, type = "response", se.fit = TRUE)


ggplot(weighted_vocal_table2, aes(x = lead_proportion, y = mean_weighted_vp)) +
  geom_point(size = 2.2, shape = 21, fill = "#0387AD", colour = "#0387AD", stroke = 1.2, alpha = .6) +
  geom_line(aes(y = fitted(beta_wleadership)), color = '#02607A', linewidth = 1.5) +
  geom_ribbon(aes(ymin = predict_wleadership$fit - 1.96 * predict_wleadership$se.fit,
                  ymax = predict_wleadership$fit + 1.96 * predict_wleadership$se.fit),
              alpha = 0.3, fill = "#0387AD") +
  labs(x = "Proportion of Flocks Led", y = "Weighted Mean Vocal Participation") +
  theme_classic() +
  theme(plot.title = element_text(size=10), axis.title = element_text(size = 16))

#Coefficients
coefficients_wlead <- summary(beta_wleadership)$coefficients$cond
intercept_wlead_logit <- coefficients_wlead["(Intercept)", "Estimate"]
slope_wlead_logit <- coefficients_wlead["lead_proportion", "Estimate"]

summary(beta_wleadership)
plogis(5.1338-1.0320)


# Transform the intercept and slope for s flocks
intercept_wlead <- inv_logit(intercept_wlead_logit)
slope_wlead <- inv_logit(intercept_wlead_logit + slope_wlead_logit) - inv_logit(intercept_wlead_logit)

# Print the results
cat("Intercept on the linear scale:", intercept_wlead, "\n")
cat("Slope on the linear scale:", slope_wlead, "\n")


#==========================================================================================================================
# Find slope, intercept, p values
# Fit the binomial regression model
binomial_weighted_leadership <- glm(mean_weighted_vp ~ lead_proportion, data = weighted_vocal_table2, family = binomial)

# Summarize the model
summary(binomial_weighted_leadership)

# Extract coefficients (slope and intercept)
weighted_leadership_coefficients <- summary(binomial_weighted_leadership)$coefficients

# Intercept (linear/odds ratio)
exp(weighted_leadership_coefficients[1, 1])

# Slope (linear/odds ratio)
exp(weighted_leadership_coefficients[2, 1])

# P-values
weighted_leadership_coefficients[, 4]

#==========================================================================================================================

# Boxplot for Leaders vs. Non-Leaders
ggplot(weighted_vocal_table2, aes(x = factor(leader_yn), y = mean_weighted_vp, fill = factor(leader_yn))) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +  # Increase the width of the box plot
  geom_violin(trim = FALSE, alpha = 0.3, width = 0.4) +    # Decrease the width of the violin plot
  labs(x = "", y = "Weighted Mean Vocal Participation") + 
  theme_classic() + 
  theme(legend.position = "none", axis.title = element_text(size = 16), axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = c("#CD5B45", "#0387ad")) +  # Purple and green colors
  scale_y_continuous(limits = c(-0.15, 1.08))

# T-test
t.test(mean_weighted_vp ~ leader_yn, data = weighted_vocal_table2)

