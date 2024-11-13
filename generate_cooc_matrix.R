setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")

# Vocal co-occurrence Analysis


# Load necessary libraries
library(dplyr)
library(tidyr)
set.seed(123)

# Read the data
data <- read.csv("vocal_participation.csv")

# Function to clean the vocalization data by removing trailing NAs
clean_vocalization_data <- function(vocalization_data) {
  # Find the last non-NA column for each row and subset the data accordingly
  last_valid_index <- apply(vocalization_data, 1, function(row) max(which(!is.na(row))))
  max_index <- max(last_valid_index)
  vocalization_data <- vocalization_data[, 1:max_index]
  vocalization_data[is.na(vocalization_data)] <- 0
  return(vocalization_data)
}

# Function to calculate z-score of co-occurrences
calculate_z_score_co_occurrences <- function(species_a, species_b, vocalization_data) {
  vocal_a <- unlist(vocalization_data %>% filter(species_ID == species_a) %>% select(-species_ID))
  vocal_b <- unlist(vocalization_data %>% filter(species_ID == species_b) %>% select(-species_ID))
  
  # Actual co-occurrences
  actual_co_occurrences <- sum(vocal_a == 1 & vocal_b == 1)
  
  # Random shuffles and co-occurrences
  shuffle_co_occurrences <- replicate(1000, {
    shuffled_a <- sample(vocal_a)
    shuffled_b <- sample(vocal_b)
    sum(shuffled_a == 1 & shuffled_b == 1)
  })
  
  # Calculate mean and standard deviation of the shuffled co-occurrences
  mean_shuffle_co_occurrences <- mean(shuffle_co_occurrences)
  sd_shuffle_co_occurrences <- sd(shuffle_co_occurrences)
  
  # Calculate z-score
  if (sd_shuffle_co_occurrences != 0) {
    z_score <- (actual_co_occurrences - mean_shuffle_co_occurrences) / sd_shuffle_co_occurrences
  } else {
    z_score <- 0
  }
  
  return(z_score)
}

# Function to calculate co-occurrence matrix for a single flock
calculate_co_occurrence_matrix <- function(flock_data) {
  # Extract relevant columns and clean vocalization data
  vocalization_data <- flock_data %>%
    select(species_ID, starts_with("X"))
  vocalization_data <- clean_vocalization_data(vocalization_data)
  
  # Get unique species IDs
  species_ids <- unique(vocalization_data$species_ID)
  
  # Initialize a matrix to store the results
  co_occurrence_matrix <- matrix(NA, nrow=length(species_ids), ncol=length(species_ids), 
                                 dimnames=list(species_ids, species_ids))
  
  # Calculate the pairwise z-score co-occurrence values
  for (i in seq_along(species_ids)) {
    for (j in seq_along(species_ids)) {
      if (i != j) {
        co_occurrence_matrix[i, j] <- calculate_z_score_co_occurrences(species_ids[i], species_ids[j], vocalization_data)
      } else {
        co_occurrence_matrix[i, j] <- NA  # Self-comparison is not needed
      }
    }
  }
  
  # Convert the matrix to a data frame for better readability
  co_occurrence_df <- as.data.frame(co_occurrence_matrix)
  
  return(co_occurrence_df)
}

# Split the data by flock_ID
flocks <- split(data, data$flock_ID)

# Calculate and save co-occurrence matrices for each flock
for (flock_id in names(flocks)) {
  flock_data <- flocks[[flock_id]]
  co_occurrence_df <- calculate_co_occurrence_matrix(flock_data)
  
  # Save the result to a CSV file
  file_name <- paste0("co_occurrence_matrix_flock_", flock_id, ".csv")
  write.csv(co_occurrence_df, file_name, row.names = TRUE)
}

#============================================================================================================================

# Making the final Matrix


# Load necessary library
library(dplyr)

# Function to read a matrix from a file and convert it to a dataframe
read_co_occurrence_matrix <- function(file) {
  df <- read.csv(file, row.names = 1, check.names = FALSE)
  df <- as.data.frame(as.matrix(df))
  df
}

# Initialize an empty list to hold the matrices
matrices <- list()

# Directory containing the CSV files
directory <- "C:/Users/ronit/OneDrive/Desktop/Thesis stuff"

# List of files to read
files <- list.files(directory, pattern = "co_occurrence_matrix_flock_\\d+\\.csv", full.names = TRUE)

# Read each matrix and store it in the list
for (file in files) {
  matrix <- read_co_occurrence_matrix(file)
  matrices[[file]] <- matrix
}

# Get all unique species names across all matrices
all_species <- unique(unlist(lapply(matrices, rownames)))

# Initialize the final matrix with NA values
final_matrix <- matrix(NA, nrow = length(all_species), ncol = length(all_species))
rownames(final_matrix) <- all_species
colnames(final_matrix) <- all_species

# Count matrix to keep track of how many times each species pair occurs
count_matrix <- matrix(0, nrow = length(all_species), ncol = length(all_species))
rownames(count_matrix) <- all_species
colnames(count_matrix) <- all_species

# Iterate through each matrix and update the final matrix
for (matrix in matrices) {
  species_in_matrix <- rownames(matrix)
  for (i in 1:length(species_in_matrix)) {
    for (j in 1:length(species_in_matrix)) {
      species_i <- species_in_matrix[i]
      species_j <- species_in_matrix[j]
      final_matrix[species_i, species_j] <- ifelse(is.na(final_matrix[species_i, species_j]),
                                                   matrix[species_i, species_j],
                                                   final_matrix[species_i, species_j] + matrix[species_i, species_j])
      count_matrix[species_i, species_j] <- count_matrix[species_i, species_j] + 1
    }
  }
}

# Average the values by dividing by the count matrix, avoiding division by zero
final_matrix <- final_matrix / count_matrix

# # Convert values where the count is 4 or less to 0
# final_matrix[count_matrix <= 4] <- 0
# final_matrix[final_matrix == 0] <- NA



# Save the final matrix to a CSV file
write.csv(final_matrix, file = "consolidated_co_occurrence_matrix.csv", row.names = TRUE)

co_occ_values =  data.frame(rowMeans(final_matrix, na.rm = T))
co_occ_values$species_ID = rownames(co_occ_values)
rownames(co_occ_values) = NULL
colnames(co_occ_values) = c("vocal_co_occ", "species_ID")
co_occ_values$vocal_co_occ[is.nan(co_occ_values$vocal_co_occ)] <- NA

write.csv(co_occ_values, "Vocal_co_occurrence_values.csv")


#==========================================================================================================================
#==========================================================================================================================

# Vocal Co-occurrence Analysis

#==========================================================================================================================
#==========================================================================================================================


library(reshape2)

# Subset the S matrix

# Ensure the column names exist in final_matrix
s_filtered_columns <- colnames(S_flocks2)[colnames(S_flocks2) %in% colnames(final_matrix)]

# Subset the matrix
s_cooc_matrix <- final_matrix[s_filtered_columns, s_filtered_columns]

# Remove the row and column
names_to_remove_s <- c("yuhina_flavicollis", "mixornis_gularis", "actinodura_cyanouroptera", "phylloscopus_inornatus")

# Find the indices of the rows and columns to keep
rows_to_keep_s <- setdiff(rownames(s_cooc_matrix), names_to_remove_s)
cols_to_keep_s <- setdiff(colnames(s_cooc_matrix), names_to_remove_s)

# Subset the matrix
s_cooc_matrix <- s_cooc_matrix[rows_to_keep_s, cols_to_keep_s]

s_melted_matrix <- melt(s_cooc_matrix)

# Create the heatmap
ggplot(data = s_melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#No labels
ggplot(data = s_melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "", x = "", y = "") + theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_text(size = 16), legend.text = element_text(size = 14))

#==========================================================================================================================

lc_filtered_columns <- colnames(LC_flocks2)[colnames(LC_flocks2) %in% colnames(final_matrix)]

# Subset the matrix
lc_cooc_matrix <- final_matrix[lc_filtered_columns, lc_filtered_columns]

# Remove the row and column
names_to_remove_lc <- c("harpactes_erythrocephalus", "tephrodornis_virgatus", 
                     "lalage_melaschistos", "heterophasia_picaoides", "dicrurus_paradiseus", "yuhina_flavicollis",
                     "alophoixus_flaveolus", "actinodura_egertoni")

# Find the indices of the rows and columns to keep
rows_to_keep_lc <- setdiff(rownames(lc_cooc_matrix), names_to_remove_lc)
cols_to_keep_lc <- setdiff(colnames(lc_cooc_matrix), names_to_remove_lc)

# Subset the matrix
lc_cooc_matrix <- lc_cooc_matrix[rows_to_keep_lc, cols_to_keep_lc]

lc_melted_matrix <- melt(lc_cooc_matrix)

# Create the heatmap
ggplot(data = lc_melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#No labels
ggplot(data = lc_melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "", x = "", y = "") + theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_text(size = 16), legend.text = element_text(size = 14))

#==========================================================================================================================

# Ensure the column names exist in final_matrix
lu_filtered_columns <- colnames(LU_flocks2)[colnames(LU_flocks2) %in% colnames(final_matrix)]

# Subset the matrix
lu_cooc_matrix <- final_matrix[lu_filtered_columns, lu_filtered_columns]


# Remove the row and column
index_to_remove <- which(rownames(lu_cooc_matrix) == "alophoixus_flaveolus")
lu_cooc_matrix <- lu_cooc_matrix[-index_to_remove, -index_to_remove]


lu_melted_matrix <- melt(lu_cooc_matrix)

# Create the heatmap
ggplot(data = lu_melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#No labels
ggplot(data = lu_melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "", x = "", y = "") + theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_text(size = 16), legend.text = element_text(size = 14))



#==========================================================================================================================
#==========================================================================================================================
#==========================================================================================================================


z_scores_s <- as.vector(s_cooc_matrix)

# Remove NA values from z_scores
z_scores_s <- z_scores_s[!is.na(z_scores_s)]

# Create a data frame for ggplot
df_s <- data.frame(z_scores_s = z_scores_s)

# Generate a sequence of values from the range of z_scores
x_values <- rnorm(1000000, mean = 0, sd = 0.5)

# Create a data frame for the normal distribution
df_normal <- data.frame(x_values = x_values)

# Plot the z-score distribution and the normal distribution
sz = ggplot() +
  geom_density(data = df_s, aes(x = z_scores_s), color = "blue", fill = "blue", alpha = 0.3) +
  geom_density(data = df_normal, aes(x = x_values), color = "red", fill = "red", alpha = 0.3) +
  labs(title = "a) Small-bodied Flocks", x = "", y = "") +
  theme_classic() + scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(limits = c(0,1.4))

#==========================================================================================================================

z_scores_lc <- as.vector(lc_cooc_matrix)

# Remove NA values from z_scores
z_scores_lc <- z_scores_lc[!is.na(z_scores_lc)]

# Create a data frame for ggplot
df_lc <- data.frame(z_scores_lc = z_scores_lc)

lcz = ggplot() +
  geom_density(data = df_lc, aes(x = z_scores_lc), color = "blue", fill = "blue", alpha = 0.3) +
  geom_density(data = df_normal, aes(x = x_values), color = "red", fill = "red", alpha = 0.3) +
  labs(title = "b) Large-bodied Canopy Flocks", x = "", y = "Density") +
  theme_classic() + scale_x_continuous(limits = c(-2,2)) +
  theme(axis.title = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,1.4))

#==========================================================================================================================

z_scores_lu <- as.vector(lu_cooc_matrix)

# Remove NA values from z_scores
z_scores_lu <- z_scores_lu[!is.na(z_scores_lu)]

# Create a data frame for ggplot
df_lu <- data.frame(z_scores_lu = z_scores_lu)

luz = ggplot() +
  geom_density(data = df_lu, aes(x = z_scores_lu), color = "blue", fill = "blue", alpha = 0.3) +
  geom_density(data = df_normal, aes(x = x_values), color = "red", fill = "red", alpha = 0.3) +
  labs(title = "c) Large-bodied Undergrowth Flocks", x = "Z-Scores", y = "") +
  theme_classic() + scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(0,1.4))

#==========================================================================================================================


sz + lcz + luz + plot_layout(ncol = 1)

t.test(df_s$z_scores_s, df_normal$x_values)
t.test(df_lc$z_scores_lc, df_normal$x_values)
t.test(df_lu$z_scores_lu, df_normal$x_values)

mean(df_s$z_scores_s)
mean(df_lc$z_scores_lc)
mean(df_lu$z_scores_lu)

sum(df_s$z_scores_s < 0)
sum(df_lc$z_scores_lc < 0)
sum(df_lu$z_scores_lu < 0)

#==========================================================================================================================

hist(co_occ_values$vocal_co_occ)




