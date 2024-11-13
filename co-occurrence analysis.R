setwd("C:/Users/ronit/OneDrive/Desktop/Thesis stuff")


library(vegan)
library(dplyr)

# Function to compute co-occurrence matrix
compute_cooccurrence <- function(mat) {
  cooc_mat <- t(mat) %*% mat
  diag(cooc_mat) <- 0
  return(cooc_mat)
}

# Function to generate random matrices with the same row and column sums
generate_random_matrix <- function(mat) {
  random_mat <- permatswap(mat, fixedmar = "both", shuffle = "samp", mtype = "count", times = 1)$perm[[1]]
  return(random_mat)
}

# Function to perform co-occurrence analysis on a subset of the data
cooccurrence_analysis <- function(subset_data, num_random_matrices = 100) {
  # Replace NA values with 0s
  subset_data[is.na(subset_data)] <- 0
  
  mat <- as.matrix(subset_data)
  
  # Compute the observed co-occurrence matrix
  observed_cooccurrence <- compute_cooccurrence(mat)
  
  # Initialize a list to store co-occurrence matrices from random matrices
  random_cooccurrences <- vector("list", num_random_matrices)
  
  # Generate random matrices and compute their co-occurrence matrices
  set.seed(123)  # For reproducibility
  for (i in 1:num_random_matrices) {
    random_matrix <- generate_random_matrix(mat)
    random_cooccurrences[[i]] <- compute_cooccurrence(random_matrix)
  }
  
  # Convert the list of random co-occurrence matrices to an array
  random_cooccurrence_array <- array(unlist(random_cooccurrences), 
                                     dim = c(nrow(observed_cooccurrence), 
                                             ncol(observed_cooccurrence), 
                                             num_random_matrices))
  
  # Calculate the mean and standard deviation of random co-occurrences for each pair
  random_cooccurrence_mean <- apply(random_cooccurrence_array, c(1, 2), mean)
  random_cooccurrence_sd <- apply(random_cooccurrence_array, c(1, 2), sd)
  
  # Compute z-scores to see how the observed co-occurrences compare to random expectations
  z_scores <- (observed_cooccurrence - random_cooccurrence_mean) / random_cooccurrence_sd
  
  return(z_scores)
}

# Load your data
data <- read.csv("/mnt/data/vocal_participation.csv")

# Select only the relevant columns (exclude count, leader)
data_selected <- data %>% select(-count, -leader)

# Split the data by flock_ID
split_data <- split(data_selected, data_selected$flock_ID)

# Perform co-occurrence analysis for each subset and store the results
results <- lapply(split_data, function(subset) {
  print(paste("Processing flock ID:", unique(subset$flock_ID)))  # Diagnostic message
  flock_id <- unique(subset$flock_ID)
  species_ids <- subset$species_ID
  subset <- subset %>% select(-flock_ID, -species_ID)  # Remove the flock_ID and species_ID columns
  z_scores <- cooccurrence_analysis(subset, num_random_matrices = 100)
  list(flock_id = flock_id, z_scores = z_scores, species_ids = species_ids)
})

# Convert results into a data frame
results_df <- do.call(rbind, lapply(results, function(result) {
  z_scores <- as.data.frame(result$z_scores)
  colnames(z_scores) <- result$species_ids
  z_scores$flock_id <- result$flock_id
  return(z_scores)
}))

# Print the results data frame
print(results_df)