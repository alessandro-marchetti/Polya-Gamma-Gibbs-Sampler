rm(list = ls())

# In this code the label "bpr0" stands for Canale-D'Angelo's Metropolis Hastings algorithm (bpr)
# when its chains are started at the zero vector

library(ggplot2)
library(dplyr)      # For data manipulation
library(coda)       # For effectiveSize

methods_code <- c("gs", "stan", "bpr0")

plot_method_labels <- c("GS", "BPR0", "Stan")

# Values for p
p_vals <- c(5, 10, 20)

# Values for n
n_vals <- c(25, 50, 100, 200)

method_labels_map <- c(gs = "GS", stan = "Stan", bpr0 = "BPR0")
target_iterations <- c(
  "gs" = 3000,
  "bpr0" = 3000,
  "stan" = 5000
)

# Set the working directory
base_dir <- "C:/Users/almar/OneDrive/Desktop/Data simulations"
setwd(base_dir)

# Data Loading and computation of ESS/iterations
all_data <- list()

for (method_idx in seq_along(methods_code)) {
  current_method_code <- methods_code[method_idx]
  current_method_label <- method_labels_map[[current_method_code]]
  current_target_iter <- target_iterations[[current_method_code]]
  
  if (is.null(current_target_iter) || current_target_iter <= 0) {
    warning(paste("Invalid or missing target_iterations for method:", current_method_code, ". Skipping this method."))
    next
  }
  
  for (p_val in p_vals) {
    for (n_val in n_vals) {
      file_name <- paste0(current_method_code, "_p", p_val, "_n", n_val, ".RData")
      var_name <- paste0("samples_", current_method_code, "_p", p_val, "_n", n_val)
      
      if (!file.exists(file_name)) {
        warning(paste("File not found:", file_name, ". Skipping this combination."))
        next
      }
      
      load(file_name) # This loads the matrix into the workspace
      mat <- get(var_name) # Retrieve the matrix by its name
      
      # Ensure matrix is not empty and has enough iterations
      if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        warning(paste("Loaded object for", file_name, "is not a valid matrix or is empty. Skipping."))
        next
      }
      if (nrow(mat) < 100) { # A threshold for reliable ESS calculation
        warning(paste("Matrix for", file_name, "has too few iterations (", nrow(mat), ") for reliable ESS calculation. Proceeding but results may be unreliable."))
      }
      
      ess_vals <- effectiveSize(mat) # This returns a vector of ESS values (one per parameter 'p')
      
      # Calculate log(ESS / fixed_iterations) for each parameter
      log_ess_per_iter <- log(ess_vals / current_target_iter)
      
      # Create a temporary data frame for this specific combination
      temp_df <- data.frame(
        log_ESS_per_iter = log_ess_per_iter,
        method = current_method_label,
        num_covariates = p_val,
        n_sample = n_val
      )
      
      all_data[[length(all_data) + 1]] <- temp_df
    }
  }
}

# Combine all temporary data frames into one master data frame
your_data <- bind_rows(all_data)

# Ensure methods are factors with a specific order for legend and plot appearance
your_data$method <- factor(your_data$method, levels = plot_method_labels)
your_data$num_covariates <- factor(your_data$num_covariates)
your_data$n_sample <- factor(your_data$n_sample)

windows()

# Define position dodging for possible overlapping boxplots
pd <- position_dodge(width = 0.8) # set width = 0 to get centered boxplots

plot_final <- ggplot(your_data, aes(x = num_covariates, y = log_ESS_per_iter, color = method, fill = method)) +
  geom_boxplot(alpha = 0.7, position = pd, outlier.shape = 1, outlier.size = 1.5) +
  stat_summary(fun = median, geom = "line", aes(group = method), size = 1, position = pd, linetype = "dashed") +
  stat_summary(fun = median, geom = "point", aes(group = method), size = 3, position = pd) +
  
  facet_wrap(~ n_sample, scales = "free_y",
             labeller = labeller(n_sample = function(x) paste0("n = ", x))) +
  
  labs(
    x = "Number of covariates",
    y = expression(log(ESS / Iterations)),
    color = "Algorithm",
    fill = "Algorithm"
  ) +
  
  # Manually set colors for methods, including BPR0
  scale_color_manual(values = c("GS" = "blue", "BPR0" = "red", "Stan" = "yellow")) +
  scale_fill_manual(values = c("GS" = "blue", "BPR0" = "red", "Stan" = "yellow")) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.5
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )

# Print the final plot
print(plot_final)

