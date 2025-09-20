rm(list = ls())

# In this code the label "bpr" stands for Canale-D'Angelo's Metropolis Hastings algorithm (bpr)
# when its chains are started at the default vector, with the maximum likelihood estimates

library(ggplot2)
library(dplyr)      # For data manipulation
library(coda)       # For effectiveSize

methods_code <- c("gs", "stan", "bpr")

plot_method_labels <- c("GS", "BPR", "Stan") # Or c("IS-Polya", "MH-Polya", "Stan") if you prefer these names

# Values for p
p_vals <- c(5, 10, 20)

# Values for n
n_vals <- c(25, 50, 100, 200)

method_labels_map <- c(gs = "GS", stan = "Stan", bpr = "BPR")
target_iterations <- c(
  "gs" = 3000,
  "bpr" = 3000,
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
  
  for (p_val in p_vals) {
    for (n_val in n_vals) {
      file_name <- paste0(current_method_code, "_p", p_val, "_n", n_val, ".RData")
      var_name <- paste0("samples_", current_method_code, "_p", p_val, "_n", n_val)
      
      load(file_name) # This loads the matrix into the workspace
      mat <- get(var_name) # Retrieve the matrix by its name
        
      ess_vals <- effectiveSize(mat) # This returns a vector of ESS values (one per parameter 'p')
        
      # Calculate log(ESS / fixed_iterations) for each parameter
      log_ess_per_iter <- log(ess_vals / current_target_iter)
        
      # Create a temporary data frame for this specific combination
      # This will create 'p_val' rows, each with one log_ess_per_iter value
      temp_df <- data.frame(
        log_ESS_per_iter = log_ess_per_iter,
        method = current_method_label, # Use the plot label
        num_covariates = p_val, # Renamed 'p' to 'num_covariates' for clarity in plot labels
        n_sample = n_val) # Renamed 'n' to 'n_sample' for clarity in plot labels

      all_data[[length(all_data) + 1]] <- temp_df
    }
  }
}

# Combine all temporary data frames into one master data frame
your_data <- bind_rows(all_data)

# 4. Factor Conversion (Crucial for correct plotting order and grouping)
# Ensure methods are factors with a specific order for legend and plot appearance
your_data$method <- factor(your_data$method, levels = plot_method_labels)
# Ensure num_covariates and n_sample are factors for discrete axes/facets
your_data$num_covariates <- factor(your_data$num_covariates)
your_data$n_sample <- factor(your_data$n_sample)


# 5. Create the Plot using ggplot2
# Use windows() or quartz() or x11() if you want the plot in a separate window
windows()

# Define position dodging for overlapping boxplots
pd <- position_dodge(width = 0.8) # set width = 0 to center the boxplots

plot_final <- ggplot(your_data, aes(x = num_covariates, y = log_ESS_per_iter, color = method, fill = method)) +
  # Boxplots (semi-transparent, with custom outlier appearance)
  geom_boxplot(alpha = 0.7, position = pd, outlier.shape = 1, outlier.size = 1.5) +
  # Dashed lines connecting medians of each method across 'num_covariates'
  stat_summary(fun = median, geom = "line", aes(group = method), size = 1, position = pd, linetype = "dashed") +
  # Points for the medians
  stat_summary(fun = median, geom = "point", aes(group = method), size = 3, position = pd) +
  
  # Faceting: Create separate panels for each 'n_sample'
  facet_wrap(~ n_sample, scales = "free_y", # Allow y-axis scale to vary per facet if needed
             labeller = labeller(n_sample = function(x) paste0("n = ", x))) + # Custom facet labels (e.g., "n = 25")
  
  # Labels for axes and legend
  labs(
    x = "Number of covariates",
    y = expression(log(ESS / Iterations)),
    color = "Algorithm", # Legend title for colors
    fill = "Algorithm"   # Legend title for fills
  ) +
  
  # Manually set colors for methods to match your preference/reference (e.g., red, green, blue)
  scale_color_manual(values = c("GS" = "blue", "BPR" = "red", "Stan" = "yellow")) +
  scale_fill_manual(values = c("GS" = "blue", "BPR" = "red", "Stan" = "yellow")) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    # Style for facet labels (text within the strip)
    strip.text = element_text(size = 12, face = "bold"),
    # Style for the background of the facet strips (the "external boxes")
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

