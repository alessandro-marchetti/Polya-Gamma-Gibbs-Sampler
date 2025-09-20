rm(list = ls())

# 1. Load necessary libraries
library(ggplot2)
library(dplyr)      # For bind_rows and data manipulation
library(tidyr)      # For pivot_longer

# 2. Define Settings and File Paths
# List of methods (as used in your file names)
methods_code <- c("gs", "stan", "bpr0") # Assuming "bpr0" is your BPR method

# Corresponding labels for plotting
plot_method_labels <- c("GS", "BPR0", "Stan") # Renamed BPR0 to BPR for clarity

# Values for 'p' (number of covariates)
p_vals <- c(5, 10, 20)

# Values for 'n' (sample size)
n_vals <- c(25, 50, 100, 200)

# Convert method code (from file names) to plot labels
method_labels_map <- c(gs = "GS", stan = "Stan", bpr0 = "BPR0")

# Where your .RData files are stored
# IMPORTANT: Change this to your actual folder path!
base_dir <- "C:/Users/almar/OneDrive/Desktop/Data simulations"
setwd(base_dir)

# 3. Data Loading and Processing Loop for Densities
all_density_data <- list() # Initialize an empty list to store density data

for (method_idx in seq_along(methods_code)) {
  current_method_code <- methods_code[method_idx]
  current_method_label <- method_labels_map[[current_method_code]]
  
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
      
      if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        warning(paste("Loaded object for", file_name, "is not a valid matrix or is empty. Skipping."))
        next
      }
      
      # For density plots, we need to gather all samples for all parameters
      # Let's consider the first 6 parameters for plotting, if available.
      # Or, if p_val is small, plot all.
      num_params_to_plot <- min(ncol(mat), 6) # Plot up to 6 parameters for density
      
      for (param_idx in 1:num_params_to_plot) {
        # Estimate density for each parameter
        # Using density() function for a kernel density estimate
        dens <- density(mat[, param_idx])
        temp_df <- data.frame(
          x = dens$x,
          y = dens$y,
          method = current_method_label,
          num_covariates = p_val,
          n_sample = n_val,
          parameter = paste0("beta_", param_idx) # Label for each parameter
        )
        all_density_data[[length(all_density_data) + 1]] <- temp_df
      }
    }
  }
}

# Combine all temporary data frames into one master data frame for densities
density_data <- bind_rows(all_density_data)

# 4. Factor Conversion (Crucial for correct plotting order and grouping)
density_data$method <- factor(density_data$method, levels = plot_method_labels)
density_data$num_covariates <- factor(density_data$num_covariates)
density_data$n_sample <- factor(density_data$n_sample)
density_data$parameter <- factor(density_data$parameter)

# 5. Create the Plot using ggplot2 for Densities
windows() # Uncomment if you want the plot in a separate window

plot_densities <- ggplot(density_data, aes(x = x, y = y, color = method, linetype = method)) +
  geom_line(size = 1) +
  labs(
    x = "Parameter Value",
    y = "Density",
    color = "Algorithm",
    linetype = "Algorithm"
  ) +
  # Manually set colors and line types
  scale_color_manual(values = c("GS" = "blue", "BPR0" = "red", "Stan" = "yellow")) +
  scale_linetype_manual(values = c("GS" = "dashed", "BPR0" = "dashed", "Stan" = "solid")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.5
    ),
 #   axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )

# Conditional faceting based on p_val and n_sample
plot_list <- list()

for (p_val in p_vals) {
  for (n_val in n_vals) {
    current_density_data <- density_data %>%
      filter(num_covariates == p_val, n_sample == n_val)
    
    # Load EP approximation object
    ep_var_name <- paste0("EP_p", p_val, "_n", n_val, ".RData")
    load(ep_var_name)
    
    mu_vec <- EP_par[[1]]
    Sigma_mat <- EP_par[[2]]
    
    # Build EP Gaussian curves
    ep_density_list <- list()
    num_params_to_plot <- min(length(mu_vec), 6)
    
    for (j in 1:num_params_to_plot) {
      mu_j <- mu_vec[j]
      sd_j <- Sigma_mat[j]
      x_range <- seq(mu_j - 4 * sd_j, mu_j + 4 * sd_j, length.out = 500)
      y_vals <- dnorm(x_range, mean = mu_j, sd = sd_j)
      ep_df <- data.frame(
        x = x_range,
        y = y_vals,
        method = "EP",
        num_covariates = factor(p_val),
        n_sample = factor(n_val),
        parameter = paste0("beta_", j)
      )
      ep_density_list[[j]] <- ep_df
    }
    
    ep_density_df <- bind_rows(ep_density_list)
    
    # Combine with existing densities
    full_density_data <- bind_rows(current_density_data, ep_density_df)
    
    # Ensure 'method' is a factor with correct level ordering
    full_density_data$method <- factor(full_density_data$method,
                                       levels = c("GS", "BPR0", "Stan", "EP"))
    
    # Plot title
    #  plot_title <- paste0("Densities for p = ", p_val, ", n = ", n_val)
    
    # Plot
    current_plot <- ggplot(full_density_data %>% 
                             filter(as.numeric(sub("beta_", "", parameter)) <= 6),
                           aes(x = x, y = y, color = method, linetype = method)) +
      geom_line(size = 1.1) +
      labs(
        #        title = plot_title,
        x = "Parameter Value",
        y = "Density",
        color = "Method",
        linetype = "Method"
      ) +
      scale_color_manual(values = c(
        "GS" = "blue",
        "BPR0" = "red",
        "Stan" = "yellow",
        "EP" = "orange"
      )) +
      scale_linetype_manual(values = c(
        "GS" = "dashed",
        "BPR0" = "dashed",
        "Stan" = "solid",
        "EP" = "solid"
      )) +
      facet_wrap(~ parameter, scales = "free", ncol = 2) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom")
    
    # Save plot to list
    plot_list[[paste0("p", p_val, "n", n_val)]] <- current_plot
  }
}


windows()
# Print all generated plots
for (plot_name in names(plot_list)) {
  print(plot_list[[plot_name]])
}
print(plot_list[["p20n200"]])

# Optional: Save the plots to files
# For example, to save the plot for p=5, n=25:
# ggsave("density_plot_p5_n25.png", plot = plot_list[["p5n25"]], width = 10, height = 8, dpi = 300)