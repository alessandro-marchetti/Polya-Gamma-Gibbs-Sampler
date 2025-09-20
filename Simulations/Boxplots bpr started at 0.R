rm(list = ls())

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define method codes and labels
methods_code <- c("gs", "stan", "bpr0")
plot_method_labels <- c("GS", "BPR0", "Stan")
method_labels_map <- c(gs = "GS", stan = "Stan", bpr0 = "BPR0")

# Parameter and sample sizes
p_vals <- c(5, 10, 20)
n_vals <- c(25, 50, 100, 200)

# Set working directory
base_dir <- "C:/Users/almar/OneDrive/Desktop/Data simulations"
setwd(base_dir)

# List to collect everything
all_box_data <- list()

# Loop to load and reshape data
for (method_idx in seq_along(methods_code)) {
  current_method_code <- methods_code[method_idx]
  current_method_label <- method_labels_map[[current_method_code]]
  
  for (p_val in p_vals) {
    for (n_val in n_vals) {
      file_name <- paste0(current_method_code, "_p", p_val, "_n", n_val, ".RData")
      var_name <- paste0("samples_", current_method_code, "_p", p_val, "_n", n_val)
      
      if (!file.exists(file_name)) {
        warning(paste("File not found:", file_name))
        next
      }
      
      load(file_name)
      mat <- get(var_name)
      
      if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        warning(paste("Invalid matrix in:", file_name))
        next
      }
      
      num_params_to_plot <- min(ncol(mat), 6)
      temp_df <- as.data.frame(mat[, 1:num_params_to_plot])
      colnames(temp_df) <- paste0("beta_", 1:num_params_to_plot)
      
      temp_long <- pivot_longer(temp_df, cols = everything(), names_to = "parameter", values_to = "value")
      temp_long$method <- current_method_label
      temp_long$num_covariates <- factor(p_val)
      temp_long$n_sample <- factor(n_val)
      
      all_box_data[[length(all_box_data) + 1]] <- temp_long
    }
  }
}

# Combine all into one data frame
box_data <- bind_rows(all_box_data)

# Set factor levels
box_data$method <- factor(box_data$method, levels = plot_method_labels)
box_data$parameter <- factor(box_data$parameter)
box_data$num_covariates <- factor(box_data$num_covariates)
box_data$n_sample <- factor(box_data$n_sample)

# Color map
color_map <- c("GS" = "blue", "BPR0" = "red", "Stan" = "yellow")

# Create plots
plot_list <- list()

for (p_val in p_vals) {
  for (n_val in n_vals) {
    current_data <- box_data %>%
      filter(num_covariates == p_val, n_sample == n_val)
    
#    plot_title <- paste0("Boxplots for p = ", p_val, ", n = ", n_val)
    
    current_plot <- ggplot(current_data, aes(x = method, y = value, fill = method)) +
      geom_boxplot(outlier.shape = 1, alpha = 0.8) +
      facet_wrap(~ parameter, scales = "free", ncol = 2) +
      labs(
#        title = plot_title,
#        x = "Method",
        y = "Parameter Value"
      ) +
      scale_fill_manual(values = color_map) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 0)
      )
    
    plot_list[[paste0("p", p_val, "_n", n_val)]] <- current_plot
  }
}
# Show example
windows()
print(plot_list[["p20_n200"]])

# Optional: save
# ggsave("boxplot_p5_n25.png", plot = plot_list[["p5_n25"]], width = 10, height = 8, dpi = 300)
