# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(scales)
library(forcats)
library(readr)

# Set working directory
setwd("path/to/your/working/directory")

########################################
# Data Loading and Preparation
########################################

## Load and prepare steroid data
steroid_data <- read.csv("path/to/data/steroids.csv",
                         row.names = 1, header = TRUE, check.names = FALSE)
# Transpose so that steroid names become column names
steroid_data <- t(steroid_data)
# Convert rownames to a column and ensure proper data types
steroid_data <- data.frame(patient = rownames(steroid_data), steroid_data,
                           row.names = NULL, check.names = FALSE)
numeric_cols <- setdiff(colnames(steroid_data), "patient")
steroid_data[numeric_cols] <- sapply(steroid_data[numeric_cols], as.numeric)

## Load and prepare RNA-seq data
RNA_seq_data <- read.csv("path/to/data/row_counts_norm.csv",
                         row.names = 1, header = TRUE, check.names = FALSE)
RNA_seq_data <- t(RNA_seq_data)
RNA_seq_data <- data.frame(patient = rownames(RNA_seq_data), RNA_seq_data,
                           row.names = NULL, check.names = FALSE)
numeric_cols <- setdiff(colnames(RNA_seq_data), "patient")
RNA_seq_data[numeric_cols] <- sapply(RNA_seq_data[numeric_cols], as.numeric)

## Load gene interaction information
gene_info <- read.csv("path/to/data/interaction_input_SteroidsDB.csv",
                      header = TRUE, check.names = FALSE)
rownames(gene_info) <- gene_info$steroid
# Remove entries with missing gene info
gene_info <- gene_info %>% filter(producing_gene != "" & receptor_gene != "")

########################################
# Function Definitions
########################################

# Calculate the average sum for the specified columns (signals) across all patients.
calc_avg_of_sums <- function(data, col_names) {
  # Trim any extraneous whitespace from the column names
  col_names <- str_trim(col_names)
  
  # Identify matching and non-matching columns
  matching_cols <- intersect(col_names, colnames(data))
  non_matching_cols <- setdiff(col_names, colnames(data))
  
  if (length(non_matching_cols) > 0) {
    message("These column names do not match any columns in the data: ", 
            toString(non_matching_cols))
  }
  
  if (length(matching_cols) > 0) {
    # Sum the specified columns for each patient
    patient_sums <- rowSums(data[, matching_cols, drop = FALSE], na.rm = TRUE)
    # Return the average sum across patients
    avg_sum <- mean(patient_sums, na.rm = TRUE)
    return(avg_sum)
  } else {
    message("No matching columns found for: ", toString(col_names))
    return(NA)
  }
}

########################################
# Data Analysis
########################################

# For each steroid in gene_info, calculate:
# 1. The average steroid concentration.
# 2. The average producing gene expression.
# 3. The average receptor gene expression.
results <- lapply(rownames(gene_info), function(steroid) {
  steroid_class <- gene_info[steroid, "class"]
  producing_genes <- unlist(strsplit(gene_info[steroid, "producing_gene"], ","))
  receptor_genes <- unlist(strsplit(gene_info[steroid, "receptor_gene"], ","))
  
  message("Processing steroid: ", steroid)
  
  avg_steroid_conc <- calc_avg_of_sums(steroid_data, steroid)
  message("Average steroid concentration for ", steroid, ": ", avg_steroid_conc)
  
  avg_producing_expr <- calc_avg_of_sums(RNA_seq_data, producing_genes)
  message("Average producing gene expression for ", steroid, ": ", avg_producing_expr)
  
  avg_receptor_expr <- calc_avg_of_sums(RNA_seq_data, receptor_genes)
  message("Average receptor gene expression for ", steroid, ": ", avg_receptor_expr)
  
  c(steroid, steroid_class, avg_steroid_conc, avg_producing_expr, avg_receptor_expr)
})

# Convert the list of results to a dataframe and set appropriate column names
df <- as.data.frame(do.call(rbind, results), stringsAsFactors = FALSE)
colnames(df) <- c("steroid", "class", "avg_steroid_conc", "avg_producing_expr", "avg_receptor_expr")

# Convert the signal columns to numeric
df <- df %>%
  mutate(across(c(avg_steroid_conc, avg_producing_expr, avg_receptor_expr), as.numeric))

# Standardize (z-score) the signals and compute the overall score
df <- df %>%
  mutate(across(c(avg_steroid_conc, avg_producing_expr, avg_receptor_expr), scale)) %>%
  mutate(overall_score = rowMeans(select(., avg_steroid_conc, avg_producing_expr, avg_receptor_expr), 
                                  na.rm = TRUE))

# Adjust overall scores to ensure they are positive (shift by the absolute minimum)
min_overall_score <- min(df$overall_score, na.rm = TRUE)
df <- df %>% mutate(overall_score = overall_score + abs(min_overall_score))

# Order steroids by class and overall score.
# (Note: 'class' is converted to a factor with a specified order.)
df <- df %>%
  mutate(class = factor(class, levels = c("sterols", "androgen", "Progestogen",
                                          "cortical hormon", "estrogen", "Vitamin"))) %>%
  arrange(class, overall_score) %>%
  mutate(steroid_order = as.numeric(factor(steroid, levels = steroid, ordered = TRUE)))

# Calculate the average overall score by steroid class
avg_scores <- df %>%
  group_by(class) %>%
  summarise(avg_score = mean(overall_score, na.rm = TRUE), .groups = "drop")

# Split the data into two subsets for plotting
df_subset1 <- df %>% filter(class %in% c("sterols", "Vitamin"))
df_subset2 <- df %>% filter(!class %in% c("sterols", "Vitamin"))

avg_scores_subset1 <- avg_scores %>% filter(class %in% c("sterols", "Vitamin"))
avg_scores_subset2 <- avg_scores %>% filter(!class %in% c("sterols", "Vitamin"))

########################################
# Plotting
########################################

# Function to create a bar plot of overall scores with class-specific average lines.
create_plot <- function(data, avg_data) {
  ggplot(data, aes(x = reorder(steroid, steroid_order), y = overall_score, fill = class)) +
    geom_bar(stat = "identity", color = "black", size = 0.25) +
    geom_hline(data = avg_data, aes(yintercept = avg_score, color = class),
               linetype = "dashed", size = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1", guide = "none") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.minor = element_blank(),
          legend.position = "bottom") +
    labs(x = "Steroid", y = "Overall Score",
         title = "Steroid Signals in TNBC Patients", fill = "Class")
}

# Generate plots for the two subsets
plot_subset1 <- create_plot(df_subset1, avg_scores_subset1)
plot_subset2 <- create_plot(df_subset2, avg_scores_subset2)

# Save the plots to a PDF file
pdf("Steroids_signals_in_TNBC.pdf", width = 6, height = 6)
print(plot_subset1)
print(plot_subset2)
dev.off()