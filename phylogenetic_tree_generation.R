# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read the Excel sheet
args <- commandArgs(trailingOnly = TRUE)

# Check for help option
if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
  cat("Usage: Rscript phylogenetic_tree_generation.R <file_path>\n")
  cat("Description:\n")
  cat("  <file_path>   Path to the Excel file to be processed. Ensure sheet contains 'Genome_Change' and 'ID' columns\n")
  cat("Options:\n")
  cat("  -h, --help    Show this help message and exit\n")
  cat("\n")
  quit(status = 0)  # Exit the script after showing help
}

# Get file path
file_path <- args[1]
data <- read_excel(file_path, sheet = 1)

# Ensure necessary columns are present
if (!all(c("Genome_Change", "ID") %in% colnames(data))) {
  stop("The sheet must contain 'Genome_Change' and 'ID' columns.")
}

# Filter and initialize samples in phylogeny -------------------------------------------

# Count the number of samples each mutation appears in
mutation_counts <- data %>%
  group_by(Genome_Change) %>%
  summarise(count = n_distinct(ID), sample_ID = paste(unique(ID), collapse = ","))

# Identify shared mutations (present in more than one sample)
shared_mutations <- mutation_counts %>%
  filter(count > 1)

# Identify the IDs present in the subsets (shared mutations)
ids_in_subsets <- unique(unlist(strsplit(shared_mutations$sample_ID, ",")))

# Filter mutation_counts to keep only rows where the sample_IDs are part of the subsets or are shared mutations
mutation_counts <- mutation_counts %>%
  filter(grepl(paste(ids_in_subsets, collapse = "|"), sample_ID) | count > 1)

# Function to find overlapping mutations -------------------------------------------

# Function to find overlaps between sets of samples, assign levels sequentially, and sort by num_samples
find_overlapping_mutations <- function(shared_mutations) {
  subsets <- unique(shared_mutations$sample_ID)
  subset_list <- strsplit(subsets, ",")
  
  # Identify the IDs present in the subsets (shared mutations)
  ids_in_subsets <- unique(unlist(strsplit(shared_mutations$sample_ID, ",")))
  
  overlap_list <- data.frame(
    subset = character(0),
    num_samples = integer(0),
    level = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(subset_list)) {
    current_subset <- subset_list[[i]]
    
    # Filter current_subset to keep only IDs present in shared_mutations
    filtered_subset <- current_subset[current_subset %in% ids_in_subsets]
    
    # Only add to overlap_list if the subset is non-empty
    if (length(filtered_subset) > 0) {
      current_subset_str <- paste(filtered_subset, collapse = ",")
      num_samples_in_subset <- length(filtered_subset)
      
      # Add the current subset with its num_samples to the dataframe
      overlap_list <- rbind(overlap_list, data.frame(
        subset = current_subset_str,
        num_samples = num_samples_in_subset,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Check if there is only one subset
  if (nrow(overlap_list) == 1) {
    overlap_list$level <- 0  # Assign level 0 for the single subset
    # Group by Genome_Change and count the number of samples (IDs) each mutation appears in
    mutation_counts <- data %>%
      group_by(Genome_Change) %>%
      summarise(count = n_distinct(ID))
    
    # Count the total number of mutations that are repeated (shared) across samples
    repeated_mutations <- nrow(filter(mutation_counts, count > 1))
    
    # Count the number of unique mutations in each sample
    unique_mutations_per_sample <- data %>%
      group_by(ID) %>%
      summarise(unique_mutations = n_distinct(Genome_Change) - sum(Genome_Change %in% mutation_counts$Genome_Change[mutation_counts$count > 1]))
    
    # Create a data frame for plotting
    plot_data <- unique_mutations_per_sample %>%
      filter(unique_mutations > 0) %>%
      mutate(y_position = seq(1, -1, length.out = n()))  # Alternate positions for top and bottom lines
    
    # Calculate total length of the plot for scaling
    max_x_position <- max(repeated_mutations + unique_mutations_per_sample$unique_mutations)
    
  } else {
    # Sort by num_samples in descending order and assign sequential levels
    overlap_list <- overlap_list %>%
      arrange(desc(num_samples)) %>%
      mutate(level = row_number() - 1)  # Assign sequential levels starting from 0
  }
  
  return(overlap_list)
}

# Data set for overlapping mutations -------------------------------------------

# Get overlapping subsets with levels sorted by num_samples and assigned sequential levels
overlapping_mutations <- find_overlapping_mutations(shared_mutations)

# Add the number of overlapping mutations to the existing DataFrame
overlapping_mutations <- overlapping_mutations %>%
  rowwise() %>%
  mutate(mut_count = sum(
    shared_mutations$sample_ID == subset
  ))

# Prepare the data for plotting
overlapping_mutations <- overlapping_mutations %>%
  arrange(level)

# Initialize x_start and x_end
overlapping_mutations$x_start <- 0
overlapping_mutations$x_end <- overlapping_mutations$mut_count[1]

# Iterate through each level to set the correct x_start and x_end
if (nrow(overlapping_mutations) > 1) {
  for (i in 2:nrow(overlapping_mutations)) {
    overlapping_mutations$x_start[i] <- overlapping_mutations$x_end[i - 1]
    overlapping_mutations$x_end[i] <- overlapping_mutations$x_start[i] + overlapping_mutations$mut_count[i]
  }
}

# Set y_position based on the level
overlapping_mutations <- overlapping_mutations %>%
  mutate(y_position = -level)  # Incremental y spacing for each line, spacing them out vertically

# Data set for unique mutations -------------------------------------------

# Create a dataset for unique mutations (mutations that don't overlap) and sample IDs with no unique mutations
unique_mutations <- mutation_counts %>%
  filter(count == 1) %>%  # Only keep mutations that appear in one sample (non-overlapping)
  group_by(sample_ID) %>%  # Group by sample ID
  summarise(mut_count = n())  # Count the number of unique mutations per sample

# Identify sample IDs that are in the subset but do not have unique mutations
sample_ids_with_no_unique <- setdiff(ids_in_subsets, unique_mutations$sample_ID)

# Add these sample IDs to unique_mutations with a mut_count of 0
unique_mutations <- bind_rows(
  unique_mutations,
  tibble(sample_ID = sample_ids_with_no_unique, mut_count = 0)
)

# Assign levels to unique mutations based on the highest level found in overlapping_mutations
unique_mutations <- unique_mutations %>%
  rowwise() %>%
  mutate(
    level = max(
      overlapping_mutations %>%
        filter(grepl(sample_ID, subset)) %>%
        pull(level),
      na.rm = TRUE
    ) + 1  # Set the level as the highest level found + 1
  )

# Sort and prepare the data for plotting unique mutations
unique_mutations <- unique_mutations %>%
  arrange(level, desc(mut_count==0), sample_ID)%>%
  group_by(level) %>%
  mutate(
    # Match x_start with overlapping_mutations at the same level or use the previous x_end if not found
    x_start = ifelse(
      level %in% overlapping_mutations$level,
      overlapping_mutations$x_start[overlapping_mutations$level == level][1],  # Take the first match if multiple
      ifelse(level > min(overlapping_mutations$level), 
             max(overlapping_mutations$x_end[overlapping_mutations$level < level], na.rm = TRUE),
             0)
    ),
    x_end = x_start + mut_count,  # Calculate x_end based on mut_count
    # Set y_position based on overlapping_mutations or create a spaced sequence for fallback
    y_position = ifelse(
      level %in% overlapping_mutations$level,
      ifelse(
        mut_count == 0,
        overlapping_mutations$y_position[overlapping_mutations$level == level][1] + 5,  # Same y_position for mut_count == 0
        {
          base_y_position <- overlapping_mutations$y_position[overlapping_mutations$level == level][1] + 2
          if (n() > 1) {
            # Spread out y_position for mut_count > 0
            seq(from = base_y_position - 0.5 * (n() - 1), by = 1, length.out = n())
          } else {
            base_y_position
          }
        }
      ),
      {
        # Fallback y_position for items not in overlapping_mutations levels
        prev_y_position <- min(overlapping_mutations$y_position[overlapping_mutations$level < level], na.rm = TRUE)
        if (!is.na(prev_y_position)) {
          seq(from = prev_y_position - 0.5, to = prev_y_position + 0.5, length.out = n())
        } else {
          seq(from = -0.5, to = 0.5, length.out = n())
        }
      }
    )
    
  ) %>%
  ungroup()

# Combined data set for plotting ----------------------------------------------------

# Combine the data for plotting
plot_data <- bind_rows(
  overlapping_mutations %>% mutate(type = "Overlapping"),
  unique_mutations %>% mutate(type = "Unique")
)

# Create data for vertical lines connecting levels
vertical_lines <- plot_data %>%
  arrange(level) %>%
  group_by(level) %>%
  mutate(
    x_position = x_start,  # Position for vertical lines should align with the start of the horizontal lines
    y_start = lag(y_position),  # Start at the y_position of the previous level
    y_end = y_position  # End at the current y_position
  ) %>%
  filter(!is.na(y_start))  # Remove any rows where y_start is NA (e.g., the first level)

# Create a dataset for concatenated labels only
concatenated_labels_data <- unique_mutations %>%
  group_by(level, x_end, y_position) %>%
  summarise(
    concatenated_labels = paste(sample_ID, collapse = "\n"),
    y_position = first(y_position),  # Ensure y_position is consistent for each group
    x_end = first(x_end)             # Keep x_end for placing the labels
  ) %>%
  ungroup()

# Set plotting parameters -----------------------------------------------------

# Plot both overlapping and unique mutations
p <- ggplot() +
  # Add vertical lines
  geom_segment(data = vertical_lines, aes(x = x_position, xend = x_position, y = y_start, yend = y_end), color = "black", linewidth = 0.5) +
  # Add horizontal lines
  geom_segment(data = plot_data, aes(x = x_start, xend = x_end, y = y_position, yend = y_position), linewidth = 0.5) +
  
  # Add "mutations" text only for non-zero mut_count values
  geom_text(
    data = plot_data %>% filter(mut_count > 0),  # Only display labels for mut_count > 0
    aes(x = (x_start + x_end) / 2, y = y_position, label = paste0(mut_count, " mut")),
    hjust = 1.2, vjust = 0.5, color = "black", size = 7, size.unit = "pt"
  ) +
  
  # Add concatenated sample_ID labels separately with independent positioning
  geom_text(
    data = concatenated_labels_data,
    aes(x = x_end, y = y_position, label = concatenated_labels),  # Use concatenated labels
    hjust = 0.5, vjust = -0.5, color = "black", size = 7, size.unit = "pt"
  ) +
  
  # Set the x and y limits to scale the figure properly
  xlim(0, max(unique_mutations$x_end) + 10) +  # Adjust for padding
  ylim(min(unique_mutations$y_position) - 2 , max(unique_mutations$y_position) + 2) +  # Adjust for vertical spacing
  
  theme_minimal() +
  ggtitle("Phylogenetic Tree") +
  xlab("Mutations (mut)") +
  ylab("") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.text.y = element_text(color = "black"),  # Show y-axis text
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  coord_flip()  # Rotate the plot 90 degrees to the left


# Display plot and save ---------------------------------------------------

# Print the plot
print(p)

# Extract the directory from the file_path
output_directory <- dirname(file_path)
# Define the output file path for the PDF
output_file_path <- file.path(output_directory, "phylogenetic_tree.pdf")

# Save the plot as a PDF with adjusted dimensions in the same directory
ggsave(output_file_path, width = 15, height = 15)
