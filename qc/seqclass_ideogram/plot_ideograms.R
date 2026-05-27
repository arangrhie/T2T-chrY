#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters
ISSUES_TRACK_HEIGHT <- 1.3  # mm
SEQCLASS_TRACK_HEIGHT <- 2  # mm
SPACE_BETWEEN_SAMPLES_MM <- 1
MAX_RULER_SIZE <- 90e6  # 90 Mb
RULER_TICK_INTERVAL <- 5e6  # 5 Mbps
RULER_HEIGHT_MM <- 2
TRACK_WIDTH_MM <- 180
SAMPLE_LABEL_OFFSET_MM <- 2  # Distance from tracks
TEXT_FONT_SIZE <- 6  # pt

# Get all unique sample names from bed files
bed_dir <- "chr1"
sample_order <- read.table("sample_order.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
all_files <- list.files(bed_dir, pattern = "\\.main\\.bed$")
available_samples <- gsub("\\.issues\\.main\\.bed$|\\.seqclass\\.main\\.bed$", "", all_files) %>% unique()
sample_names <- sample_order$ID[sample_order$ID %in% available_samples]

# Function to read a BED file
read_bed <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  df <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # Standard BED format: chrom, chromStart, chromEnd, name, score, strand, ..., itemRgb (column 9)
  if (ncol(df) >= 3) {
    colnames(df)[1:3] <- c("chrom", "start", "end")
    if (ncol(df) >= 4) colnames(df)[4] <- "name"
    if (ncol(df) >= 5) colnames(df)[5] <- "score"
    if (ncol(df) >= 6) colnames(df)[6] <- "strand"
    if (ncol(df) >= 9) {
      colnames(df)[9] <- "itemRgb"
    } else {
      df$itemRgb <- "100,100,100"  # Default gray if not present
    }
  }
  
  # Convert RGB codes to hex colors
  df$color <- sapply(df$itemRgb, function(rgb_str) {
    rgb_vals <- as.numeric(unlist(strsplit(rgb_str, ",")))
    if (length(rgb_vals) == 3) {
      rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], maxColorValue = 255)
    } else {
      "#646464"  # Default gray
    }
  })
  
  return(df %>%
    select(start, end, color, everything()) %>%
    filter(end <= MAX_RULER_SIZE))  # Filter to ruler size
}

# Prepare data for plotting
plot_data <- data.frame()
y_position <- 0
sample_labels <- data.frame()

for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  
  # Read both files
  issues_file <- file.path(bed_dir, paste0(sample_name, ".issues.main.bed"))
  seqclass_file <- file.path(bed_dir, paste0(sample_name, ".seqclass.main.bed"))
  
  issues_df <- read_bed(issues_file)
  seqclass_df <- read_bed(seqclass_file)
  
  # Process issues track (top)
  if (!is.null(issues_df)) {
    issues_df$track_type <- "issues"
    issues_df$track_y <- y_position + SEQCLASS_TRACK_HEIGHT + ISSUES_TRACK_HEIGHT / 2
    issues_df$sample <- sample_name
    issues_df$sample_order <- i
    plot_data <- rbind(plot_data, issues_df %>%
      select(start, end, color, track_type, track_y, sample, sample_order))
  }
  
  # Process seqclass track (bottom)
  if (!is.null(seqclass_df)) {
    seqclass_df$track_type <- "seqclass"
    seqclass_df$track_y <- y_position + SEQCLASS_TRACK_HEIGHT / 2
    seqclass_df$sample <- sample_name
    seqclass_df$sample_order <- i
    plot_data <- rbind(plot_data, seqclass_df %>%
      select(start, end, color, track_type, track_y, sample, sample_order))
  }
  
  # Store sample label position (SAMPLE_LABEL_OFFSET_MM physical mm to the left of tracks)
  label_x <- -(SAMPLE_LABEL_OFFSET_MM * (MAX_RULER_SIZE / TRACK_WIDTH_MM))
  sample_labels <- rbind(sample_labels, data.frame(
    x = label_x,
    y = y_position + (SEQCLASS_TRACK_HEIGHT + ISSUES_TRACK_HEIGHT) / 2,
    label = sample_name,
    stringsAsFactors = FALSE
  ))
  
  # Move to next sample position
  y_position <- y_position + SEQCLASS_TRACK_HEIGHT + ISSUES_TRACK_HEIGHT + SPACE_BETWEEN_SAMPLES_MM
}

# Create ruler
ruler_ticks <- seq(0, MAX_RULER_SIZE, by = RULER_TICK_INTERVAL)
ruler_labels <- paste0(ruler_ticks / 1e6)

ruler_data <- data.frame(
  x = ruler_ticks,
  y_start = y_position,
  y_end = y_position + RULER_HEIGHT_MM
)

# Create the plot
# Prepare track heights for each row
plot_data <- plot_data %>%
  mutate(height = ifelse(track_type == "issues", ISSUES_TRACK_HEIGHT, SEQCLASS_TRACK_HEIGHT))

p <- ggplot(plot_data, aes(xmin = start, xmax = end, ymin = track_y - height/2, ymax = track_y + height/2, fill = color)) +
  # Plot the BED regions as rectangles
  geom_rect(color = NA) +
  scale_fill_identity() +
  
  # Add ruler ticks
  geom_segment(data = ruler_data, aes(x = x, xend = x, y = y_start, yend = y_end),
               inherit.aes = FALSE, color = "black", linewidth = 0.25) +
  
  # Add ruler baseline
  geom_segment(data = data.frame(x = 0, xend = MAX_RULER_SIZE, 
                                  y = y_position, yend = y_position),
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, color = "black", linewidth = 0.25) +
  
  # Add sample labels
  geom_text(data = sample_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1, vjust = 0.5, size = TEXT_FONT_SIZE / .pt, angle = 0) +
  
  # Add ruler labels
  geom_text(data = data.frame(x = ruler_ticks, y = y_position + RULER_HEIGHT_MM + 1, label = ruler_labels),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = TEXT_FONT_SIZE / .pt, angle = 0) +
  
  # Scales and limits
  scale_x_continuous(limits = c(-15e6, MAX_RULER_SIZE), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, y_position + RULER_HEIGHT_MM + 8), expand = c(0, 0)) +
  
  # Theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(1.5, 0.5, 0.5, 0.5), "cm")
  )

# Output plot
output_file <- "seqclass_ideogram.pdf"
pdf(output_file, width = 12, height = (y_position + RULER_HEIGHT_MM + 20) / 25.4)
print(p)
dev.off()

cat("Plot saved to:", output_file, "\n")

output_file <- "seqclass_ideogram.png"
png(output_file, width = 12, height = (y_position + RULER_HEIGHT_MM + 20) / 25.4, units = "in", res = 300)
print(p)
dev.off()

cat("Plot saved to:", output_file, "\n")
