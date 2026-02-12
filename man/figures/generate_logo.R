#!/usr/bin/env Rscript

# parseMetab Logo Generator
# Creates a professional logo for the parseMetab package

library(ggplot2)
library(grid)

# Create logo
pdf(file = "parseMetab_logo.pdf", width = 8, height = 8)

# Create the plot
p <- ggplot() +
  # Background circle (light blue metabolic theme)
  annotate("circle", x = 0.5, y = 0.5, r = 0.45, fill = "#E8F4F8", color = "#2C5AA0", size = 2) +
  
  # Metabolic pathway network nodes (representing KEGG pathways)
  annotate("point", 
           x = c(0.3, 0.7, 0.5, 0.35, 0.65), 
           y = c(0.7, 0.7, 0.3, 0.4, 0.4), 
           size = c(10, 10, 14, 8, 8), 
           color = "#FF6B6B", 
           alpha = 0.85) +
  
  # Connecting lines (representing metabolic pathways)
  annotate("segment", 
           x = c(0.3, 0.7, 0.3, 0.7, 0.5), 
           y = c(0.7, 0.7, 0.7, 0.7, 0.3),
           xend = c(0.5, 0.5, 0.35, 0.65, 0.35), 
           yend = c(0.3, 0.3, 0.4, 0.4, 0.4),
           color = "#4ECDC4", 
           size = 1.5, 
           alpha = 0.6) +
  
  # DNA/Omics elements on sides
  annotate("segment", x = 0.15, y = 0.35, xend = 0.15, yend = 0.65, 
           color = "#95E1D3", size = 3, alpha = 0.7) +
  annotate("segment", x = 0.85, y = 0.35, xend = 0.85, yend = 0.65, 
           color = "#95E1D3", size = 3, alpha = 0.7) +
  
  # Analysis elements (representing data analysis)
  annotate("polygon", 
           x = c(0.5, 0.45, 0.55), 
           y = c(0.08, 0.0, 0.0),
           fill = "#A8E6CF", 
           color = "#2C5AA0", 
           alpha = 0.85, 
           size = 1) +
  
  # Package name
  annotate("text", x = 0.5, y = -0.12, 
           label = "parseMetab", 
           size = 10, 
           fontface = "bold", 
           color = "#2C5AA0") +
  
  # Subtitle
  annotate("text", x = 0.5, y = -0.18, 
           label = "Metabolic Activity Analysis", 
           size = 5, 
           color = "#4ECDC4") +
  
  # Version
  annotate("text", x = 0.5, y = -0.23, 
           label = "v1.0.0", 
           size = 4, 
           color = "#666666") +
  
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  coord_fixed(xlim = c(-0.05, 1.05), ylim = c(-0.28, 1.05))

print(p)

dev.off()

# Also create PNG version for web
png(file = "parseMetab_logo.png", width = 800, height = 800, bg = "white")
print(p)
dev.off()

cat("\n✓ Logo saved as parseMetab_logo.pdf\n")
cat("✓ Logo also saved as parseMetab_logo.png (for web use)\n\n")
