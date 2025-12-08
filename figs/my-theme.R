my_theme <- function(base_size = 8,
                     base_family = "", 
                     base_line_size = base_size / 20, 
                     base_rect_size = base_size / 20) {
  # Set color and half-line size for margins
  my_black <- "#333333"
  half_line <- base_size / 2
  
  # Base theme customization
  theme_grey(base_size = base_size,
             base_family = base_family, 
             base_line_size = base_line_size,
             base_rect_size = base_rect_size) %+replace%
    theme(
      # Axis settings
      axis.line = element_blank(),
      axis.text = element_text(colour = my_black, size = rel(0.9)),
      axis.title = element_text(colour = my_black, size = rel(1)),
      axis.ticks = element_line(colour = my_black),
      
      # Legend
      legend.key = element_rect(fill = "white", colour = NA),
      legend.text = element_text(size = rel(0.9)),
      
      # Panel
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = my_black),
      panel.grid = element_blank(),
      
      # Plot settings
      plot.margin = margin(10, 10, 10, 10),
      plot.title = element_text(
        size = rel(1.1), 
        hjust = 0, 
        vjust = 1, 
        margin = margin(b = half_line)
      ),
      
      # Strip settings (e.g., facet labels)
      strip.background = element_rect(
        fill = my_black, 
        colour = my_black, 
        linewidth = base_line_size / 2
      ),
      strip.text = element_text(
        colour = "white", 
        size = rel(0.9), 
        margin = margin(0.8 * half_line)
      ),
      
      # Text color
      text = element_text(colour = my_black),
      
      # Mark as complete theme
      complete = TRUE
    )
}
