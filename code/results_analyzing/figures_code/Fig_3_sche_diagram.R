#### SEM: Success ~ mND + mRFD ~ mfd + mpd for all species ####
### Original model ###
rm(list = ls())
setwd("D:/R projects/BSS_coexist")

library(ggplot2)
library(patchwork)
library(dplyr)

# =========================
# Colors
# =========================
col_alien   <- "#3F3F3F"  
col_native  <- "#D4D4D4"  
col_overlap <- "#7A7A7A"  
col_gold <- "#9B8325"

# =========================
# Panel function
# =========================
make_panel <- function(mu_alien, mu_native,
                       sd_alien = 0.45,
                       sd_native = 0.45,
                       amp_alien = 1,
                       amp_native = 1,
                       label = "(i)") {
  
  x <- seq(-3, 3, length.out = 2000)
  
  df <- data.frame(
    x = x,
    alien = amp_alien * dnorm(x, mean = mu_alien, sd = sd_alien),
    native = amp_native * dnorm(x, mean = mu_native, sd = sd_native)
  ) %>%
    mutate(overlap = pmin(alien, native))
  
  ggplot(df, aes(x = x)) +
    
    # Alien distribution
    geom_area(
      aes(y = alien),
      fill = col_alien,
      alpha = 0.85,
      colour = "black",
      linewidth = 0.55
    ) +
    
    # Native distribution
    geom_area(
      aes(y = native),
      fill = col_native,
      alpha = 0.80,
      colour = "#555555",
      linewidth = 0.55
    ) +
    
    # Overlap area
    geom_area(
      aes(y = overlap),
      fill = col_overlap,
      alpha = 0.90
    ) +
    
    # Half-frame axes
    geom_segment(
      aes(x = -3, xend = -3, y = 0, yend = 1.05),
      linewidth = 0.75,
      colour = "black"
    ) +
    geom_segment(
      aes(x = -3, xend = 3, y = 0, yend = 0),
      linewidth = 0.75,
      colour = "black"
    ) +
    
    # Panel label
    annotate(
      "text",
      x = -2.8,
      y = 1.03,
      label = label,
      size = 7,
      fontface = "bold",
      hjust = 0
    ) +
    
    # Axis labels
    annotate(
      "text",
      x = 0,
      y = -0.10,
      label = NA,
      size = 5.2
    ) +
    annotate(
      "text",
      x = -3.25,
      y = 0.53,
      label = NA,
      angle = 90,
      size = 4.6
    ) +
    
    coord_fixed(
      ratio = 4.5,
      xlim = c(-3, 3),
      ylim = c(0, 1.1),
      clip = "off"
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(5, 8, 15, 8),
      aspect.ratio = 1
    )
}

# =========================
# Four scenarios
# =========================
# (i)
p1 <- make_panel(
  mu_alien = -0.10,
  mu_native = 0.10,
  amp_alien = 0.90,
  amp_native = 0.90,
  label = "(i)"
)


# (ii)
p2 <- make_panel(
  mu_alien = -0.95,
  mu_native = 0.95,
  amp_alien = 0.90,
  amp_native = 0.90,
  label = "(ii)"
)

# (iii)
p3 <- make_panel(
  mu_alien = 0,
  mu_native = 0,
  amp_alien = 1.10,
  amp_native = 0.48,
  label = "(iii)"
)

# (iv)
p4 <- make_panel(
  mu_alien = -1.20,
  mu_native = 1.00,
  amp_alien = 1.10,
  amp_native = 0.48,
  label = "(iv)"
)

# =========================
# Top title and column labels
# =========================
top_lab <- ggplot() +
  annotate(
    "text",
    x = 0.5,
    y = 0.72,
    label = "Stabilizing niche differences",
    size = 8.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.28,
    y = 0.18,
    label = "Absent",
    size = 6.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.73,
    y = 0.18,
    label = "Present",
    size = 6.5,
    fontface = "bold"
  ) +
  theme_void()

# =========================
# Left title and row labels
# =========================
left_lab <- ggplot() +
  annotate(
    "text",
    x = 0.18,
    y = 0.50,
    label = "Competitive differences",
    angle = 90,
    size = 7.4,
    fontface = "bold",
    colour = col_gold
  ) +
  annotate(
    "text",
    x = 0.72,
    y = 0.76,
    label = "Absent",
    angle = 90,
    size = 6.2,
    fontface = "bold",
    colour = col_gold
  ) +
  annotate(
    "text",
    x = 0.72,
    y = 0.27,
    label = "Present",
    angle = 90,
    size = 6.2,
    fontface = "bold",
    colour = col_gold
  ) +
  theme_void()

# =========================
# Legend
# =========================
legend_plot <- ggplot() +
  
  annotate(
    "rect",
    xmin = 0.06,
    xmax = 0.11,
    ymin = 0.42,
    ymax = 0.62,
    fill = col_alien,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 0.14,
    y = 0.52,
    label = "Alien species",
    hjust = 0,
    size = 5
  ) +
  
  annotate(
    "rect",
    xmin = 0.40,
    xmax = 0.45,
    ymin = 0.42,
    ymax = 0.62,
    fill = col_native,
    colour = "#555555"
  ) +
  annotate(
    "text",
    x = 0.48,
    y = 0.52,
    label = "Native species",
    hjust = 0,
    size = 5
  ) +
  
  annotate(
    "rect",
    xmin = 0.72,
    xmax = 0.77,
    ymin = 0.42,
    ymax = 0.62,
    fill = col_overlap,
    colour = "#555555"
  ) +
  annotate(
    "text",
    x = 0.80,
    y = 0.52,
    label = "Overlap",
    hjust = 0,
    size = 5
  ) +
  
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

# =========================
# Combine figure
# =========================
main_panel <- (p1 + p2) / (p3 + p4)

final_plot <-
  top_lab /
  (left_lab + main_panel + plot_layout(widths = c(0.13, 1))) /
  legend_plot +
  plot_layout(heights = c(0.13, 1, 0.13))

final_plot


# (ii)
p2_pure <- make_panel(
  mu_alien = -0.95,
  mu_native = 0.95,
  amp_alien = 0.90,
  amp_native = 0.90,
  label = NA
)

# (iii)
p3_pure <- make_panel(
  mu_alien = 0,
  mu_native = 0,
  amp_alien = 1.10,
  amp_native = 0.48,
  label = NA
)

# (iv)
p4_pure <- make_panel(
  mu_alien = -1.20,
  mu_native = 1.00,
  amp_alien = 1.10,
  amp_native = 0.48,
  label = NA
)

# =========================
# Save
# =========================

ggsave(
  "results/figures/Fig_3_diag_nd.svg",
  p2_pure,
  width = 9,
  height = 8.5
)

ggsave(
  "results/figures/Fig_3_diag_fd.svg",
  p3_pure,
  width = 9,
  height = 8.5
)

ggsave(
  "results/figures/Fig_3_diag_nfd.svg",
  p4_pure,
  width = 9,
  height = 8.5
)

ggsave(
  "Alien_native_ND_RFD_framework.pdf",
  final_plot,
  width = 9,
  height = 8.5
)



ggsave(
  "Alien_native_ND_RFD_framework.png",
  final_plot,
  width = 9,
  height = 8.5,
  dpi = 600
)