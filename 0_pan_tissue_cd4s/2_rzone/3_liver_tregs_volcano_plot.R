library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)

initialize_color_vector <- function(df, default_color = 'white') {
  return(rep(default_color, nrow(df)))
}

assign_colors <- function(df, color_vector) {
  color_vector[df$'liver_l' > 0 &
                 df$'liver_p' < 10e-6 &
                 df$'liver_l' > 0.5] <- 'darkgreen'
  
  color_vector[df$'liver_l' < 0 &
                 df$'liver_p' < 10e-6 &
                 df$'liver_l' < -0.5] <- 'darkgray'
  
  return(color_vector)
}

set_color_names <- function(df, color_vector) {
  names(color_vector) <- ifelse(df$'liver_l' > 0 &
                                  df$'liver_p' < 10e-6 &
                                  df$'liver_l' > 0.5,
                                'Liver',
                                ifelse(df$'liver_l' < 0 &
                                         df$'liver_p' < 10e-6 &
                                         df$'liver_l' < -0.5,
                                       'Rest',
                                       'Not Significant'))
  return(color_vector)
}

plot_volcano <- function(df, color_vector) {
  EnhancedVolcano(
    df,
    title = "Liver vs. Rest",
    subtitle = "",
    lab = df$'liver_n',
    x = 'liver_l',
    y = 'liver_p',
    xlim = c(-10, 10),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-6,
    FCcutoff = 0.585,
    pointSize = 3.0,
    labSize = 5.0,
    colAlpha = 1,
    legendLabels = c('Not sig.', 'Log (base 2) FC', 'p-value', 'p-value & Log (base 2) FC'),
    legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize = 2.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    colCustom = color_vector,
    max.overlaps = 50,
    selectLab = c('Srgn', 'Tigit', 'Gzmb', 'S100a6', 'Ly6e', 'Glrx', 'Nkg7', 
                  'Treml2', 'Mgmt', 'Rragd', 'Aff3', 'Tff1', 'Lta', 'Bach2', 'Pde11a'
    )
  )
}

all <- read_excel("~/github/2_tff1/1_pyzone/1_outputs/1_degs/8_liver_treg_vs_rest.xlsx")

# subset to significant only
all <- all[all$`Is Significant` == TRUE, , drop = FALSE]

keyvals <- initialize_color_vector(all)
keyvals <- assign_colors(all, keyvals)
keyvals <- set_color_names(all, keyvals)

all_plt <- plot_volcano(all, keyvals)

# ---- FIXED LEGEND ----
all_plt <- all_plt +
  scale_colour_manual(
    values = c(Liver = "darkgreen", Rest = "darkgray"),
    breaks = c("Liver", "Rest"),
    name = NULL
  ) +
  theme(legend.position = "right")

all_plt <- all_plt +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

out_path_pdf <- "~/github/2_tff1/1_pyzone/1_outputs/0_figures/volcanoplot_liver_vs_rest.pdf"

ggsave(
  filename = out_path_pdf,
  plot = all_plt,
  device = cairo_pdf,
  width = 8, height = 7, units = "in",
  bg = "transparent"
)

all_plt

