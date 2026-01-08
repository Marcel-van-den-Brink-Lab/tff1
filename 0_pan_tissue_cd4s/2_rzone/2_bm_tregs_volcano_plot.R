library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)

initialize_color_vector <- function(df, default_color = 'white') {
  return(rep(default_color, nrow(df)))
}

assign_colors <- function(df, color_vector) {
  color_vector[df$'bm_l' > 0 &
                 df$'bm_p' < 10e-6 &
                 df$'bm_l' > 0.5] <- 'purple'
  
  color_vector[df$'bm_l' < 0 &
                 df$'bm_p' < 10e-6 &
                 df$'bm_l' < -0.5] <- 'darkgray'
  
  return(color_vector)
}

set_color_names <- function(df, color_vector) {
  names(color_vector) <- ifelse(df$'bm_l' > 0 &
                                  df$'bm_p' < 10e-6 &
                                  df$'bm_l' > 0.5,
                                'BM',
                                ifelse(df$'bm_l' < 0 &
                                         df$'bm_p' < 10e-6 &
                                         df$'bm_l' < -0.5,
                                       'Rest',
                                       'Not Significant'))
  return(color_vector)
}

plot_volcano <- function(df, color_vector) {
  EnhancedVolcano(
    df,
    title = "BM vs. Rest",
    subtitle = "",
    lab = df$'bm_n',
    x = 'bm_l',
    y = 'bm_p',
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
    selectLab = c('S100a10', 'Gimap7', 'Ccnd2', 'Tnfrsf4', 'Hif1a', 'S100a11', 
                  'Ctsw', 'Jun', 'Tff1', 'Ccr7', 'Lrrc32', 'Abca1'
                  )
  )
}

all <- read_excel("~/github/2_tff1/1_pyzone/1_outputs/1_degs/7_bm_treg_vs_rest.xlsx")

# subset to significant only
all <- all[all$`Is Significant` == TRUE, , drop = FALSE]

keyvals <- initialize_color_vector(all)
keyvals <- assign_colors(all, keyvals)
keyvals <- set_color_names(all, keyvals)

all_plt <- plot_volcano(all, keyvals)

# ---- FIXED LEGEND ----
all_plt <- all_plt +
  scale_colour_manual(
    values = c(BM = "purple", Rest = "darkgray"),
    breaks = c("BM", "Rest"),
    name = NULL
  ) +
  theme(legend.position = "right")

all_plt <- all_plt +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

out_path_pdf <- "~/github/2_tff1/1_pyzone/1_outputs/0_figures/volcanoplot_bm_vs_rest.pdf"

ggsave(
  filename = out_path_pdf,
  plot = all_plt,
  device = cairo_pdf,
  width = 8, height = 7, units = "in",
  bg = "transparent"
)

all_plt

