
# themes
ggplot2::theme_set(
  ggplot2::theme_light() +
    theme(
      strip.text = element_text(colour = 'grey20'),
      strip.background = element_rect(
        fill = 'grey90',  linewidth = 0.2),
      panel.border = element_rect(colour = 'grey90'),
      #panel.grid = element_blank(),
      panel.spacing = unit(0, "lines")
    ))


## Themes -------------------------------------------------------------------------------
theme_pcoa <- theme(
  legend.position = "right",
  #panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(color = "white", fill = "white"),
  panel.border = element_rect(color = "black"),
  plot.background = element_rect(color = "white", fill = "white")
)

# Barplots (currently only used by 0.3.1 meteo)
theme_barplot <- theme(
  legend.title = element_text(size=18),
  legend.text = element_text(size=16, face = "italic"),
  plot.title = element_text(size = 20),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(color = "white", fill = "white"),
  panel.border = element_rect(color = "black"),
  plot.background = element_rect(color = "white", fill = "white"),
  legend.key.height = grid::unit(0.45, "cm"))
