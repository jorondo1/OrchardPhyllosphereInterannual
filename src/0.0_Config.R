
ggplot2::theme_set(ggplot2::theme_light())

strip_theme <- list(
  strip.text = element_text(colour = 'grey20'),
  strip.background = element_rect(
    fill = 'grey90', colour = 'grey90', linewidth = 0.1),
  panel.border = element_rect(colour = 'grey90'),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "lines")   # default is 0.5, shrink gutters between panels
)
