# Mechanical I/O wrappers for the beta-diversity scripts (16S/ITS_beta_div_
# main/exploration.R): figure saving (size/dpi/bg) and stat-table saving
# (kable HTML). No plotting or stats logic lives here.

pacman::p_load(kableExtra, writexl, update = FALSE)

# PERMANOVA/dispersion table -> styled kable HTML.
# as.data.frame() before rownames_to_column() matters: adonis2()'s return
# class (anova.cca/anova/data.frame) makes rownames_to_column() silently
# return row NUMBERS instead of the real term names otherwise.
#
# rownames_to: NULL for tables that already have meaningful columns and no
# informative rownames (e.g. arrow coordinates, correlation stats); default
# "Parameter" is for PERMANOVA-style tables.
save_stat_kable <- function(path, tbl, caption = NULL, rownames_to = "Parameter"){
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tbl <- as.data.frame(tbl)
  if (!is.null(rownames_to)) tbl <- tibble::rownames_to_column(tbl, rownames_to)
  tbl %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) %>%
    kableExtra::kable("html", align = "l", caption = caption) %>%
    kableExtra::kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
    kableExtra::save_kable(file = path)
  invisible(path)
}

# Named list of stat tables -> one .xlsx, one sheet per element (sheet = list
# name). Same table prep as save_stat_kable().
save_stat_xlsx <- function(path, tbl_list, rownames_to = "Parameter"){
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sheets <- purrr::map(tbl_list, function(tbl){
    tbl <- as.data.frame(tbl)
    if (!is.null(rownames_to)) tbl <- tibble::rownames_to_column(tbl, rownames_to)
    dplyr::mutate(tbl, dplyr::across(where(is.numeric), ~ signif(.x, 3)))
  })
  writexl::write_xlsx(sheets, path)
  invisible(path)
}

# mm -> inches, for grDevices::pdf()/pdf_device() calls that take width/height in inches.
mm_to_in <- function(mm) mm / 25.4
