# Stat-table savers for the beta-diversity scripts (HTML, xlsx)

pacman::p_load(kableExtra, writexl, update = FALSE)

# Stat table -> styled HTML
# - as.data.frame() first: on adonis2() output, rownames_to_column() returns row numbers otherwise
# - rownames_to = NULL for tables without meaningful rownames
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

# Named list of stat tables -> one .xlsx, one sheet per element
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
