# other functions
extract_si <- function(sdisca, si = "SI", fromSGSVD = FALSE) {
  if ("list" %in% class(sdisca)) {
    if (fromSGSVD) {
      return(sdisca$SI[[si]])
    }else {
      return(sdisca$sparsity$SI[[si]])
    }
  } else if (is.na(sdisca)) {
    return(NA)
  } else {
    stop("Something is a foot!")
  }
}

pivot_the_tab <- function(dat) {
  df <- data.frame(V = unname(t(data.frame(dat))))
  df %>% 
    tidyr::pivot_longer(starts_with("V"), names_to = "k", names_prefix = "V\\.") %>%
    select(value) %>% purrr::as_vector() %>% unname()
}