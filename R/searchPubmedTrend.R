#' @export
searchPubmedTrend <- function (term, add_term = NULL, period) 
{
    if (!is.null(add_term)) {
        supp <- paste0("AND ABSTRACT:", add_term) %>% paste(collapse = " ")
    }
    else {
        supp <- NULL
    }
    res <- lapply(term, function(i) {
        europepmc::epmc_hits_trend(paste0("ABSTRACT:", i, " ", 
            supp), period = period, synonym = FALSE)
    })
    names(res) <- term
    return(res)
}
