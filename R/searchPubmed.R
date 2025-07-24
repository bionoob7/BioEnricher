#' @export
searchPubmed <- function (term, add_term = NULL, num = 100) 
{
    if (!is.null(add_term)) {
        supp <- paste0("AND ABSTRACT:", add_term) %>% paste(collapse = " ")
    }
    else {
        supp <- NULL
    }
    res <- lapply(term, function(i) {
        europepmc::epmc_search(query = paste0("ABSTRACT:", i, 
            " ", supp), limit = num, verbose = F) %>% as.data.frame()
    })
    names(res) <- term
    return(res)
}
