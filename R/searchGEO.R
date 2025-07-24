#' @export
searchGEO <- function (searchterm, minnum = 0, maxnum = 1000) 
{
    searchterm <- stringr::str_trim(searchterm, side = "both") %>% 
        stringr::str_replace_all(pattern = "\\s+", replacement = "+")
    Terms <- stringr::str_split(searchterm, "\\+") %>% unlist()
    esearch <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=searchterm&retmax=maxnum&retstart=minnum&usehistory=y"
    esearch <- stringr::str_replace(esearch, "searchterm", searchterm) %>% 
        stringr::str_replace("maxnum", "0") %>% stringr::str_replace("minnum", 
        as.character(minnum))
    esearch_res <- httr::GET(esearch) %>% XML::xmlParse() %>% 
        XML::xmlToList()
    esummary <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&version=2.0&retmax=maxnum&retstart=minnum&usehistory=y&query_key=KEY&WebEnv=WEBENV"
    esummary <- stringr::str_replace(esummary, "KEY", esearch_res$QueryKey) %>% 
        stringr::str_replace("WEBENV", esearch_res$WebEnv) %>% 
        stringr::str_replace("minnum", as.character(minnum)) %>% 
        stringr::str_replace("maxnum", as.character(maxnum))
    esummary_res <- httr::GET(esummary) %>% xml2::read_xml() %>% 
        xml2::as_list()
    esummary_list <- esummary_res$eSummaryResult$DocumentSummarySet
    esummary_list <- esummary_list[2:length(esummary_list)]
    Accession <- c()
    Date <- c()
    Title <- c()
    Summary <- c()
    Platform <- c()
    Organism <- c()
    Seqtype <- c()
    Suppfile <- c()
    Nsamples <- c()
    if (length(esummary_list) == 0) {
        res <- NULL
    }
    else {
        for (i in seq_len(length(esummary_list))) {
            doc <- esummary_list[i]$DocumentSummary
            Accession <- tryCatch({
                append(Accession, ifelse(is.null(doc$Accession[[1]]), 
                  NA, doc$Accession[[1]]))
            }, error = function(e) {
                append(Accession, NA)
            })
            Date <- tryCatch({
                append(Date, ifelse(is.null(doc$PDAT[[1]]), NA, 
                  doc$PDAT[[1]]))
            }, error = function(e) {
                append(Date, NA)
            })
            Title <- tryCatch({
                append(Title, ifelse(is.null(doc$title[[1]]), 
                  NA, doc$title[[1]]))
            }, error = function(e) {
                append(Title, NA)
            })
            Summary <- tryCatch({
                append(Summary, ifelse(is.null(doc$summary[[1]]), 
                  NA, doc$summary[[1]]))
            }, error = function(e) {
                append(Summary, NA)
            })
            Platform <- tryCatch({
                append(Platform, ifelse(is.null(doc$GPL[[1]]), 
                  NA, paste0("GPL", doc$GPL)))
            }, error = function(e) {
                append(Platform, NA)
            })
            Organism <- tryCatch({
                append(Organism, ifelse(is.null(doc$taxon[[1]]), 
                  NA, doc$taxon[[1]]))
            }, error = function(e) {
                append(Organism, NA)
            })
            Seqtype <- tryCatch({
                append(Seqtype, ifelse(is.null(doc$gdsType[[1]]), 
                  NA, doc$gdsType[[1]]))
            }, error = function(e) {
                append(Seqtype, NA)
            })
            Suppfile <- tryCatch({
                append(Suppfile, ifelse(is.null(doc$suppFile[[1]]), 
                  NA, doc$suppFile[[1]]))
            }, error = function(e) {
                append(Suppfile, NA)
            })
            Nsamples <- tryCatch({
                append(Nsamples, ifelse(is.null(doc$n_samples[[1]]), 
                  NA, doc$n_samples[[1]]))
            }, error = function(e) {
                append(Nsamples, NA)
            })
        }
        res <- data.frame(Accession = Accession, Date = Date, 
            Organism = Organism, Title = Title, Summary = Summary, 
            Seqtype = Seqtype, Suppfile = Suppfile, Nsamples = Nsamples, 
            Platform = Platform)
        for (n in tolower(Terms)) {
            res <- res[grepl(n, tolower(res$Summary)), ]
        }
    }
    return(res)
}
