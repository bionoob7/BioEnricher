#' @export
lzq_getEF <- function (res) 
{
    res@result$EnrichmentFactor <- apply(res@result, 1, function(x) {
        GeneRatio <- eval(parse(text = x["GeneRatio"]))
        BgRatio <- eval(parse(text = x["BgRatio"]))
        EF <- round(GeneRatio/BgRatio, 2)
    })
    return(res)
}
