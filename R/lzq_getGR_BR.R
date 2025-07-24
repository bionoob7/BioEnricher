#' @export
lzq_getGR_BR <- function (res) 
{
    res@result$GeneRatio <- apply(res@result, 1, function(x) {
        eval(parse(text = x["GeneRatio"]))
    })
    res@result$BgRatio <- apply(res@result, 1, function(x) {
        eval(parse(text = x["BgRatio"]))
    })
    return(res)
}
