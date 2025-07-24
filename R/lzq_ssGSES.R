#' @export
lzq_ssGSES <- function (exp, gene.list, method = "ssgsea") 
{
    if (method == "ssgsea") {
        score <- as.data.frame(t(GSVA::gsva(as.matrix(exp), gene.list, 
            method == "ssgsea")))
    }
    if (method == "gsva") {
        score <- as.data.frame(t(GSVA::gsva(as.matrix(exp), gene.list, 
            method = "gsva", kcdf = "Gaussian", abs.ranking = F, 
            min.sz = 1, max.sz = Inf, mx.diff = T, verbose = T)))
    }
    if (method == "zscore") {
        score <- as.data.frame(t(GSVA::gsva(as.matrix(exp), gene.list, 
            method = "zscore")))
    }
    if (method == "plage") {
        score <- as.data.frame(t(GSVA::gsva(as.matrix(exp), gene.list, 
            method = "plage")))
    }
    return(score)
}
