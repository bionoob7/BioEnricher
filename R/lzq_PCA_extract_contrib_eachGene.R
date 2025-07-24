#' @export
lzq_PCA_extract_contrib_eachGene <- function (expr, scale = T, axes = 1) 
{
    if (scale) {
        expr <- t(scale(t(expr)))
    }
    pca <- FactoMineR::PCA(expr, graph = F)
    k <- factoextra::fviz_contrib(pca, choice = "ind", axes = 1)
    d <- k$data
    d <- d[order(d$contrib, decreasing = T), ]
    d2 <- d$contrib
    names(d2) <- d$name
    return(d2)
}
