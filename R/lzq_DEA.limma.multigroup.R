#' @export
lzq_DEA.limma.multigroup <- function (expr, groups, voom = F, Select.P = "FDR", cutoff.P = 0.05, 
    cutoff.logFC = 1) 
{
    lzq_limma <- ifelse(voom, lzq_DEA.limma_voom, lzq_DEA.limma)
    res <- purrr::map(unique(groups), function(x) {
        gg <- ifelse(groups == x, "A", "B")
        cat(paste0(x, " vs Others:"))
        d <- lzq_limma(expr, gg, contrasts = "A-B", Select.P, 
            cutoff.P, cutoff.logFC)
        cat("\n")
        return(d)
    })
    names(res) <- unique(groups)
    return(res)
}
