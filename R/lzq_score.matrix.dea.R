#' @export
lzq_score.matrix.dea <- function (score.matrix, groups, control.group, Select.P = "FDR", 
    cutoff.P = 0.05, cutoff.logFC = 2, ...) 
{
    groups <- ifelse(groups == control.group, "A", "B")
    res <- lzq_DEA.limma(score.matrix, groups, contrasts = "B-A", 
        Select.P = Select.P, cutoff.P = cutoff.P, cutoff.logFC = cutoff.logFC)
    p <- lzq_volcano(res, Select.P = Select.P, cutoff.P = cutoff.P, 
        cutoff.logFC = cutoff.logFC, ...)
    print(p)
    return(list(res = res, plot = p))
}
