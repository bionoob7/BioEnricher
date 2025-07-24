#' @export
lzq_DEA.limma <- function (expr, groups, contrasts = "Tumor-Normal", Select.P = "FDR", 
    cutoff.P = 0.05, cutoff.logFC = 1) 
{
    k <- strsplit(contrasts, "-")[[1]]
    groups <- ifelse(groups == k[1], "V1", "V2")
    design <- stats::model.matrix(~0 + factor(groups))
    colnames(design) <- levels(factor(groups))
    rownames(design) <- colnames(expr)
    contrast.matrix <- limma::makeContrasts(V1 - V2, levels = design)
    fit <- limma::lmFit(expr, design)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    DEG <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "logFC", 
        adjust.method = "fdr")
    DEG <- tibble::rownames_to_column(DEG, "Gene")
    if (Select.P == "NP") {
        DEG$Type <- ifelse(DEG$P.Value >= cutoff.P, "NoSig", 
            ifelse(DEG$logFC > cutoff.logFC, "Up", ifelse(DEG$logFC < 
                -cutoff.logFC, "Down", "NoSig")))
    }
    if (Select.P == "FDR") {
        DEG$Type <- ifelse(DEG$adj.P.Val >= cutoff.P, "NoSig", 
            ifelse(DEG$logFC > cutoff.logFC, "Up", ifelse(DEG$logFC < 
                -cutoff.logFC, "Down", "NoSig")))
    }
    print(table(DEG$Type))
    return(DEG)
}
