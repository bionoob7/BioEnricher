#' @export
lzq_DEA.limma_voom <- function (counts, groups, contrasts, Select.P = "FDR", cutoff.P = 0.05, 
    cutoff.logFC = 1) 
{
    k <- strsplit(contrasts, "-")[[1]]
    groups <- ifelse(groups == k[1], "V1", "V2")
    design <- stats::model.matrix(~0 + factor(groups))
    colnames(design) <- levels(factor(groups))
    rownames(design) <- colnames(counts)
    norm <- limma::voom(counts, design)
    fit <- limma::lmFit(norm, design, method = "ls")
    contrast <- limma::makeContrasts(V1 - V2, levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast)
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
