#' @export
lzq_updateSymbolforDL <- function (data, unmapGene_keep = F) 
{
    ann <- HGNChelper::checkGeneSymbols(rownames(data), unmapped.as.na = unmapGene_keep) %>% 
        stats::na.omit() %>% dplyr::select(-Approved)
    colnames(ann) <- c("ID", "Symbol")
    if (sum(grepl("///", ann$Symbol)) > 0) {
        ann1 <- ann[!grepl("///", ann$Symbol), ]
        ann2 <- ann[grepl("///", ann$Symbol), ]
        ann2 <- lapply(ann2$Symbol, function(x) {
            l <- unlist(strsplit(x, " /// "))
            tmp <- ann2[ann2$Symbol == x, ] %>% dplyr::distinct(ID, 
                Symbol, .keep_all = T)
            tmp2 <- lapply(l, function(d) {
                kk <- tmp
                kk$Symbol <- d
                return(kk)
            }) %>% Reduce(rbind, .)
            return(tmp2)
        }) %>% Reduce(rbind, .) %>% dplyr::distinct(ID, Symbol, 
            .keep_all = T)
        ann <- rbind(ann1, ann2) %>% dplyr::distinct(ID, Symbol, 
            .keep_all = T)
    }
    data2 <- merge(ann, data, by.x = 1, by.y = 0)[, -1] %>% dplyr::group_by(Symbol) %>% 
        dplyr::summarise_all(mean) %>% tibble::column_to_rownames("Symbol")
    return(data2)
}
