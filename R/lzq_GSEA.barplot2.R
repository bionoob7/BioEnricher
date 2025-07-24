#' @export
lzq_GSEA.barplot2 <- function (enrich.obj, Selct.P = "FDR", cutoff.P = 0.05, types = c("Positive", 
    "Negative"), type.colors = c("#ED6355", "#3E94B5"), pos.top.pathway.num = 10, 
    neg.top.pathway.num = 10, bar.width = 0.6, add.bar.border = T, 
    x.limit.fold = 1.05, label.size = 3.5, legend.position = "bottom", 
    use.Chinese = F, appid = "20231122001888718", key = "5GpDqe8F3pmXfnOkEKGQ") 
{
    if (Selct.P == "BP") {
        r1 <- enrich.obj@result %>% dplyr::filter(pvalue < cutoff.P, 
            NES > 0) %>% dplyr::arrange(dplyr::desc(NES)) %>% 
            dplyr::mutate(Type = types[1])
        pos.top.pathway.num <- ifelse(nrow(r1) >= pos.top.pathway.num, 
            pos.top.pathway.num, nrow(r1))
        if (pos.top.pathway.num == 0) {
            stop("No significant positive term was enriched in enrich.obj!")
        }
        r1 <- r1[seq_len(pos.top.pathway.num), ]
        r2 <- enrich.obj@result %>% dplyr::filter(pvalue < cutoff.P, 
            NES < 0) %>% dplyr::arrange(NES) %>% dplyr::mutate(Type = types[2])
        neg.top.pathway.num <- ifelse(nrow(r2) >= neg.top.pathway.num, 
            neg.top.pathway.num, nrow(r2))
        if (neg.top.pathway.num == 0) {
            stop("No significant negative term was enriched in enrich.obj!")
        }
        r2 <- r2[seq_len(neg.top.pathway.num), ]
    }
    if (Selct.P == "FDR") {
        r1 <- enrich.obj@result %>% dplyr::filter(p.adjust < 
            cutoff.P, NES > 0) %>% dplyr::arrange(dplyr::desc(NES)) %>% 
            dplyr::mutate(Type = types[1])
        pos.top.pathway.num <- ifelse(nrow(r1) >= pos.top.pathway.num, 
            pos.top.pathway.num, nrow(r1))
        if (pos.top.pathway.num == 0) {
            stop("No significant positive term was enriched in enrich.obj!")
        }
        r1 <- r1[seq_len(pos.top.pathway.num), ]
        r2 <- enrich.obj@result %>% dplyr::filter(p.adjust < 
            cutoff.P, NES < 0) %>% dplyr::arrange(NES) %>% dplyr::mutate(Type = types[2])
        neg.top.pathway.num <- ifelse(nrow(r2) >= neg.top.pathway.num, 
            neg.top.pathway.num, nrow(r2))
        if (neg.top.pathway.num == 0) {
            stop("No significant negative term was enriched in enrich.obj!")
        }
        r2 <- r2[seq_len(neg.top.pathway.num), ]
    }
    rr <- rbind(r1[, intersect(colnames(r1), colnames(r2))], 
        r2[, intersect(colnames(r1), colnames(r2))]) %>% dplyr::mutate(Description = Hmisc::capitalize(Description))
    if (use.Chinese) {
        rr$Description <- purrr::map_vec(rr$Description, ~lzq_translate(.x, 
            appid = appid, key = key))
        showtext::showtext_auto()
    }
    rr %<>% dplyr::mutate(Type = factor(Type, types)) %<>% dplyr::arrange(desc(NES)) %<>% 
        dplyr::distinct(Description, .keep_all = T) %>% dplyr::mutate(Description = factor(Description, 
        rev(Description)))
    ggplot2::ggplot(rr, aes(NES, Description, fill = Type)) + 
        ggplot2::geom_bar(stat = "identity", width = bar.width, 
            color = ifelse(add.bar.border, "black", NA)) + ggplot2::scale_x_continuous(limits = c(-max(abs(rr$NES)) * 
        x.limit.fold, max(abs(rr$NES)) * x.limit.fold), name = "NES") + 
        ggplot2::theme_classic(base_line_size = 0.9) + ggplot2::theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.text.x = element_text(size = 10, 
            colour = "black"), axis.title.x = element_text(size = 13, 
            colour = "black", face = "bold"), legend.position = legend.position, 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.title = element_blank(), legend.text = element_text(size = 13, 
            colour = "black")) + ggplot2::geom_text(data = subset(rr, 
        Type == types[2]), aes(x = 0, y = Description, label = paste0(" ", 
        Description)), size = label.size, hjust = 0) + ggplot2::geom_text(data = subset(rr, 
        Type == types[1]), aes(x = -0.1, y = Description, label = Description), 
        size = label.size, hjust = 1) + ggplot2::scale_fill_manual(values = type.colors)
}
