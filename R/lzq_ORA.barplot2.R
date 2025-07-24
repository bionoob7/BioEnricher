#' @export
lzq_ORA.barplot2 <- function (enrich.obj1, enrich.obj2, Selct.P = "FDR", cutoff.P = 0.05, 
    obj.types = c("Up", "Down"), obj.type.colors = c("#ED6355", 
        "#3E94B5"), obj1.top.pathway.num = 10, obj2.top.pathway.num = 10, 
    bar.width = 0.6, add.bar.border = T, x.limit.fold = 1.05, 
    label.size = 3.5, legend.position = "bottom", use.Chinese = F, 
    appid = "20231122001888718", key = "5GpDqe8F3pmXfnOkEKGQ") 
{
    if (Selct.P == "BP") {
        r1 <- enrich.obj1@result %>% dplyr::filter(pvalue < cutoff.P) %>% 
            dplyr::arrange(pvalue) %>% dplyr::mutate(Type = obj.types[1])
        obj1.top.pathway.num <- ifelse(nrow(r1) >= obj1.top.pathway.num, 
            obj1.top.pathway.num, nrow(r1))
        if (obj1.top.pathway.num == 0) {
            stop("No significant term was enriched in enrich.obj1!")
        }
        r1 <- r1[seq_len(obj1.top.pathway.num), ]
        r2 <- enrich.obj2@result %>% dplyr::filter(pvalue < cutoff.P) %>% 
            dplyr::arrange(pvalue) %>% dplyr::mutate(Type = obj.types[2])
        obj2.top.pathway.num <- ifelse(nrow(r2) >= obj2.top.pathway.num, 
            obj2.top.pathway.num, nrow(r2))
        if (obj2.top.pathway.num == 0) {
            stop("No significant term was enriched in enrich.obj2!")
        }
        r2 <- r2[seq_len(obj2.top.pathway.num), ]
    }
    if (Selct.P == "FDR") {
        r1 <- enrich.obj1@result %>% dplyr::filter(p.adjust < 
            cutoff.P) %>% dplyr::arrange(p.adjust) %>% dplyr::mutate(Type = obj.types[1])
        obj1.top.pathway.num <- ifelse(nrow(r1) >= obj1.top.pathway.num, 
            obj1.top.pathway.num, nrow(r1))
        if (obj1.top.pathway.num == 0) {
            stop("No significant term was enriched in enrich.obj1!")
        }
        r1 <- r1[1:obj1.top.pathway.num, ]
        r2 <- enrich.obj2@result %>% dplyr::filter(p.adjust < 
            cutoff.P) %>% dplyr::arrange(p.adjust) %>% dplyr::mutate(Type = obj.types[2])
        obj2.top.pathway.num <- ifelse(nrow(r2) >= obj2.top.pathway.num, 
            obj2.top.pathway.num, nrow(r2))
        if (obj2.top.pathway.num == 0) {
            stop("No significant term was enriched in enrich.obj2!")
        }
        r2 <- r2[1:obj2.top.pathway.num, ]
    }
    rr <- rbind(r1[, intersect(colnames(r1), colnames(r2))], 
        r2[, intersect(colnames(r1), colnames(r2))]) %>% dplyr::mutate(Description = Hmisc::capitalize(Description))
    if (use.Chinese) {
        rr$Description <- purrr::map_vec(rr$Description, ~lzq_translate(.x, 
            appid = appid, key = key))
        showtext::showtext_auto()
    }
    if (Selct.P == "BP") {
        rr %<>% dplyr::mutate(log10P = log10(pvalue)) %<>% dplyr::mutate(log10P = ifelse(Type == 
            obj.types[1], -log10P, log10P))
        x.lab <- bquote(~-Log[10] ~ italic("P-value"))
    }
    if (Selct.P == "FDR") {
        rr %<>% dplyr::mutate(log10P = log10(p.adjust)) %<>% 
            dplyr::mutate(log10P = ifelse(Type == obj.types[1], 
                -log10P, log10P))
        x.lab <- bquote(~-Log[10] ~ "FDR")
    }
    rr %<>% dplyr::mutate(Type = factor(Type, obj.types)) %<>% 
        dplyr::arrange(desc(log10P)) %<>% dplyr::distinct(Description, 
        .keep_all = T) %>% dplyr::mutate(Description = factor(Description, 
        rev(Description)))
    ggplot(rr, aes(log10P, Description, fill = Type)) + geom_bar(stat = "identity", 
        width = bar.width, color = ifelse(add.bar.border, "black", 
            NA)) + scale_x_continuous(limits = c(-max(abs(rr$log10P)) * 
        x.limit.fold, max(abs(rr$log10P)) * x.limit.fold), labels = abs, 
        name = x.lab) + theme_classic(base_line_size = 0.9) + 
        theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
            axis.text.x = element_text(size = 10, colour = "black"), 
            axis.title.x = element_text(size = 13, colour = "black", 
                face = "bold"), legend.position = legend.position, 
            legend.background = element_blank(), legend.key = element_blank(), 
            legend.title = element_blank(), legend.text = element_text(size = 13, 
                colour = "black")) + geom_text(data = subset(rr, 
        Type == obj.types[2]), aes(x = 0, y = Description, label = paste0(" ", 
        Description)), size = label.size, hjust = 0) + geom_text(data = subset(rr, 
        Type == obj.types[1]), aes(x = -0.1, y = Description, 
        label = Description), size = label.size, hjust = 1) + 
        scale_fill_manual(values = obj.type.colors)
}
