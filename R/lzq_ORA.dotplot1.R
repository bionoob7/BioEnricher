#' @export
lzq_ORA.dotplot1 <- function (enrich.obj, x = "GeneRatio", show.term.num = 15, color.by = "p.adjust", 
    colors = c("#003c30", "#01665e", "#35978f", "#80cdc1", "#c7eae5", 
        "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005"), 
    color.title = color.by, size.by = "Count", size.range = c(3, 
        8), size.title = size.by, y.label.position = "right", 
    title = NULL, legend.position = "right", theme.plot = theme_bw(base_rect_size = 1.5), 
    use.Chinese = F, appid = "20231122001888718", key = "5GpDqe8F3pmXfnOkEKGQ") 
{
    enrich.obj <- lzq_getGR_BR(enrich.obj)
    x.lab <- x
    if (x == "pvalue") {
        enrich.obj@result$Sig <- -log10(enrich.obj@result$pvalue)
        x.lab <- bquote(~-Log[10] ~ italic("P-value"))
        x <- "Sig"
    }
    if (x == "p.adjust") {
        enrich.obj@result$Sig <- -log10(enrich.obj@result$p.adjust)
        x.lab <- bquote(~-Log[10] ~ "FDR")
        x <- "Sig"
    }
    if (x == "EnrichmentFactor") {
        x.lab <- "Enrichment factor"
    }
    if (color.by == "pvalue") {
        enrich.obj@result$SigL <- -log10(enrich.obj@result$pvalue)
        color.title <- bquote(~-Log[10] ~ italic("P-value"))
        color.by <- "SigL"
    }
    if (color.by == "p.adjust") {
        enrich.obj@result$SigL <- -log10(enrich.obj@result$p.adjust)
        color.title <- bquote(~-Log[10] ~ "FDR")
        color.by <- "SigL"
    }
    if (size.by == "pvalue") {
        enrich.obj@result$SigL2 <- -log10(enrich.obj@result$pvalue)
        size.title <- bquote(~-Log[10] ~ italic("P-value"))
        size.by <- "SigL2"
    }
    if (size.by == "p.adjust") {
        enrich.obj@result$SigL2 <- -log10(enrich.obj@result$p.adjust)
        size.title <- bquote(~-Log[10] ~ "FDR")
        size.by <- "SigL2"
    }
    show.term.num <- ifelse(nrow(enrich.obj@result) >= show.term.num, 
        show.term.num, nrow(enrich.obj@result))
    dd <- enrich.obj@result %>% dplyr::arrange(pvalue) %>% .[seq_len(show.term.num), 
        ] %>% dplyr::arrange(desc(get(x)))
    if (use.Chinese) {
        dd$Description <- purrr::map_vec(dd$Description, ~lzq_translate(.x, 
            appid = appid, key = key))
        dd <- dplyr::mutate(dd, Description = factor(Description, 
            rev(Description)))
        showtext::showtext_auto()
    }
    else {
        dd <- dplyr::mutate(dd, Description = factor(Description, 
            rev(Description)))
    }
    suppressWarnings(ggplot(dd, aes_string(x, "Description", 
        fill = color.by, size = size.by)) + geom_point(shape = 21, 
        color = "black") + theme.plot + labs(fill = color.title, 
        size = size.title, x = x.lab, y = NULL, title = title) + 
        scale_fill_gradientn(colours = colors) + scale_size(range = size.range) + 
        scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position) + 
        theme(axis.text.x = element_text(size = 10, colour = "black"), 
            axis.title.x = element_text(size = 13, colour = "black", 
                face = "bold"), axis.text.y = element_text(size = 13, 
                colour = "black"), axis.title = element_blank(), 
            panel.background = element_blank(), legend.text = element_text(size = 11, 
                colour = "black"), legend.title = element_text(size = 13, 
                colour = "black", face = "bold"), legend.background = element_blank(), 
            legend.key = element_blank(), legend.position = legend.position, 
            plot.title = element_text(hjust = 0.5, size = 14, 
                colour = "black", face = "bold")))
}
