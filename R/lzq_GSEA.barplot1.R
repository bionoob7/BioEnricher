#' @export
lzq_GSEA.barplot1 <- function (enrich.obj, type = "Positive", show.term.num = 15, 
    Selct.P = "FDR", cutoff.P = 0.05, colors = c("#003c30", "#01665e", 
        "#35978f", "#80cdc1", "#c7eae5", "#f6e8c3", "#dfc27d", 
        "#bf812d", "#8c510a", "#543005"), add.bar.border = T, 
    bar.width = 0.6, y.label.position = "right", title = NULL, 
    legend.position = "right", theme.plot = theme_bw(base_rect_size = 1.5), 
    use.Chinese = F, appid = "20231122001888718", key = "5GpDqe8F3pmXfnOkEKGQ") 
{
    if (Selct.P == "BP") {
        if (type %in% c("Positive", "positive", "pos", "p", "po", 
            "P", "Po", "Pos")) {
            r <- enrich.obj@result %>% dplyr::filter(pvalue < 
                cutoff.P, NES > 0) %>% dplyr::mutate(sig = -log10(pvalue)) %>% 
                dplyr::arrange(dplyr::desc(sig))
        }
        else {
            r <- enrich.obj@result %>% dplyr::filter(pvalue < 
                cutoff.P, NES < 0) %>% dplyr::mutate(sig = -log10(pvalue)) %>% 
                dplyr::arrange(dplyr::desc(sig))
        }
        x.lab <- bquote(~-Log[10] ~ italic("P-value"))
    }
    if (Selct.P == "FDR") {
        if (type %in% c("Positive", "positive", "pos", "p", "po", 
            "P", "Po", "Pos")) {
            r <- enrich.obj@result %>% dplyr::filter(p.adjust < 
                cutoff.P, NES > 0) %>% dplyr::mutate(sig = -log10(p.adjust)) %>% 
                dplyr::arrange(dplyr::desc(sig))
        }
        else {
            r <- enrich.obj@result %>% dplyr::filter(p.adjust < 
                cutoff.P, NES < 0) %>% dplyr::mutate(sig = -log10(p.adjust)) %>% 
                dplyr::arrange(dplyr::desc(sig))
        }
        x.lab <- bquote(~-Log[10] ~ "FDR")
    }
    show.term.num <- ifelse(nrow(r) >= show.term.num, show.term.num, 
        nrow(r))
    if (show.term.num == 0) {
        stop(paste0("No significant ", tolower(type), " term was enriched in enrich.obj!"))
    }
    r <- r[seq_len(show.term.num), ]
    if (use.Chinese) {
        r$Description <- purrr::map_vec(r$Description, ~lzq_translate(.x, 
            appid = appid, key = key))
        showtext::showtext_auto()
    }
    r %<>% dplyr::mutate(Description = factor(Description, rev(Description)))
    if (!type %in% c("Positive", "positive", "pos", "p", "po", 
        "P", "Po", "Pos")) {
        color.title <- "abs(NES)"
    }
    else {
        color.title <- "NES"
    }
    ggplot(r, aes(sig, Description, fill = abs(NES))) + geom_bar(stat = "identity", 
        width = bar.width, color = ifelse(add.bar.border, "black", 
            NA)) + labs(fill = color.title, x = x.lab, y = NULL, 
        title = title) + scale_fill_gradientn(colours = colors) + 
        scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position) + 
        theme.plot + theme(axis.text.y = element_text(size = 13, 
        colour = "black"), axis.title.y = element_blank(), axis.text.x = element_text(size = 10, 
        colour = "black"), axis.title.x = element_text(size = 13, 
        colour = "black", face = "bold"), plot.title = element_text(size = 14, 
        colour = "black", face = "bold", hjust = 0.5), legend.position = legend.position, 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 13, colour = "black", 
            face = "bold"), legend.text = element_text(size = 11, 
            colour = "black"))
}
