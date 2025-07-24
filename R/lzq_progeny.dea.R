#' @export
lzq_progeny.dea <- function (progeny.res, groups, control.group, theme.plot = theme_classic(base_line_size = 0.8)) 
{
    groups <- ifelse(groups == control.group, T, F)
    result <- apply(progeny.res, 2, function(x) {
        broom::tidy(stats::lm(x ~ !groups)) %>% filter(term == 
            "!groupsTRUE") %>% dplyr::select(-term)
    }) %>% Reduce(rbind, .)
    result <- cbind(pathways = colnames(progeny.res), result)
    t_value <- stats::qt(1 - 0.05/2, nrow(progeny.res) - 2)
    p <- ggplot(result, aes(y = stats::reorder(pathways, statistic), 
        x = statistic)) + geom_vline(xintercept = t_value, lty = 2, 
        col = "grey70", lwd = 0.65) + geom_vline(xintercept = -t_value, 
        lty = 2, col = "grey70", lwd = 0.65) + geom_segment(aes(x = 0, 
        xend = statistic, y = stats::reorder(pathways, statistic), 
        yend = stats::reorder(pathways, statistic)), color = "grey80", 
        size = 0.8) + geom_point(aes(color = pathways, size = -log10(p.value))) + 
        theme_bw() + labs(x = "Statistic", y = NULL) + scale_color_manual(values = paletteer::paletteer_d("ggsci::category20_d3", 
        n = 14)) + scale_size(range = c(2, 6)) + theme.plot + 
        theme(axis.text.x = element_text(size = 10, colour = "black"), 
            axis.title.x = element_text(size = 13, colour = "black", 
                face = "bold"), axis.text.y = element_text(size = 12, 
                colour = "black"), panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white"), 
            plot.title = element_text(hjust = 0.5, size = 14, 
                colour = "black", face = "bold"), legend.position = "none")
    print(p)
    return(list(res = result, plot = p))
}
