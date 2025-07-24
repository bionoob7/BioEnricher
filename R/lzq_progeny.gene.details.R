#' @export
lzq_progeny.gene.details <- function (dea.table, pathway, organism = "Human", top = 100, 
    y.lab = bquote(~Log[2] ~ "(Fold change)"), colors = c("#3E94B5", 
        "grey70", "#ED6355"), point.size = 2, label.size = 4, 
    theme.plot = theme_classic(base_line_size = 0.8)) 
{
    colnames(dea.table) <- c("Gene", "Stat")
    if (organism == "Human") {
        model <- progeny::model_human_full
    }
    if (organism == "Mouse") {
        model <- progeny::model_mouse_full
    }
    model %<>% dplyr::group_by(pathway) %>% dplyr::slice_min(order_by = p.value, 
        n = top)
    model <- model[model$pathway == pathway, ]
    model2 <- merge(model, dea.table, by = 1) %>% dplyr::mutate(color = "A") %>% 
        dplyr::mutate(color = ifelse(weight > 0 & Stat > 0, "C", 
            color)) %>% dplyr::mutate(color = ifelse(weight > 
        0 & Stat <= 0, "B", color)) %>% dplyr::mutate(color = ifelse(weight <= 
        0 & Stat > 0, "B", color)) %>% dplyr::mutate(color = ifelse(weight <= 
        0 & Stat <= 0, "A", color))
    model2$Gene2 <- ifelse(model2$color %in% c("A", "C"), as.character(model2$gene), 
        NA)
    p <- suppressWarnings(ggplot(model2, aes(x = weight, y = Stat, 
        color = color)) + geom_vline(xintercept = 0, lty = 2, 
        col = "grey70", lwd = 0.65) + geom_hline(yintercept = 0, 
        lty = 2, col = "grey70", lwd = 0.65) + geom_point(size = point.size) + 
        labs(x = "Weight (PROGENy)", y = y.lab, title = pathway) + 
        theme.plot + scale_colour_manual(values = colors) + ggrepel::geom_label_repel(mapping = aes(label = Gene2), 
        size = label.size, box.padding = unit(0.35, "lines"), 
        point.padding = unit(0.3, "lines"), max.overlaps = 20) + 
        theme(axis.text = element_text(size = 10, colour = "black"), 
            axis.title.x = element_text(size = 13, colour = "black", 
                face = "bold"), axis.title.y = element_text(size = 13, 
                colour = "black", face = "bold"), panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white"), 
            plot.title = element_text(hjust = 0.5, size = 14, 
                colour = "black", face = "bold"), legend.position = "none"))
    print(suppressWarnings(p))
    return(list(res = model2[, seq_len(5)], plot = p))
}
