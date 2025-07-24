#' @export
lzq_gseaplot <- function (GSEA.result, Pathway.ID, heatbar = T, rank = T, line.color = "#41A98E", 
    rank.colors = viridis::viridis(10), heatbar.colors = c(rev(RColorBrewer::brewer.pal(5, 
        "Blues")), RColorBrewer::brewer.pal(5, "Reds")), add.x.ann = T, 
    x.lab = "Gene ranks", line.y.lab = "Enrichment score", rank.y.lab = "logFC", 
    statistic.position = c(0.5, 0.2), statistic.face = "italic", 
    statistic.size = 3.5, rel.heights = c(1.5, 0.2, 1), theme.plot = theme_bw(base_rect_size = 1.5)) 
{
    gsdata <- enrichplot:::gsInfo(GSEA.result, Pathway.ID)
    label <- paste0("NES = ", sprintf("%.3f", GSEA.result@result$NES[GSEA.result@result$ID == 
        Pathway.ID]), "\nFDR = ", ifelse(GSEA.result@result$p.adjust[GSEA.result@result$ID == 
        Pathway.ID] < 0.001, format(GSEA.result@result$p.adjust[GSEA.result@result$ID == 
        Pathway.ID], scientific = T, digit = 3), sprintf("%.3f", 
        GSEA.result@result$p.adjust[GSEA.result@result$ID == 
            Pathway.ID])))
    plotlist <- list()
    plotlist[["line"]] <- ggplot(gsdata, aes(x)) + geom_line(aes(y = runningScore), 
        linewidth = 1, color = line.color) + scale_x_continuous(expand = c(0, 
        0)) + scale_y_continuous(expand = c(0, 0.01)) + labs(x = NULL, 
        y = line.y.lab, title = Hmisc::capitalize(GSEA.result@result$Description[GSEA.result@result$ID == 
            Pathway.ID])) + theme.plot + annotate("text", x = nrow(gsdata) * 
        statistic.position[1], y = min(gsdata$runningScore) + 
        (max(gsdata$runningScore) - min(gsdata$runningScore)) * 
            statistic.position[2], label = label, hjust = 0, 
        fontface = statistic.face, size = statistic.size) + theme(plot.title = element_text(hjust = 0.5, 
        size = 13, colour = "black", face = "bold"), panel.grid = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10, 
            colour = "black"), axis.title.y = element_text(size = 13, 
            colour = "black", face = "bold"), axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), plot.margin = margin(t = 0.2, 
            r = 0.2, b = -0.07, l = 0.2, unit = "cm"))
    if (heatbar) {
        plotlist[["heatbar"]] <- ggplot(gsdata, aes(x)) + geom_linerange(aes(ymin = ymin, 
            ymax = ymax), color = "grey30") + labs(x = NULL, 
            y = NULL) + theme.plot + theme(legend.position = "none", 
            panel.grid = element_blank(), plot.margin = margin(t = -0.1, 
                b = 0, unit = "cm"), axis.ticks = element_blank(), 
            axis.text = element_blank(), axis.line.x = element_blank()) + 
            scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
            0))
        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) {
            inv <- inv + 1
        }
        col <- heatbar.colors
        ymin <- min(plotlist[["heatbar"]]$data$ymin)
        yy <- max(plotlist[["heatbar"]]$data$ymax - plotlist[["heatbar"]]$data$ymin) * 
            0.3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
            xmax = xmax, col = col[unique(inv)])
        plotlist[["heatbar"]] <- plotlist[["heatbar"]] + geom_rect(aes(xmin = xmin, 
            xmax = xmax, ymin = ymin, ymax = 0, fill = I(col)), 
            data = d, alpha = 0.9, inherit.aes = FALSE)
    }
    if (rank) {
        plotlist[["rank"]] <- ggplot(gsdata, aes(x)) + labs(y = rank.y.lab) + 
            scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
            0)) + geom_segment(aes(x = x, xend = x, y = geneList, 
            yend = 0, color = geneList)) + scale_color_gradientn(colours = rank.colors) + 
            theme.plot + theme(plot.title = element_text(hjust = 0.5, 
            size = 13, colour = "black", face = "bold"), panel.grid = element_blank(), 
            legend.position = "none", plot.margin = margin(t = -0.17, 
                r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
    }
    n <- length(plotlist)
    if (add.x.ann) {
        plotlist[[n]] <- plotlist[[n]] + labs(x = x.lab) + theme(axis.line.x = element_line(), 
            axis.ticks.x = element_line(), axis.text.x = element_text(size = 10, 
                colour = "black"), axis.title.x = element_text(size = 13, 
                colour = "black", face = "bold"))
    }
    if (n == 3) {
        print(cowplot::plot_grid(plotlist = plotlist, ncol = 1, 
            align = "v", rel_heights = rel.heights))
    }
    if (n == 2 & heatbar) {
        rel.heights <- rel.heights[seq_len(2)]
        print(cowplot::plot_grid(plotlist = plotlist, ncol = 1, 
            align = "v", rel_heights = rel.heights))
    }
    if (n == 2 & rank) {
        rel.heights <- rel.heights[seq_len(2)]
        print(cowplot::plot_grid(plotlist = plotlist, ncol = 1, 
            align = "v", rel_heights = rel.heights))
    }
}
