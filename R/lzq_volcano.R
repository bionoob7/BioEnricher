#' Volcano plot
#'
#' Creates a volcano plot for visualizing differential expression results,
#' showing the relationship between fold change and statistical significance.
#'
#' @param DEG Data frame containing differential expression results
#' @param logFC.Ncol Column number containing log fold change values. Default is 2.
#' @param Select.P Type of p-value to use ("FDR" or "pvalue"). Default is "FDR".
#' @param P.Ncol Column number containing p-values. Default is 6.
#' @param DEG.type.Ncol Column number containing differential expression type. Default is 8.
#' @param cutoff.P P-value cutoff for significance. Default is 0.05.
#' @param cutoff.logFC Log fold change cutoff for significance. Default is 1.
#' @param colors Vector of colors for down-regulated, non-significant, and up-regulated genes. Default is c("#3E94B5", "#E3E3E3", "#ED6355").
#' @param color.levels Vector of level names corresponding to colors. Default is c("Down", "NoSig", "Up").
#' @param point.maxsize Maximum point size. Default is 4.
#' @param point.alpha Point transparency. Default is 0.8.
#' @param intercept.width Width of significance threshold lines. Default is 0.65.
#' @param Gene.Ncol Column number containing gene names. Default is 1.
#' @param Select.genes Vector of specific genes to label. Default is NULL.
#' @param label.size Size of gene labels. Default is 4.
#' @param legend.position Position of legend. Default is "bottom".
#' @param theme.plot ggplot2 theme to apply. Default is theme_bw(base_rect_size = 1.5).
#'
#' @return A ggplot2 object representing the volcano plot.
#'
#' @examples
#' \dontrun{
#' # Basic volcano plot
#' volcano_plot <- lzq_volcano(DEG = my_deg_results)
#' 
#' # Volcano plot with specific genes labeled
#' volcano_plot <- lzq_volcano(DEG = my_deg_results, 
#'                            Select.genes = c("TP53", "BRCA1"))
#' }
#'
#' @export
lzq_volcano <- function (DEG, logFC.Ncol = 2, Select.P = "FDR", P.Ncol = 6, 
    DEG.type.Ncol = 8, cutoff.P = 0.05, cutoff.logFC = 1, colors = c("#3E94B5", 
        "#E3E3E3", "#ED6355"), color.levels = c("Down", "NoSig", 
        "Up"), point.maxsize = 4, point.alpha = 0.8, intercept.width = 0.65, 
    Gene.Ncol = 1, Select.genes = NULL, label.size = 4, legend.position = "bottom", 
    theme.plot = theme_bw(base_rect_size = 1.5)) 
{
    tmp <- DEG
    colnames(tmp)[c(logFC.Ncol, P.Ncol, DEG.type.Ncol)] <- c("logFC", 
        "P", "Type")
    tmp$Type <- factor(tmp$Type, levels = color.levels)
    if (!is.null(Select.genes)) {
        colnames(tmp)[Gene.Ncol] <- "Gene"
        tmp$Select_label <- ifelse(tmp$Gene %in% Select.genes, 
            tmp$Gene, NA)
        tmp2 <- tmp[stats::complete.cases(tmp$Select_label), 
            ]
    }
    p1 <- ggplot2::ggplot(tmp, aes(logFC, -log10(P))) + ggplot2::geom_point(alpha = point.alpha, 
        aes(color = Type, size = abs(logFC))) + ggplot2::geom_vline(xintercept = c(-cutoff.logFC, 
        cutoff.logFC), lty = 2, col = "grey20", lwd = intercept.width) + 
        ggplot2::geom_hline(yintercept = -log10(cutoff.P), lty = 2, 
            col = "grey20", lwd = intercept.width) + ggplot2::scale_color_manual(values = colors) + 
        ggplot2::scale_size_area(max_size = point.maxsize) + 
        theme.plot + ggplot2::labs(x = bquote(~Log[2] ~ "(Fold change)"), 
        y = switch(Select.P, NP = bquote(~-Log[10] ~ italic("P-value")), 
            FDR = bquote(~-Log[10] ~ "FDR"))) + ggplot2::xlim(-max(abs(tmp$logFC)) * 
        1.05, max(abs(tmp$logFC)) * 1.05) + ggplot2::theme(axis.text = element_text(size = 10, 
        colour = "black"), axis.title.x = element_text(size = 13, 
        colour = "black", face = "bold"), axis.title.y = element_text(size = 13, 
        colour = "black", face = "bold"), panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5, 
            size = 14, colour = "black", face = "bold"), legend.position = legend.position, 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.title = element_blank(), legend.text = element_text(size = 13, 
            colour = "black")) + ggplot2::guides(size = "none", 
        color = guide_legend(order = 0, override.aes = list(size = point.maxsize, 
            alpha = 1)))
    if (!is.null(Select.genes)) {
        p2 <- p1 + ggplot2::geom_point(data = tmp2, alpha = 1, 
            size = point.maxsize, shape = 1, stroke = 1, color = "black") + 
            ggrepel::geom_text_repel(data = tmp2, mapping = aes(label = Select_label), 
                size = label.size, box.padding = unit(0.35, "lines"), 
                point.padding = unit(0.3, "lines"), max.overlaps = 20)
    }
    if (!is.null(Select.genes)) {
        return(list(volcano = p1, Volcano_with_labels = p2))
    }
    else {
        return(p1)
    }
}
