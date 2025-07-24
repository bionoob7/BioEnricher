#' @export
lzq_PCAplot <- function (expr, scale = T, Group, levels, cols = c("#3E94B5", 
    "#ED6355"), rect_size = 1.5, point_size = 3, point_alpha = 1, 
    point_stroke = 1, point_border_col = "black", ellipse_level = 0.99, 
    ellipse_fill_alpha = 0.2, ellipse_linewidth = 0.8, legend_position = "right") 
{
    if (scale) {
        expr <- t(scale(t(expr)))
    }
    pca <- FactoMineR::PCA(t(expr), graph = F)
    pca2 <- as.data.frame(pca$ind$coord)
    pca2$Group <- factor(Group, levels = levels)
    pca2 <- pca2[order(pca2$Group, decreasing = T), ]
    pp <- ggplot(pca2, aes(Dim.1, Dim.2)) + theme_bw(base_rect_size = rect_size) + 
        stat_ellipse(aes(color = Group), geom = "polygon", fill = NA, 
            linewidth = ellipse_linewidth, level = ellipse_level) + 
        stat_ellipse(aes(fill = Group), geom = "polygon", alpha = ellipse_fill_alpha, 
            level = ellipse_level) + geom_point(aes(fill = Group), 
        shape = 21, size = point_size, color = point_border_col, 
        alpha = point_alpha, stroke = point_stroke) + labs(x = "PC1", 
        y = "PC2") + scale_fill_manual(values = cols) + scale_color_manual(values = cols) + 
        theme(axis.text = element_text(size = 10, colour = "black"), 
            axis.title.x = element_text(size = 13, colour = "black", 
                face = "bold"), axis.title.y = element_text(size = 13, 
                colour = "black", face = "bold"), panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white"), 
            plot.title = element_text(hjust = 0.5, size = 14, 
                colour = "black", face = "bold"), legend.position = legend_position, 
            legend.background = element_blank(), legend.key = element_blank(), 
            legend.title = element_blank(), legend.text = element_text(size = 13, 
                colour = "black"))
    print(pp)
    return(list(PCA_res = pca, PCA_plot = pp))
}
