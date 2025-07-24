#' @export
lzq_WGCNA.cor.plot <- function (WGCNA.res, module, trait, point.size = 2, theme.plot = theme_bw(base_rect_size = 1.5)) 
{
    genes <- WGCNA.res$module.gene.list[[module]]
    cat(paste0("The ", module, " module includes ", length(genes), 
        " genes."))
    MM <- WGCNA.res$MM[genes, paste0("ME", module)]
    GS <- WGCNA.res$GS[genes, trait]
    fit <- stats::cor.test(MM, GS)
    lb <- paste0("r = ", signif(fit$estimate, 3), "\n", "p = ", 
        signif(fit$p.value, 3))
    if (fit$estimate > 0) {
        gt <- geom_text(aes(min(MM), max(GS)), label = lb, size = 4, 
            vjust = 1, hjust = 0, fontface = "italic")
    }
    else {
        gt <- geom_text(aes(max(MM), max(GS)), label = lb, size = 4, 
            vjust = 1, hjust = 1, fontface = "italic")
    }
    data <- data.frame(MM = MM, GS = GS)
    ggplot(data, aes(MM, GS)) + geom_point(color = module, size = point.size) + 
        gt + labs(x = "Module membership", y = "Gene significance", 
        title = paste0(Hmisc::capitalize(module), " module (MM vs GS)")) + 
        theme.plot + theme(panel.grid = element_blank(), plot.title = element_text(size = 14, 
        hjust = 0.5, colour = "black", face = "bold"), axis.text.x = element_text(size = 10, 
        colour = "black"), axis.text.y = element_text(size = 10, 
        colour = "black"), axis.title.x = element_text(size = 13, 
        colour = "black", face = "bold"), axis.title.y = element_text(size = 13, 
        colour = "black", face = "bold"))
}
