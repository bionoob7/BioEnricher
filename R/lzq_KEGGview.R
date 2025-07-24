#' @export
lzq_KEGGview <- function (gene.data = NULL, gene.type = "SYMBOL", pathway.id, 
    species = "hsa", figure.suffix = "") 
{
    pathview::pathview(gene.data = gene.data, gene.idtype = gene.type, 
        pathway.id = gsub("\\D", "", pathway.id), species = species, 
        kegg.native = T, out.suffix = figure.suffix, low = list(gene = "#0a9396", 
            cpd = "blue"), mid = list(gene = "#e9d8a6", cpd = "gray"), 
        high = list(gene = "#bb3e03", cpd = "yellow"), )
    imgpng <- png::readPNG(paste0(pathway.id, ".", figure.suffix, 
        ".png"))
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::rasterImage(imgpng, 0, 0, 1, nrow(imgpng)/ncol(imgpng))
}
