#' @export
lzq_interaction.1toM.plot <- function (data, point.colors = c("#00b4d8", "#e76f51"), point.shapes = c("circle", 
    "square"), point.sizes = c(15, 10), line.width = 4, line.color = "#BEBEBE99", 
    label.colors = c("black", "black"), label.sizes = c(1, 0.8), 
    one.circle.layout = F, seed = 123) 
{
    colnames(data) <- c("source", "target", "score")
    normalize <- function(x) {
        x2 <- (x - min(x))/(max(x) - min(x)) * 0.9 + 0.1
        return(x2)
    }
    data$score <- normalize(abs(data$score)) * line.width
    net <- igraph::graph_from_data_frame(data, directed = FALSE)
    igraph::V(net)$color <- ifelse(igraph::V(net)$name == data$source[1], 
        point.colors[1], point.colors[2])
    igraph::E(net)$width <- data$score
    igraph::E(net)$color <- line.color
    igraph::V(net)$frame.color <- igraph::V(net)$color
    igraph::V(net)$size <- ifelse(igraph::V(net)$name == data$source[1], 
        point.sizes[1], point.sizes[2])
    igraph::V(net)$shape <- ifelse(igraph::V(net)$name == data$source[1], 
        point.shapes[1], point.shapes[2])
    igraph::V(net)$label.cex <- ifelse(igraph::V(net)$name == 
        data$source[1], label.sizes[1], label.sizes[2])
    igraph::V(net)$label.color <- ifelse(igraph::V(net)$name == 
        data$source[1], label.colors[1], label.colors[2])
    igraph::V(net)$label.family <- "Arial"
    set.seed(seed)
    if (one.circle.layout) {
        layout <- matrix(ncol = 2, nrow = igraph::vcount(net))
        center_x <- 0
        center_y <- 0
        angle <- 2 * pi/24
        for (i in 1:igraph::vcount(net)) {
            if (igraph::V(net)$name[i] == data$source[1]) {
                layout[i, ] <- c(center_x, center_y)
            }
            else {
                layout[i, ] <- c(center_x + cos((i - 2) * angle), 
                  center_y + sin((i - 2) * angle))
            }
        }
        plot(net, layout = layout, vertex.size = igraph::V(net)$size, 
            edge.arrow.size = 0.5)
    }
    else {
        plot(net, vertex.size = igraph::V(net)$size, edge.arrow.size = 0.5)
    }
}
