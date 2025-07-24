#' @export
lzq_mfuzz <- function (data, cluster.num = 8, plot.ncol = ceiling(cluster.num/2), 
    line.colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", 
        "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#f46d43", 
        "#d73027", "#a50026"), line.width = 0.3, add.median.line = T, 
    median.line.color = "firebrick3", median.line.width = 0.7, 
    theme.plot = theme_classic(base_line_size = 0.8), y.lab = "Expression Change", 
    y.lab.fontface = "bold", x.lab = NULL, x.text.angle = 0, 
    x.text.fontface = "bold", facet.strip.background = element_blank(), 
    facet.strip.text.fontface = "bold", facet.scales = NULL, 
    legend.direction = "vertical", legend.height = 1, legend.width = 0.3, 
    legend.position = ifelse(legend.direction == "vertical", 
        "right", "bottom"), seed = 123) 
{
    mfuzz_class <- lzq_quiet(methods::new("ExpressionSet", exprs = as.matrix(data)) %>% 
        Mfuzz::filter.std(., min.std = 0, visu = F) %>% Mfuzz::standardise(.))
    set.seed(seed)
    cluster_num <- cluster.num
    m <- Mfuzz::mestimate(mfuzz_class)
    cl <- Mfuzz::mfuzz(mfuzz_class, c = cluster_num, m = m)
    raw_cluster_anno <- cbind(data[names(cl$cluster), ], cluster = cl$cluster)
    norm_cluster_anno <- cbind(mfuzz_class@assayData$exprs, cluster = cl$cluster)
    mem <- cbind(cl$membership, cluster2 = cl$cluster) %>% as.data.frame() %>% 
        dplyr::mutate(gene = rownames(.))
    membership_info <- lapply(1:(ncol(mem) - 2), function(x) {
        ms <- mem %>% dplyr::filter(cluster2 == x)
        res <- data.frame(membership = ms[[x]], gene = ms$gene, 
            cluster2 = ms$cluster2)
    }) %>% do.call("rbind", .)
    dnorm <- cbind(mfuzz_class@assayData$exprs, cluster = cl$cluster) %>% 
        as.data.frame() %>% dplyr::mutate(gene = rownames(.))
    final_res <- merge(dnorm, membership_info, by = "gene")
    final_res <- final_res[, -ncol(final_res)]
    df <- reshape2::melt(final_res, id.vars = c("cluster", "gene", 
        "membership"), variable.name = "group", value.name = "norm_value")
    df$cluster_name <- paste("Cluster ", df$cluster, sep = "")
    df$cluster_name <- factor(df$cluster_name, levels = paste0("Cluster ", 
        sort(unique(df$cluster))))
    df$membership[which.max(df$membership)] <- 1
    df$membership[which.min(df$membership)] <- 0
    p <- ggplot(df, aes(x = group, y = norm_value)) + geom_line(aes(color = membership, 
        group = gene), linewidth = line.width) + theme.plot + 
        labs(y = y.lab, x = x.lab) + scale_color_gradientn(colours = line.colors) + 
        theme(axis.ticks = element_line(linewidth = 0.5), axis.ticks.length = unit(0.15, 
            "cm"), axis.text.x = element_text(angle = x.text.angle, 
            hjust = ifelse(x.text.angle != 0, 1, 0.5), size = 12, 
            color = "black", face = x.text.fontface), axis.text.y = element_text(colour = "black", 
            size = 10), axis.title = element_text(size = 13, 
            face = y.lab.fontface, colour = "black"), strip.background = facet.strip.background, 
            strip.text = element_text(size = 14, face = facet.strip.text.fontface, 
                colour = "black"), legend.direction = legend.direction, 
            legend.title = element_blank(), legend.text = element_text(colour = "black", 
                size = 10), legend.background = element_blank(), 
            legend.position = legend.position, legend.key.height = unit(legend.height, 
                "cm"), legend.key.width = unit(legend.width, 
                "cm")) + facet_wrap(~cluster_name, ncol = plot.ncol, 
        scales = facet.scales)
    if (add.median.line) {
        p <- p + geom_line(stat = "summary", fun = "median", 
            colour = median.line.color, linewidth = median.line.width, 
            aes(group = 1))
    }
    print(p)
    res <- data.frame(Gene = names(cl$cluster), Cluster = cl$cluster, 
        row.names = NULL)
    return(list(res = res, mfuzz.cl = cl, plot = p))
}
