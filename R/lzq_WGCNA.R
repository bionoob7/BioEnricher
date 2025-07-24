#' @export
lzq_WGCNA <- function (exp, pheno, net.type = "signed", powers = c(c(1:10), 
    seq(from = 12, to = 20, by = 2)), max.block.size = ncol(exp), 
    min.module.size = 30, merge.cut.height = 0.25, cor.heatmap.x.lab.angle = 45, 
    cor.heatmap.sig = "text", to.cytoscape = F, plot = T) 
{
    if (!all(rownames(exp) == rownames(pheno))) {
        stop("The expression and phenotype data must have the same row names!")
    }
    file_conn <- file("WGCNA-details.txt", open = "w")
    cat("+++ Checking genes and samples across this dataset...")
    cat("\n")
    cat("+++ Checking genes and samples across this dataset...", 
        file = file_conn, sep = "\n")
    gsg <- goodSamplesGenes(exp, verbose = FALSE)
    if (!gsg$allOK) {
        if (sum(!gsg$goodGenes) > 0) {
            cat(paste("+++ Removing genes:", paste(names(exp)[!gsg$goodGenes], 
                collapse = ", ")))
            cat("\n")
            cat(paste("+++ Removing genes:", paste(names(exp)[!gsg$goodGenes], 
                collapse = ", ")), file = file_conn)
            cat("\n", file = file_conn, sep = "\n")
        }
        if (sum(!gsg$goodSamples) > 0) {
            cat(paste("+++ Removing samples:", paste(rownames(exp)[!gsg$goodSamples], 
                collapse = ", ")))
            cat("\n")
            cat(paste("+++ Removing samples:", paste(rownames(exp)[!gsg$goodSamples], 
                collapse = ", ")), file = file_conn)
            cat("\n", file = file_conn, sep = "\n")
        }
        exp <- exp[gsg$goodSamples, gsg$goodGenes]
    }
    else {
        cat("+++ Genes and samples are good!")
        cat("\n")
        cat("+++ Genes and samples are good!", file = file_conn, 
            sep = "\n")
        cat("\n", file = file_conn, sep = "\n")
    }
    sampleTree <- stats::hclust(stats::dist(exp), method = "complete")
    par(cex = 0.6)
    par(mar = c(1, 5, 2, 0))
    cat("+++ Detecting sample outliers...")
    cat("\n")
    cat("+++ Detecting sample outliers...", file = file_conn, 
        sep = "\n")
    plot(sampleTree, main = "Sample clustering to detect outliers", 
        sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    height_cutoff <- as.numeric(readline(prompt = "+++ Height cutoff value: "))
    cat(paste0("+++ Height cutoff value: ", height_cutoff), file = file_conn, 
        sep = "\n")
    abline(h = height_cutoff, col = "red")
    clust <- cutreeStatic(sampleTree, cutHeight = height_cutoff, 
        minSize = 1)
    keepSamples <- (clust == 1)
    cat(paste0("+++ Based on the height cutoff you inputted, ", 
        sum(clust != 1), " outliers were removed from this dataset!"))
    cat("\n")
    cat(paste0("+++ Based on the height cutoff you inputted, ", 
        sum(clust != 1), " outliers were removed from this dataset!"), 
        file = file_conn, sep = "\n")
    cat("\n", file = file_conn, sep = "\n")
    exp2 <- exp[keepSamples, ]
    pheno2 <- pheno[rownames(exp2), ]
    sampleTree2 <- stats::hclust(stats::dist(exp2), method = "complete")
    plot(sampleTree2, main = "New sample clustering after remove outliers", 
        sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    traitColors <- numbers2colors(pheno2, signed = FALSE)
    plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(pheno2), 
        main = "Sample dendrogram and trait heatmap")
    cat("+++ Picking soft threshold to construct scale-free network...")
    cat("\n")
    cat("+++ Picking soft threshold to construct scale-free network...", 
        file = file_conn, sep = "\n")
    sft <- suppressWarnings(pickSoftThreshold(exp2, powerVector = powers, 
        networkType = net.type, verbose = FALSE))
    cat("\n")
    tmp <- sft$fitIndices
    power_s <- Inf
    while (power_s > max(powers)) {
        R.squre_cutoff <- as.numeric(readline(prompt = "+++ R squre cutoff value (0.85): "))
        power_s <- min(tmp$Power[tmp$SFT.R.sq >= R.squre_cutoff])
    }
    cat(paste0("+++ R squre cutoff value: ", R.squre_cutoff), 
        file = file_conn, sep = "\n")
    cat(paste0("+++ Based on the R squre cutoff you inputted, ", 
        power_s, " is the most appropriate power for constructing scale-free network!"), 
        file = file_conn, sep = "\n")
    cat(paste0("+++ Based on the R squre cutoff you inputted, ", 
        power_s, " is the most appropriate power!"))
    cat("\n")
    tmp$V1 <- ifelse(tmp$Power == power_s, tmp$SFT.R.sq, NA)
    tmp$V2 <- ifelse(tmp$Power == power_s, tmp$mean.k., NA)
    p1 <- ggplot(tmp) + geom_hline(yintercept = R.squre_cutoff, 
        lty = 2, col = "grey20", lwd = 0.65) + geom_point(aes(Power, 
        V1), data = na.omit(tmp), size = 9, shape = 22, alpha = 0.8, 
        fill = "#EEDF7C", color = "#EEDF7C") + geom_point(aes(Power, 
        SFT.R.sq), size = 2, shape = 21, fill = "#ED6355", color = "#ED6355") + 
        geom_point(aes(Power, SFT.R.sq), size = 4, shape = 21, 
            stroke = 1) + labs(x = "Soft threshold (Power)", 
        y = "Topology model fitness R^2", title = "Scale independence") + 
        theme_bw(base_rect_size = 1.5) + theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5, colour = "black", 
            face = "bold"), axis.text.x = element_text(size = 10, 
            colour = "black"), axis.text.y = element_text(size = 10, 
            colour = "black"), axis.title.x = element_text(size = 13, 
            colour = "black", face = "bold"), axis.title.y = element_text(size = 13, 
            colour = "black", face = "bold"))
    p2 <- ggplot(tmp) + geom_point(aes(Power, V2), data = na.omit(tmp), 
        size = 9, shape = 22, alpha = 0.8, fill = "#EEDF7C", 
        color = "#EEDF7C") + geom_point(aes(Power, mean.k.), 
        size = 2, shape = 21, fill = "#3E94B5", color = "#3E94B5") + 
        geom_point(aes(Power, mean.k.), size = 4, shape = 21, 
            stroke = 1) + labs(x = "Soft threshold (Power)", 
        y = "Mean connectivity (Node degree)", title = "Scale-free distribuion") + 
        theme_bw(base_rect_size = 1.5) + theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5, colour = "black", 
            face = "bold"), axis.text.x = element_text(size = 10, 
            colour = "black"), axis.text.y = element_text(size = 10, 
            colour = "black"), axis.title.x = element_text(size = 13, 
            colour = "black", face = "bold"), axis.title.y = element_text(size = 13, 
            colour = "black", face = "bold"))
    print(p1 + p2)
    if (plot) {
        ggsave(filename = "Pick-softThreshold.pdf", width = 7.4, 
            height = 3.9)
    }
    cat("+++ Calculating network...")
    cat("\n")
    net <- suppressWarnings(blockwiseModules(exp2, power = power_s, 
        maxBlockSize = max.block.size, corType = "pearson", networkType = net.type, 
        TOMType = net.type, minModuleSize = min.module.size, 
        mergeCutHeight = merge.cut.height, numericLabels = TRUE, 
        saveTOMs = TRUE, saveTOMFileBase = "Tom", verbose = FALSE))
    moduleColors <- labels2colors(net$colors)
    table(moduleColors)
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
        "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
        guideHang = 0.05)
    if (plot) {
        dev.copy2pdf(file = "Cluster-Dendrogram.pdf", width = 5, 
            height = 3)
    }
    MEs <- moduleEigengenes(exp2, colors = moduleColors)$eigengenes
    plotEigengeneNetworks(orderMEs(MEs), "Eigengene adjacency heatmap", 
        marDendro = c(3, 3, 2, 4), marHeatmap = c(3, 4, 2, 2), 
        plotDendrograms = T, xLabelsAngle = 90)
    if (to.cytoscape) {
        cat("+++ Exporting network to Cytoscape...")
        cat("\n")
        load(net$TOMFiles[1], verbose = F)
        TOM <- as.matrix(TOM)
        dissTOM <- 1 - TOM
        genes <- names(net$colors[net$blockGenes[[1]]])
        dimnames(TOM) <- list(genes, genes)
        cyt <- exportNetworkToCytoscape(TOM, edgeFile = "cyto.edges.txt", 
            nodeFile = "cyto.nodes.txt", weighted = TRUE, threshold = 0, 
            nodeNames = genes, nodeAttr = moduleColors[net$blockGenes[[1]]])
    }
    cat("+++ Extracting module genes...")
    cat("\n")
    moduleGenes <- sapply(unique(moduleColors), function(x) {
        names(net$colors)[which(moduleColors == x)]
    })
    cat("+++ Calculating module and trait relationships...")
    cat("\n")
    moduleTraitCor <- cor(MEs, pheno2, use = "p")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(pheno2))
    if (cor.heatmap.sig == "text") {
        textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", 
            signif(moduleTraitPvalue, 1), ")", sep = "")
        dim(textMatrix) <- dim(moduleTraitCor)
        annM <- textMatrix
    }
    if (cor.heatmap.sig == "signif") {
        sig <- apply(moduleTraitPvalue, 2, function(x) {
            ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 
                0.05, "*", "ns")))
        })
        sigMatrix <- paste(signif(moduleTraitCor, 2), "\n", sig, 
            sep = "")
        dim(sigMatrix) <- dim(moduleTraitCor)
        annM <- sigMatrix
    }
    if (cor.heatmap.sig == "none") {
        annM <- moduleTraitCor
        annM[annM < 2] <- ""
    }
    colors <- substring(rownames(moduleTraitCor), 3)
    names(colors) <- rownames(moduleTraitCor)
    left <- ComplexHeatmap::HeatmapAnnotation(which = "row", 
        Cluster = rownames(moduleTraitCor), border = FALSE, show_legend = FALSE, 
        show_annotation_name = FALSE, simple_anno_size = grid::unit(4, 
            "mm"), col = list(Cluster = colors))
    ch <- ComplexHeatmap::Heatmap(moduleTraitCor, name = " ", 
        border = TRUE, col = circlize::colorRamp2(seq(from = -1, 
            to = 1, length.out = 50), blueWhiteRed(50)), show_row_names = TRUE, 
        row_names_side = "right", row_names_gp = grid::gpar(fontsize = 12), 
        show_column_names = TRUE, column_names_side = "bottom", 
        column_names_gp = grid::gpar(fontsize = 12), column_names_rot = cor.heatmap.x.lab.angle, 
        column_names_centered = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, 
        left_annotation = left, color_space = "RGB", heatmap_legend_param = list(border = TRUE, 
            legend_height = unit(5, "cm"), labels_gp = grid::gpar(fontsize = 11)), 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.text(annM[i, j], x, y, gp = grid::gpar(fontsize = 10))
        })
    print(ch)
    if (plot) {
        dev.copy2pdf(file = "cor-heatmap.pdf", width = 5, height = 4)
    }
    cat("+++ Calculating module membership and gene significance...")
    cat("\n")
    geneModuleMembership <- cor(exp2, orderMEs(MEs), use = "p")
    geneSignificanceCor <- cor(exp2, pheno2, use = "p")
    close(file_conn)
    cat("+++ Done!")
    cat("\n")
    return(list(soft.pick.res = sft, module.gene.list = moduleGenes, 
        GS = geneSignificanceCor, MM = geneModuleMembership, 
        plots = list(soft.pick = p1 + p2, cor.heatmap = ch)))
}
