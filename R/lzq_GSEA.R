#' @export
lzq_GSEA <- function (genes, gene.type = "SYMBOL", enrich.type, organism = "Human", 
    GO.ont = "BP", GO.simplify = T, KEGG.use.internal.data = F, 
    MsigDB.category = "H", CMAP.min.Geneset.Size = 3, pvalue.cutoff = 0.05, 
    padjust.method = "BH", min.Geneset.Size = 10, max.Geneset.Size = 1000) 
{
    d <- as.data.frame(genes)
    if (gene.type == "SYMBOL") {
        cat("+++ Updating gene symbols...")
        cat("\n")
        d2 <- suppressWarnings(lzq_updateSymbolforDL(d, unmapGene_keep = T)) %>% 
            dplyr::arrange(dplyr::desc(genes))
    }
    if (organism == "Human") {
        OrgDb <- "org.Hs.eg.db"
    }
    if (organism == "Mouse") {
        OrgDb <- "org.Mm.eg.db"
    }
    if (gene.type != "ENTREZID") {
        cat(paste0("+++ Transforming ", gene.type, " to ENTREZID..."))
        cat("\n")
        d3 <- suppressWarnings(clusterProfiler::bitr(rownames(d2), 
            fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))
        d3 <- merge(d3, d2, by.x = 1, by.y = 0) %>% dplyr::arrange(dplyr::desc(genes))
        genes <- d3$genes
        names(genes) <- d3$ENTREZID
    }
    else {
        genes <- d2$genes
        names(genes) <- rownames(d2)
    }
    if (enrich.type == "GO") {
        cat(paste0("+++ Performing GO-", GO.ont, " enrichment..."))
        cat("\n")
        ego <- suppressWarnings(suppressMessages(clusterProfiler::gseGO(geneList = genes, 
            OrgDb = OrgDb, keyType = "ENTREZID", ont = GO.ont, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            eps = 1e-100, seed = T, verbose = F)))
        ego <- DOSE::setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(ego), " significant terms were detected..."))
        cat("\n")
        if (GO.simplify & !is.null(ego)) {
            cat("+++ Symplifying GO results...")
            cat("\n")
            simple_ego <- clusterProfiler::simplify(ego)
            res <- list(GO = ego, simplyGO = simple_ego)
        }
        else {
            res <- ego
        }
    }
    if (enrich.type == "KEGG") {
        cat("+++ Performing KEGG enrichment...")
        cat("\n")
        if (organism == "Human") {
            KEGG.organism <- "hsa"
        }
        if (organism == "Mouse") {
            KEGG.organism <- "mmu"
        }
        ekegg <- suppressWarnings(suppressMessages(clusterProfiler::gseKEGG(geneList = genes, 
            organism = KEGG.organism, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, 
            pAdjustMethod = padjust.method, eps = 1e-100, seed = T, 
            verbose = F, use_internal_data = KEGG.use.internal.data)))
        res <- ekegg
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "MKEGG") {
        cat("+++ Performing Module KEGG enrichment...")
        cat("\n")
        if (organism == "Human") {
            KEGG.organism <- "hsa"
        }
        if (organism == "Mouse") {
            KEGG.organism <- "mmu"
        }
        emkegg <- suppressWarnings(suppressMessages(clusterProfiler::gseMKEGG(geneList = genes, 
            organism = KEGG.organism, keyType = "kegg", minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, 
            pAdjustMethod = padjust.method, eps = 1e-100, seed = T, 
            verbose = F)))
        res <- emkegg
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "WikiPathways" | enrich.type == "WP") {
        cat("+++ Performing WikiPathways enrichment...")
        cat("\n")
        if (organism == "Human") {
            WikiPathways.organism <- "Homo sapiens"
        }
        if (organism == "Mouse") {
            WikiPathways.organism <- "Mus musculus"
        }
        eWP <- suppressWarnings(clusterProfiler::gseWP(geneList = genes, 
            organism = WikiPathways.organism, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, 
            pAdjustMethod = padjust.method, eps = 1e-100, seed = T, 
            verbose = F))
        res <- eWP
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "Reactome" | enrich.type == "RP") {
        cat("+++ Performing Reactome pathways enrichment...")
        cat("\n")
        if (organism == "Human") {
            Reactome.organism <- "human"
        }
        if (organism == "Mouse") {
            Reactome.organism <- "mouse"
        }
        eRP <- suppressWarnings(ReactomePA::gsePathway(geneList = genes, 
            organism = Reactome.organism, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, 
            pAdjustMethod = padjust.method, eps = 1e-100, seed = T, 
            verbose = F))
        res <- eRP
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "DO") {
        cat("+++ Performing Disease Ontoloty enrichment...")
        cat("\n")
        eDO <- suppressWarnings(DOSE::gseDO(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F))
        res <- eDO
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "CGN") {
        cat("+++ Performing Cancer Gene Network enrichment...")
        cat("\n")
        eNCG <- suppressWarnings(DOSE::gseNCG(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F))
        res <- eNCG
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "DisGeNET") {
        cat("+++ Performing DisGeNET enrichment...")
        cat("\n")
        eDGN <- suppressWarnings(DOSE::gseDGN(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F))
        res <- eDGN
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "CellMarker" | enrich.type == "CM") {
        cat("+++ Performing CellMarker enrichment...")
        cat("\n")
        if (organism == "Human") {
            cell_marker_data <- suppressMessages(vroom::vroom("http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt"))
        }
        if (organism == "Mouse") {
            cell_marker_data <- suppressMessages(vroom::vroom("http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt"))
        }
        cells <- cell_marker_data %>% dplyr::select(cellName, 
            geneID) %>% dplyr::mutate(geneID = strsplit(geneID, 
            ", ")) %>% tidyr::unnest(cols = geneID)
        eCM <- suppressWarnings(clusterProfiler::GSEA(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F, TERM2GENE = cells))
        res <- eCM
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "MsigDB") {
        cat(paste0("+++ Performing MsigDB-", MsigDB.category, 
            " enrichment..."))
        cat("\n")
        if (organism == "Human") {
            MsigDB.organism <- "Homo sapiens"
        }
        if (organism == "Mouse") {
            MsigDB.organism <- "Mus musculus"
        }
        mall <- msigdbr::msigdbr(species = MsigDB.organism)
        if (MsigDB.category == "All" | MsigDB.category == "ALL") {
            mg <- mall %>% dplyr::select(gs_name, entrez_gene)
        }
        else {
            mg <- msigdbr::msigdbr(species = MsigDB.organism, 
                category = MsigDB.category) %>% dplyr::select(gs_name, 
                entrez_gene)
        }
        eMSIG <- suppressWarnings(clusterProfiler::GSEA(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F, TERM2GENE = mg))
        res <- eMSIG
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "CMAP") {
        if (organism != "Human") {
            stop("CMAP only supports organism = Human!")
        }
        cat("+++ Performing CMAP enrichment...")
        cat("\n")
        eCMAP <- suppressWarnings(clusterProfiler::GSEA(geneList = genes, 
            minGSSize = CMAP.min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F, TERM2GENE = CMAPfromDSEATM[, 
                seq_len(2)]))
        res <- eCMAP
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    cat("+++ Done!")
    return(res)
}
