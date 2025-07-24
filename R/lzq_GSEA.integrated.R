#' @export
lzq_GSEA.integrated <- function (genes, gene.type = "SYMBOL", organism = "Human", GO.ont = "BP", 
    KEGG.use.internal.data = F, perform.WikiPathways = F, perform.Reactome = F, 
    perform.MsigDB = F, MsigDB.category = "H", perform.disease.ontoloty = F, 
    perform.Cancer.Gene.Network = F, perform.DisGeNET = F, perform.CellMarker = F, 
    perform.CMAP = T, pvalue.cutoff = 0.05, padjust.method = "BH", 
    min.Geneset.Size = 10, max.Geneset.Size = 1000, CMAP.min.Geneset.Size = 3) 
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
    cat(paste0("+++ Performing GO-", GO.ont, " enrichment..."))
    cat("\n")
    ego <- suppressWarnings(suppressMessages(clusterProfiler::gseGO(geneList = genes, 
        OrgDb = OrgDb, keyType = "ENTREZID", ont = GO.ont, pvalueCutoff = pvalue.cutoff, 
        pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
        maxGSSize = max.Geneset.Size, eps = 1e-100, seed = T, 
        verbose = F)))
    ego <- DOSE::setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID")
    cat("+++ Symplifying GO results...")
    cat("\n")
    if (!is.null(ego)) {
        simple_ego <- clusterProfiler::simplify(ego)
    }
    else {
        simple_ego <- NULL
    }
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
    ekegg <- DOSE::setReadable(ekegg, OrgDb = OrgDb, keyType = "ENTREZID")
    cat("+++ Performing Module KEGG enrichment...")
    cat("\n")
    emkegg <- suppressWarnings(suppressMessages(clusterProfiler::gseMKEGG(geneList = genes, 
        organism = KEGG.organism, keyType = "kegg", minGSSize = min.Geneset.Size, 
        maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, 
        pAdjustMethod = padjust.method, eps = 1e-100, seed = T, 
        verbose = F)))
    emkegg <- DOSE::setReadable(emkegg, OrgDb = OrgDb, keyType = "ENTREZID")
    res <- list(GO = ego, simplyGO = simple_ego, KEGG = ekegg, 
        Module.KEGG = emkegg)
    if (perform.WikiPathways) {
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
        res[["WikiPathways"]] <- DOSE::setReadable(eWP, OrgDb = OrgDb, 
            keyType = "ENTREZID")
    }
    if (perform.Reactome) {
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
        res[["Reactome"]] <- DOSE::setReadable(eRP, OrgDb = OrgDb, 
            keyType = "ENTREZID")
    }
    if (perform.disease.ontoloty) {
        cat("+++ Performing Disease Ontoloty enrichment...")
        cat("\n")
        eDO <- suppressWarnings(DOSE::gseDO(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F))
        res[["DiseaseOntology"]] <- DOSE::setReadable(eDO, OrgDb = OrgDb, 
            keyType = "ENTREZID")
    }
    if (perform.Cancer.Gene.Network) {
        cat("+++ Performing Cancer Gene Network enrichment...")
        cat("\n")
        eNCG <- suppressWarnings(DOSE::gseNCG(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = T, verbose = F))
        res[["CancerGeneNetwork"]] <- DOSE::setReadable(eNCG, 
            OrgDb = OrgDb, keyType = "ENTREZID")
    }
    if (perform.DisGeNET) {
        cat("+++ Performing DisGeNET enrichment...")
        cat("\n")
        eDGN <- suppressWarnings(DOSE::gseDGN(geneList = genes, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
            eps = 1e-100, seed = TRUE, verbose = FALSE))
        res[["DisGeNET"]] <- DOSE::setReadable(eDGN, OrgDb = OrgDb, 
            keyType = "ENTREZID")
    }
    if (perform.CellMarker) {
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
            eps = 1e-100, seed = TRUE, verbose = FALSE, TERM2GENE = cells))
        res[["CellMarker.Enrichment"]] <- list(Enrich.results = DOSE::setReadable(eCM, 
            OrgDb = OrgDb, keyType = "ENTREZID"), Genesets = cells)
    }
    if (perform.MsigDB) {
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
        res[[paste0("MsigDB.", MsigDB.category)]] <- list(Enrich.results = DOSE::setReadable(eMSIG, 
            OrgDb = OrgDb, keyType = "ENTREZID"), Genesets = mg)
    }
    if (perform.CMAP) {
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
        res[["CMAP"]] <- list(Enrich.results = DOSE::setReadable(eCMAP, 
            OrgDb = OrgDb, keyType = "ENTREZID"), Genesets = CMAPfromDSEATM)
    }
    nsig <- ifelse(is.null(nrow(res$GO)), 0, nrow(res$GO)) + 
        ifelse(is.null(nrow(res$KEGG)), 0, nrow(res$KEGG)) + 
        ifelse(is.null(nrow(res$Module.KEGG)), 0, nrow(res$Module.KEGG)) + 
        ifelse(is.null(nrow(res$WikiPathways)), 0, nrow(res$WikiPathways)) + 
        ifelse(is.null(nrow(res$Reactome)), 0, nrow(res$Reactome)) + 
        ifelse(is.null(nrow(res$DiseaseOntology)), 0, nrow(res$DiseaseOntology)) + 
        ifelse(is.null(nrow(res$CancerGeneNetwork)), 0, nrow(res$CancerGeneNetwork)) + 
        ifelse(is.null(nrow(res$DisGeNET)), 0, nrow(res$DisGeNET)) + 
        ifelse(is.null(nrow(res$CellMarker.Enrichment$Enrich.results)), 
            0, nrow(res$CellMarker.Enrichment$Enrich.results)) + 
        ifelse(is.null(nrow(res[[paste0("MsigDB.", MsigDB.category)]]$Enrich.results)), 
            0, nrow(res[[paste0("MsigDB.", MsigDB.category)]]$Enrich.results)) + 
        ifelse(is.null(nrow(res$CMAP$Enrich.results)), 0, nrow(res$CMAP$Enrich.results))
    cat(paste0("+++ ", nsig, " significant terms were detected..."))
    cat("\n")
    cat("+++ Done!")
    return(res)
}
