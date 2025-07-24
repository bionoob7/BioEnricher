#' @export
lzq_ORA.integrated <- function (genes, background.genes = NULL, gene.type = "SYMBOL", 
    organism = "Human", GO.ont = "BP", KEGG.use.internal.data = F, 
    perform.WikiPathways = F, perform.Reactome = F, perform.MsigDB = F, 
    MsigDB.category = "H", perform.disease.ontoloty = F, perform.Cancer.Gene.Network = F, 
    perform.DisGeNET = F, perform.CellMarker = F, perform.CMAP = T, 
    pvalue.cutoff = 0.05, qvalue.cutoff = 0.05, padjust.method = "BH", 
    min.Geneset.Size = 10, max.Geneset.Size = 1000, CMAP.min.Geneset.Size = 3) 
{
    if (gene.type == "SYMBOL") {
        cat("+++ Updating gene symbols...")
        cat("\n")
        genes <- suppressWarnings(lzq_updateSymbol(genes = genes, 
            unmapGene_keep = T)[[2]])
        if (!is.null(background.genes)) {
            background.genes <- suppressWarnings(lzq_updateSymbol(genes = stats::na.omit(background.genes), 
                unmapGene_keep = T)[[2]])
        }
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
        genes <- suppressWarnings(clusterProfiler::bitr(genes, 
            fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))[, 
            "ENTREZID"]
        if (!is.null(background.genes)) {
            background.genes <- suppressWarnings(clusterProfiler::bitr(background.genes, 
                fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))[, 
                "ENTREZID"]
        }
    }
    cat(paste0("+++ Performing GO-", GO.ont, " enrichment..."))
    cat("\n")
    ego <- suppressWarnings(suppressMessages(clusterProfiler::enrichGO(gene = genes, 
        OrgDb = OrgDb, keyType = "ENTREZID", ont = GO.ont, pvalueCutoff = pvalue.cutoff, 
        qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
        minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
        readable = T, universe = background.genes)))
    if (!is.null(ego)) {
        ego <- lzq_getEF(ego)
    }
    cat("+++ Symplifying GO results...")
    cat("\n")
    if (!is.null(ego)) {
        simple_ego <- clusterProfiler::simplify(ego)
        if (!is.null(simple_ego)) {
            simple_ego <- lzq_getEF(simple_ego)
        }
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
    ekegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichKEGG(gene = genes, 
        organism = KEGG.organism, pvalueCutoff = pvalue.cutoff, 
        qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
        minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
        use_internal_data = KEGG.use.internal.data, universe = background.genes)))
    ekegg <- DOSE::setReadable(ekegg, OrgDb = OrgDb, keyType = "ENTREZID")
    if (!is.null(ekegg)) {
        ekegg <- lzq_getEF(ekegg)
    }
    cat("+++ Performing Module KEGG enrichment...")
    cat("\n")
    emkegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichMKEGG(gene = genes, 
        organism = KEGG.organism, pvalueCutoff = pvalue.cutoff, 
        qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
        minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
        universe = background.genes)))
    emkegg <- DOSE::setReadable(emkegg, OrgDb = OrgDb, keyType = "ENTREZID")
    if (!is.null(emkegg)) {
        emkegg <- lzq_getEF(emkegg)
    }
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
        eWP <- suppressWarnings(clusterProfiler::enrichWP(gene = genes, 
            organism = WikiPathways.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            universe = background.genes))
        eWP <- DOSE::setReadable(eWP, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eWP)) {
            res[["WikiPathways"]] <- lzq_getEF(eWP)
        }
        else {
            res[["WikiPathways"]] <- eWP
        }
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
        eRP <- suppressWarnings(ReactomePA::enrichPathway(gene = genes, 
            organism = Reactome.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            readable = T, universe = background.genes))
        eRP <- DOSE::setReadable(eRP, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eRP)) {
            res[["Reactome"]] <- lzq_getEF(eRP)
        }
        else {
            res[["Reactome"]] <- eRP
        }
    }
    if (perform.disease.ontoloty) {
        cat("+++ Performing Disease Ontoloty enrichment...")
        cat("\n")
        eDO <- suppressWarnings(DOSE::enrichDO(gene = genes, 
            ont = "DO", pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, readable = T, universe = background.genes))
        eDO <- DOSE::setReadable(eDO, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eDO)) {
            res[["DiseaseOntology"]] <- lzq_getEF(eDO)
        }
        else {
            res[["DiseaseOntology"]] <- eDO
        }
    }
    if (perform.Cancer.Gene.Network) {
        cat("+++ Performing Cancer Gene Network enrichment...")
        cat("\n")
        eNCG <- suppressWarnings(DOSE::enrichNCG(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, readable = T, universe = background.genes))
        eNCG <- DOSE::setReadable(eNCG, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eNCG)) {
            res[["CancerGeneNetwork"]] <- lzq_getEF(eNCG)
        }
        else {
            res[["CancerGeneNetwork"]] <- eNCG
        }
    }
    if (perform.DisGeNET) {
        cat("+++ Performing DisGeNET enrichment...")
        cat("\n")
        eDGN <- suppressWarnings(DOSE::enrichDGN(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, readable = T, universe = background.genes))
        eDGN <- DOSE::setReadable(eDGN, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eDGN)) {
            res[["DisGeNET"]] <- lzq_getEF(eDGN)
        }
        else {
            res[["DisGeNET"]] <- eDGN
        }
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
        eCM <- suppressWarnings(clusterProfiler::enricher(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, TERM2GENE = cells, 
            universe = background.genes))
        eCM <- DOSE::setReadable(eCM, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eCM)) {
            res[["CellMarker.Enrichment"]] <- list(Enrich.results = lzq_getEF(eCM), 
                Genesets = cells)
        }
        else {
            res[["CellMarker.Enrichment"]] <- list(Enrich.results = eCM, 
                Genesets = cells)
        }
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
        eMSIG <- suppressWarnings(clusterProfiler::enricher(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, TERM2GENE = mg, universe = background.genes))
        eMSIG <- DOSE::setReadable(eMSIG, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eMSIG)) {
            res[[paste0("MsigDB.", MsigDB.category)]] <- list(Enrich.results = lzq_getEF(eMSIG), 
                Genesets = mg)
        }
        else {
            res[[paste0("MsigDB.", MsigDB.category)]] <- list(Enrich.results = eMSIG, 
                Genesets = mg)
        }
    }
    if (perform.CMAP) {
        if (organism != "Human") {
            stop("CMAP only supports organism = Human!")
        }
        cat("+++ Performing CMAP enrichment...")
        cat("\n")
        eCMAP <- suppressWarnings(clusterProfiler::enricher(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = CMAP.min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, TERM2GENE = CMAPfromDSEATM[, 
                seq_len(2)], universe = background.genes))
        eCMAP <- DOSE::setReadable(eCMAP, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eCMAP)) {
            res[["CMAP"]] <- list(Enrich.results = lzq_getEF(eCMAP), 
                Genesets = CMAPfromDSEATM)
        }
        else {
            res[["CMAP"]] <- list(Enrich.results = eCMAP, Genesets = CMAPfromDSEATM)
        }
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
