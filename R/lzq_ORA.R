#' Over-representative analysis
#'
#' Performs over-representation analysis (ORA) for gene sets using various databases
#' including GO, KEGG, WikiPathways, Reactome, MsigDB, and others.
#'
#' @param genes Vector of genes for enrichment analysis
#' @param background.genes Vector of background genes. If NULL, uses all genes in the database. Default is NULL.
#' @param gene.type Type of gene identifiers ("SYMBOL", "ENTREZID", etc.). Default is "SYMBOL".
#' @param enrich.type Type of enrichment analysis to perform (e.g., "GO", "KEGG", "WikiPathways")
#' @param organism Organism for analysis ("Human", "Mouse", etc.). Default is "Human".
#' @param GO.ont GO ontology type ("BP", "CC", "MF"). Default is "BP".
#' @param GO.simplify Logical, whether to simplify GO results. Default is TRUE.
#' @param KEGG.use.internal.data Logical, whether to use internal KEGG data. Default is FALSE.
#' @param MsigDB.category MsigDB category to use. Default is "H".
#' @param CMAP.min.Geneset.Size Minimum gene set size for CMAP analysis. Default is 3.
#' @param pvalue.cutoff P-value cutoff for significance. Default is 0.05.
#' @param qvalue.cutoff Q-value cutoff for significance. Default is 0.05.
#' @param padjust.method Method for p-value adjustment. Default is "BH".
#' @param min.Geneset.Size Minimum gene set size. Default is 10.
#' @param max.Geneset.Size Maximum gene set size. Default is 1000.
#'
#' @return An enrichResult object containing the over-representation analysis results.
#'
#' @examples
#' \dontrun{
#' # GO enrichment analysis
#' ora_results <- lzq_ORA(genes = my_genes, enrich.type = "GO")
#' 
#' # KEGG pathway analysis
#' kegg_results <- lzq_ORA(genes = my_genes, enrich.type = "KEGG")
#' 
#' # With custom background
#' ora_custom <- lzq_ORA(genes = my_genes, 
#'                      background.genes = background_genes,
#'                      enrich.type = "GO")
#' }
#'
#' @export
lzq_ORA <- function (genes, background.genes = NULL, gene.type = "SYMBOL", 
    enrich.type, organism = "Human", GO.ont = "BP", GO.simplify = T, 
    KEGG.use.internal.data = F, MsigDB.category = "H", CMAP.min.Geneset.Size = 3, 
    pvalue.cutoff = 0.05, qvalue.cutoff = 0.05, padjust.method = "BH", 
    min.Geneset.Size = 10, max.Geneset.Size = 1000) 
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
    if (enrich.type == "GO") {
        cat(paste0("+++ Performing GO-", GO.ont, " enrichment..."))
        cat("\n")
        ego <- suppressWarnings(suppressMessages(clusterProfiler::enrichGO(gene = genes, 
            OrgDb = OrgDb, keyType = "ENTREZID", ont = GO.ont, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, universe = background.genes, 
            readable = T)))
        if (!is.null(ego)) {
            ego <- lzq_getEF(ego)
        }
        cat(paste0("+++ ", nrow(ego), " significant terms were detected..."))
        cat("\n")
        if (GO.simplify & !is.null(ego)) {
            cat("+++ Symplifying GO results...")
            cat("\n")
            simple_ego <- clusterProfiler::simplify(ego)
            if (!is.null(simple_ego)) {
                simple_ego <- lzq_getEF(simple_ego)
            }
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
        ekegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichKEGG(gene = genes, 
            organism = KEGG.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            universe = background.genes, use_internal_data = KEGG.use.internal.data)))
        res <- ekegg
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
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
        emkegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichMKEGG(gene = genes, 
            organism = KEGG.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            universe = background.genes)))
        res <- emkegg
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
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
        eWP <- suppressWarnings(clusterProfiler::enrichWP(gene = genes, 
            organism = WikiPathways.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            universe = background.genes))
        res <- eWP
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
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
        eRP <- suppressWarnings(ReactomePA::enrichPathway(gene = genes, 
            organism = Reactome.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
            minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            readable = T, universe = background.genes))
        res <- eRP
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "DO") {
        cat("+++ Performing Disease Ontoloty enrichment...")
        cat("\n")
        eDO <- suppressWarnings(DOSE::enrichDO(gene = genes, 
            ont = "DO", pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, readable = T, universe = background.genes))
        res <- eDO
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "CGN") {
        cat("+++ Performing Cancer Gene Network enrichment...")
        cat("\n")
        eNCG <- suppressWarnings(DOSE::enrichNCG(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, readable = T, universe = background.genes))
        res <- eNCG
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "DisGeNET") {
        cat("+++ Performing DisGeNET enrichment...")
        cat("\n")
        eDGN <- suppressWarnings(DOSE::enrichDGN(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, readable = T, universe = background.genes))
        res <- eDGN
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
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
        eCM <- suppressWarnings(clusterProfiler::enricher(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, TERM2GENE = cells, 
            universe = background.genes))
        res <- eCM
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
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
        eMSIG <- suppressWarnings(clusterProfiler::enricher(gene = genes, 
            pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, TERM2GENE = mg, universe = background.genes))
        res <- eMSIG
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type == "CMAP") {
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
        res <- eCMAP
        if (!is.null(res)) {
            res <- lzq_getEF(res)
        }
        cat(paste0("+++ ", nrow(res), " significant terms were detected..."))
        cat("\n")
    }
    if (enrich.type != "GO") {
        res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    }
    cat("+++ Done!")
    return(res)
}
