#' @export
listEnrichMethod <- function () 
{
    print(c("GO", "KEGG", "MKEGG", "WikiPathways", "Reactome", 
        "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP"))
}
