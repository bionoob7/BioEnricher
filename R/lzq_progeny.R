#' @export
lzq_progeny <- function (exp, scale = T, organism = "Human", top = 100, perm = 1, 
    z_scores = F, get_nulldist = F, assay_name = "RNA", return_assay = F) 
{
    pathways <- progeny::progeny(expr = exp, scale = scale, organism = organism, 
        top = top, verbose = F, z_scores = z_scores, get_nulldist = get_nulldist, 
        perm = perm, assay_name = assay_name, return_assay = return_assay)
    return(pathways)
}
