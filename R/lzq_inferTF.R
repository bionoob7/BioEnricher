#' @export
lzq_inferTF <- function (exp, organism = "Human", use.cancer.regulons = F, confidences = c("A", 
    "B", "C")) 
{
    if (organism == "Human") {
        if (use.cancer.regulons) {
            regulons <- dorothea::dorothea_hs
        }
        else {
            regulons <- dorothea::dorothea_hs_pancancer
        }
    }
    if (organism == "Mouse") {
        if (use.cancer.regulons) {
            regulons <- dorothea::dorothea_mm
        }
        else {
            regulons <- dorothea::dorothea_mm_pancancer
        }
    }
    regulons <- regulons %>% dplyr::filter(confidence %in% confidences)
    tf_activities <- suppressWarnings(dorothea::run_viper(input = exp, 
        regulons = regulons, options = list(method = "scale", 
            minsize = 4, eset.filter = F, cores = 1, verbose = F))) %>% 
        as.data.frame()
    return(tf_activities)
}
