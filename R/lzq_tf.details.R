#' @export
lzq_tf.details <- function (tf.genes = NULL, targets = NULL, organism = "Human", 
    use.cancer.regulons = F, confidences = c("A", "B", "C")) 
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
    if (is.null(tf.genes)) {
        if (is.null(targets)) {
            stop("You must specify one or multiple TFs or targets!")
        }
        res <- as.data.frame(regulons[regulons$target %in% targets, 
            ])
    }
    if (is.null(targets)) {
        if (is.null(tf.genes)) {
            stop("You must specify one or multiple TFs or targets!")
        }
        res <- as.data.frame(regulons[regulons$tf %in% tf.genes, 
            ])
    }
    if (!is.null(targets) & !is.null(tf.genes)) {
        res <- as.data.frame(regulons[regulons$tf %in% tf.genes, 
            ])
        res <- res[res$target %in% targets, ]
    }
    return(res)
}
