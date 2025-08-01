% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lzq_getGR_BR.R
\name{lzq_getGR_BR}
\alias{lzq_getGR_BR}
\title{Get numeric GeneRatio and BgRatio from enrichResult}
\usage{
lzq_getGR_BR(res)
}
\arguments{
\item{res}{An enrichResult object from ORA or GSEA analysis}
}
\value{
The same enrichResult object with numeric GeneRatio and BgRatio columns
}
\description{
Converts the GeneRatio and BgRatio columns from character fractions (e.g., "5/100")
to numeric values (e.g., 0.05) in an enrichResult object. This is useful for
downstream analysis and plotting that requires numeric values.
}
\examples{
\dontrun{
# Convert ratios to numeric values
ora_results <- lzq_ORA(genes = my_genes, enrich.type = "GO")
ora_numeric <- lzq_getGR_BR(ora_results)
}

}
