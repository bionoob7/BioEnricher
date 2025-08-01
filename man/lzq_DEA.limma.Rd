% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lzq_DEA.limma.R
\name{lzq_DEA.limma}
\alias{lzq_DEA.limma}
\title{Perform differential expression analysis with limma}
\usage{
lzq_DEA.limma(
  expr,
  groups,
  contrasts = "Tumor-Normal",
  Select.P = "FDR",
  cutoff.P = 0.05,
  cutoff.logFC = 1
)
}
\arguments{
\item{expr}{Expression matrix with genes as rows and samples as columns}

\item{groups}{Vector indicating group membership for each sample}

\item{contrasts}{String specifying the contrast to test (e.g., "Tumor-Normal"). Default is "Tumor-Normal".}

\item{Select.P}{Type of p-value to use for filtering ("FDR" or "pvalue"). Default is "FDR".}

\item{cutoff.P}{P-value cutoff for significance. Default is 0.05.}

\item{cutoff.logFC}{Log fold change cutoff for significance. Default is 1.}
}
\value{
A data frame containing differential expression results with columns for
log fold change, p-values, adjusted p-values, and significance classification.
}
\description{
Conducts differential expression analysis using the limma package, which is
designed for analyzing gene expression data from microarray and RNA-seq experiments.
This function performs statistical testing to identify differentially expressed genes
between two conditions.
}
\examples{
\dontrun{
# Example expression matrix and groups
expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(expr_matrix) <- paste0("Gene", 1:100)
groups <- rep(c("Normal", "Tumor"), each = 5)

# Perform differential expression analysis
deg_results <- lzq_DEA.limma(expr = expr_matrix, 
                            groups = groups,
                            contrasts = "Tumor-Normal")

# With custom cutoffs
deg_strict <- lzq_DEA.limma(expr = expr_matrix,
                           groups = groups,
                           cutoff.P = 0.01,
                           cutoff.logFC = 1.5)
}

}
