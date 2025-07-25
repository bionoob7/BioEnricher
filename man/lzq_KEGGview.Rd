% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lzq_KEGGview.R
\name{lzq_KEGGview}
\alias{lzq_KEGGview}
\title{KEGG pathway visualization}
\usage{
lzq_KEGGview(
  gene.data = NULL,
  gene.type = "SYMBOL",
  pathway.id,
  species = "hsa",
  figure.suffix = ""
)
}
\arguments{
\item{gene.data}{Named numeric vector of gene expression values (fold changes or scores)}

\item{gene.type}{Type of gene identifiers in gene.data ("SYMBOL", "ENTREZID", etc.). Default is "SYMBOL".}

\item{pathway.id}{KEGG pathway identifier (e.g., "hsa04110" for cell cycle)}

\item{species}{KEGG species code (e.g., "hsa" for human, "mmu" for mouse). Default is "hsa".}

\item{figure.suffix}{Suffix to add to output figure filename. Default is "".}
}
\value{
A list containing:
\itemize{
\item plot.data.gene - Gene data mapped to pathway
\item image - PNG image of the pathway with expression data overlaid
}
Also saves pathway diagram as PNG file.
}
\description{
Visualizes gene expression data on KEGG pathway diagrams. Maps gene expression
values onto pathway components to show which genes are up- or down-regulated
in the context of biological pathways.
}
\examples{
\dontrun{
# Visualize cell cycle pathway with expression data
gene_expr <- c("TP53" = -2.5, "CCND1" = 1.8, "RB1" = -1.2)
kegg_plot <- lzq_KEGGview(gene.data = gene_expr, 
                         pathway.id = "hsa04110")

# Custom species and suffix
kegg_plot <- lzq_KEGGview(gene.data = gene_expr,
                         pathway.id = "mmu04110",
                         species = "mmu",
                         figure.suffix = "treatment_vs_control")
}

}
