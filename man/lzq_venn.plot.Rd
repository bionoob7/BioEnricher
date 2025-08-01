% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lzq_venn.plot.R
\name{lzq_venn.plot}
\alias{lzq_venn.plot}
\title{Plot venn diagram}
\usage{
lzq_venn.plot(
  data.list,
  set.names = NULL,
  fill.colors = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF"),
  fill.alpha = 0.5,
  border.color = "black",
  border.alpha = 1,
  border.width = 1,
  border.linetype = 1,
  set.name.colors = "black",
  set.name.size = 6,
  text.color = "black",
  text.size = 5.5
)
}
\arguments{
\item{data.list}{List of character vectors, each containing identifiers for one set}

\item{set.names}{Character vector of names for each set. If NULL, uses list names. Default is NULL.}

\item{fill.colors}{Vector of colors for filling each set. Default is c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF").}

\item{fill.alpha}{Transparency of fill colors (0-1). Default is 0.5.}

\item{border.color}{Color of set borders. Default is "black".}

\item{border.alpha}{Transparency of borders (0-1). Default is 1.}

\item{border.width}{Width of borders. Default is 1.}

\item{border.linetype}{Line type for borders. Default is 1.}

\item{set.name.colors}{Color of set names. Default is "black".}

\item{set.name.size}{Size of set name text. Default is 6.}

\item{text.color}{Color of overlap numbers. Default is "black".}

\item{text.size}{Size of overlap number text. Default is 5.5.}
}
\value{
A ggplot2 object representing the Venn diagram.
}
\description{
Creates Venn diagrams to visualize overlaps between multiple sets of genes
or other identifiers. Supports up to 4 sets with customizable colors and styling.
}
\examples{
\dontrun{
# Basic Venn diagram with 3 gene sets
gene_sets <- list(
  Set1 = c("TP53", "BRCA1", "EGFR", "MYC"),
  Set2 = c("TP53", "BRCA2", "EGFR", "RB1"),
  Set3 = c("MYC", "EGFR", "APC", "KRAS")
)
venn_plot <- lzq_venn.plot(data.list = gene_sets)

# Customized Venn diagram
venn_plot <- lzq_venn.plot(data.list = gene_sets,
                          set.names = c("Tumor", "Normal", "Metastasis"),
                          fill.colors = c("red", "blue", "green"))
}

}
