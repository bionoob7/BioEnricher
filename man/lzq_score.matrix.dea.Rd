% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lzq_score.matrix.dea.R
\name{lzq_score.matrix.dea}
\alias{lzq_score.matrix.dea}
\title{Perform differential analysis for the score matrix}
\usage{
lzq_score.matrix.dea(
  score.matrix,
  groups,
  control.group,
  Select.P = "FDR",
  cutoff.P = 0.05,
  cutoff.logFC = 2,
  ...
)
}
\arguments{
\item{score.matrix}{Matrix of scores with features as rows and samples as columns}

\item{groups}{Vector indicating group membership for each sample}

\item{control.group}{Name of the control group for comparison}

\item{Select.P}{Type of p-value to use ("FDR" or "pvalue"). Default is "FDR".}

\item{cutoff.P}{P-value cutoff for significance. Default is 0.05.}

\item{cutoff.logFC}{Log fold change cutoff for significance. Default is 2.}

\item{...}{Additional arguments passed to volcano plot function}
}
\value{
A list containing:
\itemize{
\item res - Data frame with differential analysis results
\item plot - Volcano plot visualization
}
}
\description{
Performs differential analysis on a score matrix (e.g., pathway scores, TF activities)
between different groups. Combines statistical testing with volcano plot visualization
to identify significantly different scores between conditions.
}
\examples{
\dontrun{
# Differential analysis of pathway scores
pathway_dea <- lzq_score.matrix.dea(score.matrix = pathway_scores,
                                    groups = sample_groups,
                                    control.group = "Normal")

# Access results and plot
results <- pathway_dea$res
volcano_plot <- pathway_dea$plot
}

}
