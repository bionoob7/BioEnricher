% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lzq_quiet.R
\name{lzq_quiet}
\alias{lzq_quiet}
\title{Makes the function run without Verbose}
\usage{
lzq_quiet(..., messages = FALSE, cat = FALSE)
}
\arguments{
\item{...}{Expression or function call to execute quietly}

\item{messages}{Logical, whether to suppress messages. Default is FALSE.}

\item{cat}{Logical, whether to suppress cat() output. Default is FALSE.}
}
\value{
The result of the evaluated expression, with output suppressed
}
\description{
Suppresses output and messages from functions to run them quietly.
Useful for reducing console clutter when running functions that produce
verbose output or status messages.
}
\examples{
\dontrun{
# Run a function quietly
result <- lzq_quiet(some_verbose_function(data))

# Suppress both output and messages
result <- lzq_quiet(function_call(), messages = TRUE, cat = TRUE)
}

}
