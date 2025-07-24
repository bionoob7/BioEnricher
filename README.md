# BioEnricher

Integrate Analysis and Visualization for Bioinformatic Enrichment Analyzer

## Description

This package addresses two main issues:

1. **Seamless Integration for Enrichment Analysis**: Encompasses diverse functionalities including:
   - GO, KEGG, WikiPathways, Reactome
   - MsigDB, Disease Ontology, Cancer Gene Network
   - DisGeNET, CellMarker, and CMAP (drugs)
   - Transcription factor activity inference
   - PROGENy cancer pathways
   - Gene information, PubMed records, and GEO metadata search

2. **Advanced Visualization Functions**: Streamlines the process for faster and more convenient data presentation.

## Installation

You can install the development version from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install BioEnricher
devtools::install_github("Bionoob7/BioEnricher")
```

## Dependencies

This package requires R >= 4.3.0 and depends on several Bioconductor packages. Make sure you have Bioconductor installed:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install()
```

## Usage

```r
library(BioEnricher)

# Example usage will be added here
```

## Author

Zaoqu Liu (<liuzaoqu@163.com>)

## Maintainer

Keqiang Ma (<bionoob777@outlook.com>)

## License

MIT License

## Issues

Please report issues at: https://github.com/bionoob7/BioEnricher/issues