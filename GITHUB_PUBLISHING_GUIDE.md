# GitHub Publishing Guide for BioEnricher

## Prerequisites

Before publishing to GitHub, ensure you have:
- ✅ Git installed on your system
- ✅ GitHub account created
- ✅ Package passes local tests

## Step-by-Step Publishing Process

### Phase 1: Test Your Package

1. **Run the quick test**:
   ```r
   source("quick_test.R")
   ```

2. **Run comprehensive tests**:
   ```r
   source("test_and_publish.R")
   ```

3. **Fix any errors** before proceeding to GitHub

### Phase 2: Prepare for GitHub

1. **Update DESCRIPTION file**:
   - Replace `YOUR-USERNAME` with your actual GitHub username
   - Update the URL and BugReports fields:
   ```
   URL: https://github.com/YOUR-USERNAME/BioEnricher
   BugReports: https://github.com/YOUR-USERNAME/BioEnricher/issues
   ```

2. **Verify essential files exist**:
   - ✅ DESCRIPTION
   - ✅ NAMESPACE  
   - ✅ LICENSE
   - ✅ README.md
   - ✅ .gitignore
   - ✅ R/ directory with function files
   - ✅ man/ directory with documentation

### Phase 3: Initialize Git Repository

Open terminal/command prompt in your package directory and run:

```bash
# Initialize git repository
git init

# Add all files
git add .

# Make initial commit
git commit -m "Initial commit: BioEnricher R package with comprehensive documentation"
```

### Phase 4: Create GitHub Repository

1. Go to [GitHub.com](https://github.com)
2. Click the **"New"** button (green button) or **"+"** → **"New repository"**
3. Fill in repository details:
   - **Repository name**: `BioEnricher`
   - **Description**: `Integrate Analysis and Visualization for Bioinformatic Enrichment Analyzer`
   - **Visibility**: Public (recommended for R packages)
   - **DO NOT** initialize with README, .gitignore, or license (you already have these)
4. Click **"Create repository"**

### Phase 5: Connect Local Repository to GitHub

After creating the GitHub repository, run these commands:

```bash
# Set main branch
git branch -M main

# Add GitHub remote (replace YOUR-USERNAME with your actual username)
git remote add origin https://github.com/YOUR-USERNAME/BioEnricher.git

# Push to GitHub
git push -u origin main
```

### Phase 6: Verify Publication

1. **Check GitHub repository**: Visit your repository URL to ensure all files uploaded correctly

2. **Test installation from GitHub**:
   ```r
   # In a fresh R session
   devtools::install_github("YOUR-USERNAME/BioEnricher")
   library(BioEnricher)
   
   # Test basic functionality
   listEnrichMethod()
   ?gene.info
   ```

## Installation Instructions for Users

Once published, users can install your package with:

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install BioEnricher from GitHub
devtools::install_github("YOUR-USERNAME/BioEnricher")

# Load the package
library(BioEnricher)

# View available enrichment methods
listEnrichMethod()
```

## Troubleshooting

### Common Issues:

1. **Authentication Error**: 
   - Use GitHub Personal Access Token instead of password
   - Or use SSH key authentication

2. **Package Check Errors**:
   - Run `devtools::check()` and fix all errors before publishing
   - Warnings are acceptable but errors must be fixed

3. **Documentation Missing**:
   - Run `devtools::document()` to generate man/ files
   - Ensure all functions have proper roxygen documentation

4. **Installation Fails**:
   - Check that all dependencies are listed in DESCRIPTION
   - Verify NAMESPACE file is correct

### Getting Help:

- Check R package development guide: https://r-pkgs.org/
- GitHub documentation: https://docs.github.com/
- Ask questions on Stack Overflow with tags: [r], [r-package]

## Success! 🎉

Once published, your BioEnricher package will be:
- ✅ Publicly available on GitHub
- ✅ Installable via `devtools::install_github()`
- ✅ Fully documented with help files
- ✅ Ready for community use and contributions

Remember to:
- Update version numbers for new releases
- Tag releases on GitHub for version tracking
- Consider submitting to CRAN for wider distribution