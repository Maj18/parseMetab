# parseMetab News & Release History

## parseMetab 1.0.3 (Release - February 2026)

🎉 **Major Release - Production Ready**

**DOI**: [10.5281/zenodo.18730007](https://doi.org/10.5281/zenodo.18730007)
**Improvements:**
- Updated getSystemNetwork, now we calculate correlation within each cancer types and then take the average
- Updated diffEffectBoxplot_bySystem, now all features are plotted.


## parseMetab 1.0.2 (Release - February 2026)

🎉 **Major Release - Production Ready**

**DOI**: [10.5281/zenodo.18647837](https://doi.org/10.5281/zenodo.18647837)

**Improvements:**
- Enhanced roxygen2 documentation block with:
  - Clearer function description explaining nested list structure
  - Detailed `@param` descriptions for `INDIR`, `pattern`, and `meta`
  - Comprehensive `@return` documentation
  - Example usage in `@examples`
  - Added `@export` and `@importFrom` tags (dplyr, magrittr)

- Added runtime validation:
  - Check that `INDIR` exists and is a valid directory
  - Validate `meta` parameter contains required columns (`cancer`, `datasetID`)
  
- Improved error handling:
  - Check for file existence before attempting to read
  - Enhanced warning messages when files are missing or unreadable
  - Better error reporting with condition messages

### 2. Added Documentation for `getCircularBarplot()` Function

---


## parseMetab 1.0.1 (Release - February 2026)

🎉 **Major Release - Production Ready**

**DOI**: [10.5281/zenodo.18615953](https://doi.org/10.5281/zenodo.18615953)

### New Features
- **Co-expression Network Analysis** (`getSystemNetwork()`) - Construct and visualize metabolic system-specific co-expression networks with hub gene identification
- **Survival Analysis Integration** - Kaplan-Meier and Cox proportional hazards regression for survival studies
- Enhanced network visualization with hub highlighting and connectivity-based node sizing
- Improved KEGG pathway integration and metabolic database access

### Improvements
- Comprehensive roxygen2 documentation for all functions
- Enhanced README with detailed quick-start examples and troubleshooting guide
- Professional logo and branding assets
- Standardized DESCRIPTION file following CRAN best practices
- Better error handling and data validation
- Improved multi-cohort analysis pipelines
- Enhanced visualization aesthetics and publication-readiness

### Documentation
- Complete function reference documentation
- Updated vignettes with practical examples
- Improved inline code comments and roxygen examples
- Professional package logo (SVG/PNG formats)

### Dependencies
- Added explicit version requirements for all imports
- Added development tools to Suggests (devtools, roxygen2, testthat, knitr, rmarkdown)

### Backward Compatibility
- All 0.1.0 functions maintained with full backward compatibility
- Minor parameter adjustments for consistency

---

## parseMetab 0.1.0 (Initial Release - January 2026)

Initial development version with core functionality:
- Data parsing and processing for CellFie outputs
- GSVA and CellFie-based metabolic activity scoring
- Limma-based differential testing (including paired designs)
- Multi-cohort analysis support for TCGA and similar datasets
- Comprehensive visualization suite:
  - Dotplots for multi-cohort comparisons
  - Volcano plots for differential results
  - Heatmaps for task scores
  - Boxplots stratified by cancer type and biological system
  - Circular barplots for sample size summaries
- Effect size and differential activity visualizations
- KEGG pathway integration
- Early co-expression network analysis
- Basic survival analysis capabilities
