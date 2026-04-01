# tessera <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/tessera)](https://CRAN.R-project.org/package=tessera)
[![R-CMD-check](https://github.com/jtonini/tessera/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jtonini/tessera/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**TESSERA** (Topology-Enhanced Similarity and Structure for Entity Regime
Analysis) detects regime boundaries in complex systems by analyzing how
similarity network topology shifts over time.

TESSERA's core insight: **similarity network topology carries predictive
information that direct classifiers miss.** By computing multi-measure
entity similarity, constructing networks, detecting communities, and
quantifying divergence across temporal bins, TESSERA reveals regime
boundaries (ecotones) that traditional feature-space methods overlook.

## Installation

```r
# From CRAN (when available)
install.packages("tessera")

# Development version from GitHub
# install.packages("pak")
pak::pak("jtonini/tessera")
```

## Quick Start

```r
library(tessera)

# Full pipeline with ecology preset
result <- tessera(
  coral_data,
  time_col = "year",
  preset = "ecology"
)

# Inspect results
result
plot(result)

# Step-by-step for full control
char <- tessera_characterize(data, time_col = "year", n_bins = 4)
sim  <- tessera_similarity(char, measures = c("bray_curtis", "simpson"))
net  <- tessera_network(sim, method = "knn", k = 12)
com  <- tessera_communities(net, algorithm = "leiden")
div  <- tessera_divergence(com, method = "nmi")
```

## Key Features

- **7 similarity measures**: Simpson, Cosine, Bray-Curtis, Jaccard,
  Pearson, Euclidean, Mahalanobis
- **Sparse by default**: kNN representation enables analysis of large
  datasets (10k+ entities) on standard hardware
- **Domain presets**: `"ecology"`, `"manufacturing"`, `"hpc"` — sensible
  defaults for each domain
- **Inspectable pipeline**: every step returns an S3 object with `print()`,
  `summary()`, and `plot()` methods
- **Quantile-based binning**: handles uneven temporal distributions
  (critical for long-span observational data)

## The TESSERA Framework

| Step | Function | What it does |
|------|----------|-------------|
| 1. Characterize | `tessera_characterize()` | Temporal binning, feature prep |
| 2a. Discover | `tessera_similarity()` | Multi-measure pairwise similarity |
| 2b. Discover | `tessera_network()` | Similarity network construction |
| 3a. Quantify | `tessera_communities()` | Community detection (Leiden/Louvain/spectral/NJ) |
| 3b. Quantify | `tessera_divergence()` | Regime divergence scoring |

Steps 4-6 (Predict, Test, Generalize) involving GNN, LSTM, and ensemble
prediction are available in the Python implementation and planned for a
future R release using R torch.

## Validated Domains

TESSERA has been validated across three domains:

- **Ecology**: Global Coral Bleaching Database (van Woesik et al. 2022)
- **Manufacturing**: UCI SECOM semiconductor dataset
- **HPC**: Synthetic + NOMAD monitoring data

## Citation

```
Tonini, J. & Parish, C. (2026). TESSERA: Topology-Enhanced Similarity and
Structure for Entity Regime Analysis. Nature Computational Science. (submitted)
```

## License

MIT
