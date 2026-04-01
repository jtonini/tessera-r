#' tessera: Topology-Enhanced Similarity and Structure for Entity Regime Analysis
#'
#' Implements the TESSERA framework for detecting regime boundaries in complex
#' systems using similarity network topology. The package provides a pipeline
#' of functions that compute multi-measure entity similarity, construct and
#' analyze networks, identify community structure, and quantify regime
#' divergence across temporal bins.
#'
#' @section Pipeline functions:
#' The core workflow follows TESSERA's 6-step methodology (Phase 1 covers
#' Steps 1-3):
#'
#' \describe{
#'   \item{[tessera_characterize()]}{Step 1: Characterize entities via temporal
#'     binning and feature preparation}
#'   \item{[tessera_similarity()]}{Step 2a: Compute pairwise entity similarity
#'     using multiple measures}
#'   \item{[tessera_network()]}{Step 2b: Construct similarity networks via kNN,
#'     threshold, or weighted methods}
#'   \item{[tessera_communities()]}{Step 3a: Detect community structure using
#'     Leiden, spectral, or Neighbor-Joining}
#'   \item{[tessera_divergence()]}{Step 3b: Quantify regime divergence across
#'     temporal bins}
#' }
#'
#' @section Convenience wrapper:
#' [tessera()] runs the full pipeline with sensible defaults or named presets.
#'
#' @section Key design choices:
#' \describe{
#'   \item{Sparse by default}{Similarity matrices use sparse kNN representation
#'     via [Matrix::sparseMatrix()], enabling analysis of large datasets
#'     (10k+ entities) on standard hardware.}
#'   \item{Quantile-based binning}{Temporal bins use quantile boundaries by
#'     default, ensuring balanced entity counts across bins even with uneven
#'     temporal distributions.}
#'   \item{Inspectable objects}{Every pipeline step returns an S3 object with
#'     `print()`, `summary()`, and `plot()` methods.}
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix
#'   cluster_louvain components membership modularity layout_with_fr
#'   layout_with_drl V E vcount ecount degree is_igraph
#' @importFrom Matrix sparseMatrix t nnzero
#' @importFrom stats dist cor cov quantile sd complete.cases cutree as.dist
#' @importFrom graphics barplot image par plot
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom utils object.size
#' @importFrom ape nj
#' @importFrom vegan vegdist
#' @docType package
#' @name tessera-package
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # Check for optional dependencies and note availability
  has_leiden <- requireNamespace("leiden", quietly = TRUE)
  has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
  has_kernlab <- requireNamespace("kernlab", quietly = TRUE)

  assign(".tessera_env", new.env(parent = emptyenv()), envir = parent.env(environment()))
  .tessera_env$has_leiden <- has_leiden
  .tessera_env$has_ggplot2 <- has_ggplot2
  .tessera_env$has_kernlab <- has_kernlab
}

# Internal environment for package state
.tessera_env <- NULL

# Suppress R CMD check NOTEs for ggplot2 .data pronoun
utils::globalVariables(".data")
