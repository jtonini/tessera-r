#' Run the Full TESSERA Pipeline
#'
#' Convenience wrapper that runs Steps 1–3 of the TESSERA framework in a single
#' call: characterize entities, compute similarity, build networks, detect
#' communities, and quantify regime divergence. For fine-grained control, use
#' the individual pipeline functions.
#'
#' @param data A data.frame of entity features with a temporal column.
#' @param time_col Character. Name of the temporal column.
#' @param entity_col Character or `NULL`. Name of the entity ID column.
#' @param preset Character or `NULL`. Named parameter preset:
#'   `"ecology"` — optimized for ecological observational data (quantile bins,
#'     Bray-Curtis + Simpson measures, Leiden communities);
#'   `"manufacturing"` — optimized for sensor/quality data (fixed bins,
#'     Mahalanobis + Euclidean, spectral communities);
#'   `"hpc"` — optimized for HPC job data (quantile bins, cosine + Pearson,
#'     Louvain communities);
#'   `NULL` — uses defaults.
#' @param measures Character vector. Similarity measures (see [tessera_similarity()]).
#' @param n_bins Integer. Number of temporal bins.
#' @param bin_method Character. `"quantile"` or `"fixed"`.
#' @param network_method Character. Network construction method.
#' @param k Integer. kNN parameter.
#' @param community_algorithm Character. Community detection algorithm.
#' @param resolution Numeric. Resolution parameter for community detection.
#' @param divergence_method Character. Divergence measure.
#' @param divergence_reference Character. Divergence reference strategy.
#' @param scale Logical. Standardize features. Default `TRUE`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return A `tessera_result` object containing all intermediate results:
#'   `$characterize`, `$similarity`, `$network`, `$communities`, `$divergence`.
#'
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   year = sample(2000:2020, 500, replace = TRUE),
#'   temp = rnorm(500, 25, 3),
#'   depth = runif(500, 0, 50),
#'   salinity = rnorm(500, 35, 1)
#' )
#' result <- tessera(df, time_col = "year")
#' result
#' plot(result)
#'
#' @export
tessera <- function(data,
                     time_col = NULL,
                     entity_col = NULL,
                     preset = NULL,
                     measures = NULL,
                     n_bins = NULL,
                     bin_method = NULL,
                     network_method = NULL,
                     k = NULL,
                     community_algorithm = NULL,
                     resolution = NULL,
                     divergence_method = NULL,
                     divergence_reference = NULL,
                     scale = TRUE,
                     verbose = TRUE) {

  # Apply preset defaults (user-specified params override)
  defaults <- .get_preset_defaults(preset)
  if (is.null(measures)) measures <- defaults$measures
  if (is.null(n_bins)) n_bins <- defaults$n_bins
  if (is.null(bin_method)) bin_method <- defaults$bin_method
  if (is.null(network_method)) network_method <- defaults$network_method
  if (is.null(k)) k <- defaults$k
  if (is.null(community_algorithm)) community_algorithm <- defaults$community_algorithm
  if (is.null(resolution)) resolution <- defaults$resolution
  if (is.null(divergence_method)) divergence_method <- defaults$divergence_method
  if (is.null(divergence_reference)) divergence_reference <- defaults$divergence_reference

  t0 <- proc.time()

  # Step 1: Characterize
  if (verbose) message("Step 1/5: Characterizing entities...")
  char <- tessera_characterize(
    data, time_col = time_col, entity_col = entity_col,
    n_bins = n_bins, bin_method = bin_method, scale = scale
  )
  if (verbose) message(sprintf("  %d entities, %d bins, %d features",
                                nrow(char$data), nrow(char$bins),
                                length(char$features)))

  # Step 2a: Similarity
  if (verbose) message("Step 2/5: Computing similarity (", paste(measures, collapse = ", "), ")...")
  sim <- tessera_similarity(
    char, measures = measures, k = k
  )

  # Step 3: Network
  if (verbose) message("Step 3/5: Building networks (", network_method, ")...")
  net <- tessera_network(
    sim, method = network_method, k = k
  )

  # Step 4: Communities
  if (verbose) message("Step 4/5: Detecting communities (", community_algorithm, ")...")
  com <- tessera_communities(
    net, algorithm = community_algorithm, resolution = resolution
  )

  # Step 5: Divergence
  if (verbose) message("Step 5/5: Computing regime divergence (", divergence_method, ")...")
  div <- tessera_divergence(
    com, method = divergence_method, reference = divergence_reference
  )

  elapsed <- (proc.time() - t0)["elapsed"]
  if (verbose) message(sprintf("Done in %.1f seconds.", elapsed))

  structure(
    list(
      characterize = char,
      similarity = sim,
      network = net,
      communities = com,
      divergence = div,
      preset = preset,
      elapsed = elapsed
    ),
    class = "tessera_result"
  )
}


#' Get preset parameter defaults
#' @noRd
.get_preset_defaults <- function(preset = NULL) {
  base <- list(
    measures = "simpson",
    n_bins = 4L,
    bin_method = "quantile",
    network_method = "knn",
    k = 15L,
    community_algorithm = "louvain",  # safe default (no extra deps)
    resolution = 1.0,
    divergence_method = "nmi",
    divergence_reference = "adjacent"
  )

  if (is.null(preset)) return(base)

  switch(preset,
    ecology = {
      base$measures <- c("bray_curtis", "simpson")
      base$bin_method <- "quantile"
      base$community_algorithm <- "leiden"
      base$k <- 15L
      base
    },
    manufacturing = {
      base$measures <- c("mahalanobis", "euclidean")
      base$bin_method <- "fixed"
      base$community_algorithm <- "spectral"
      base$k <- 10L
      base
    },
    hpc = {
      base$measures <- c("cosine", "pearson")
      base$bin_method <- "quantile"
      base$community_algorithm <- "louvain"
      base$k <- 12L
      base
    },
    stop("Unknown preset: '", preset, "'. Options: ecology, manufacturing, hpc",
         call. = FALSE)
  )
}


# --- S3 methods for tessera_result ---

#' @export
print.tessera_result <- function(x, ...) {
  cat("TESSERA Analysis Result\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  if (!is.null(x$preset)) {
    cat(sprintf("Preset: %s\n", x$preset))
  }
  cat(sprintf("Entities: %d | Bins: %d | Features: %d\n",
              nrow(x$characterize$data),
              nrow(x$characterize$bins),
              length(x$characterize$features)))
  cat(sprintf("Measures: %s\n", paste(x$similarity$measures, collapse = ", ")))
  cat(sprintf("Network: %s | Communities: %s\n",
              x$network$method, x$communities$algorithm))

  # Divergence summary
  div <- x$divergence
  if (is.matrix(div$scores)) {
    cat(sprintf("Divergence: %s (all-pairs, mean=%.4f)\n",
                div$method, mean(div$scores[upper.tri(div$scores)])))
  } else {
    cat(sprintf("Divergence: %s (mean=%.4f, max=%.4f)\n",
                div$method, mean(div$scores), max(div$scores)))
  }
  cat(sprintf("Elapsed: %.1f seconds\n", x$elapsed))
  invisible(x)
}


#' @export
summary.tessera_result <- function(object, ...) {
  print(object)
  cat("\n--- Bin Details ---\n")
  print(object$characterize$bins, row.names = FALSE)
  cat("\n--- Divergence Scores ---\n")
  print(object$divergence)
  invisible(object)
}


#' @export
plot.tessera_result <- function(x, which = "divergence", ...) {
  switch(which,
    divergence = plot(x$divergence, ...),
    network    = plot(x$network, ...),
    communities = plot(x$communities, ...),
    similarity = plot(x$similarity, ...),
    bins       = plot(x$characterize, ...),
    all = {
      oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      on.exit(par(oldpar))
      plot(x$characterize, ...)
      plot(x$similarity, ...)
      plot(x$communities, ...)
      plot(x$divergence, ...)
    },
    stop("Unknown plot type: '", which, "'. ",
         "Options: divergence, network, communities, similarity, bins, all",
         call. = FALSE)
  )
  invisible(x)
}
