#' Construct Similarity Networks
#'
#' Step 2b of the TESSERA pipeline. Converts similarity matrices into graph
#' objects using kNN, threshold, or weighted methods. The resulting network
#' encodes the topological structure that community detection will partition.
#'
#' @param x A `tessera_similarity` object.
#' @param method Character. Network construction method:
#'   `"knn"` (default) — connect each entity to its k nearest neighbors;
#'   `"threshold"` — connect entities with similarity above a threshold;
#'   `"weighted"` — fully connected weighted graph (only for small datasets).
#' @param k Integer. Number of neighbors for kNN method. Default 10. If the
#'   similarity object was already sparsified with a larger k, this further
#'   prunes edges.
#' @param threshold Numeric in \[0, 1\]. Similarity threshold for the threshold
#'   method. Default 0.5.
#' @param measure Character. Which similarity measure to use for network
#'   construction. Default uses the first available measure.
#'
#' @return A `tessera_network` object containing:
#' \describe{
#'   \item{graph}{An [igraph::igraph] graph object with similarity as edge weights}
#'   \item{method}{Construction method used}
#'   \item{measure}{Similarity measure used}
#'   \item{params}{List of parameters}
#' }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(200), nrow = 50, ncol = 4)
#' sim <- tessera_similarity(mat, measures = "cosine", k = 10)
#' net <- tessera_network(sim, method = "knn", k = 8)
#' summary(net)
#'
#' @export
tessera_network <- function(x,
                             method = c("knn", "threshold", "weighted"),
                             k = 10L,
                             threshold = 0.5,
                             measure = NULL) {
  method <- match.arg(method)

  if (!inherits(x, "tessera_similarity")) {
    stop("`x` must be a tessera_similarity object.", call. = FALSE)
  }

  # Handle per-bin similarity objects
  if (!is.null(x$by_bin)) {
    results <- lapply(x$by_bin, function(bin_sim) {
      tessera_network(bin_sim, method = method, k = k,
                      threshold = threshold, measure = measure)
    })
    names(results) <- names(x$by_bin)
    return(structure(
      list(
        by_bin = results,
        method = method,
        measure = results[[1]]$measure,
        n_bins = length(results),
        params = list(k = k, threshold = threshold)
      ),
      class = "tessera_network"
    ))
  }

  # Select measure
  if (is.null(measure)) {
    measure <- x$measures[1]
  }
  if (!measure %in% names(x$matrices)) {
    stop("Measure '", measure, "' not found. Available: ",
         paste(names(x$matrices), collapse = ", "), call. = FALSE)
  }

  sim_mat <- x$matrices[[measure]]
  if (inherits(sim_mat, "sparseMatrix")) {
    sim_mat <- as.matrix(sim_mat)
  }
  n <- nrow(sim_mat)

  # Build adjacency matrix based on method
  adj <- switch(method,
    knn = {
      k_use <- min(as.integer(k), n - 1L)
      a <- matrix(0, n, n)
      for (i in seq_len(n)) {
        row_vals <- sim_mat[i, ]
        row_vals[i] <- -Inf
        top_k <- order(row_vals, decreasing = TRUE)[seq_len(k_use)]
        a[i, top_k] <- sim_mat[i, top_k]
      }
      # Symmetrize
      pmax(a, t(a))
    },
    threshold = {
      a <- sim_mat
      a[a < threshold] <- 0
      diag(a) <- 0
      a
    },
    weighted = {
      a <- sim_mat
      diag(a) <- 0
      a[a < 0] <- 0
      a
    }
  )

  # Construct igraph object
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Store entity indices as vertex attribute
  igraph::V(g)$entity_idx <- seq_len(n)

  structure(
    list(
      graph = g,
      method = method,
      measure = measure,
      n_entities = n,
      params = list(k = k, threshold = threshold)
    ),
    class = "tessera_network"
  )
}


#' @export
print.tessera_network <- function(x, ...) {
  cat("TESSERA Network\n")
  if (!is.null(x$by_bin)) {
    cat(sprintf("  Bins: %d\n", x$n_bins))
    cat(sprintf("  Method: %s | Measure: %s\n", x$method, x$measure))
    for (bname in names(x$by_bin)) {
      bn <- x$by_bin[[bname]]
      cat(sprintf("  %s: %d vertices, %d edges (density %.3f)\n",
                  bname,
                  igraph::vcount(bn$graph),
                  igraph::ecount(bn$graph),
                  igraph::ecount(bn$graph) / (igraph::vcount(bn$graph) *
                    (igraph::vcount(bn$graph) - 1) / 2)))
    }
  } else {
    g <- x$graph
    cat(sprintf("  Vertices: %d | Edges: %d\n",
                igraph::vcount(g), igraph::ecount(g)))
    cat(sprintf("  Density: %.4f\n",
                igraph::ecount(g) / (igraph::vcount(g) *
                  (igraph::vcount(g) - 1) / 2)))
    cat(sprintf("  Method: %s | Measure: %s\n", x$method, x$measure))
    deg <- igraph::degree(g)
    cat(sprintf("  Degree: mean=%.1f, min=%d, max=%d\n",
                mean(deg), min(deg), max(deg)))
  }
  invisible(x)
}


#' @export
summary.tessera_network <- function(object, ...) {
  print(object)
  if (is.null(object$by_bin)) {
    g <- object$graph
    w <- igraph::E(g)$weight
    cat(sprintf("\n  Edge weights: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
                mean(w), sd(w), min(w), max(w)))
    cat(sprintf("  Components: %d\n", igraph::components(g)$no))
  }
  invisible(object)
}


#' @export
plot.tessera_network <- function(x, layout = "auto", ...) {
  if (!is.null(x$by_bin)) {
    message("Plotting first bin. Use $by_bin[[n]] to plot specific bins.")
    x <- x$by_bin[[1]]
  }

  g <- x$graph

  if (layout == "auto") {
    if (igraph::vcount(g) < 200) {
      l <- igraph::layout_with_fr(g)
    } else {
      l <- igraph::layout_with_drl(g)
    }
  }

  # Scale edge widths by weight
  w <- igraph::E(g)$weight
  w_scaled <- 0.5 + 2.5 * (w - min(w)) / max(max(w) - min(w), 1e-10)

  plot(g, layout = l,
       vertex.size = 3,
       vertex.label = NA,
       edge.width = w_scaled,
       edge.color = grDevices::adjustcolor("#4A90D9", alpha.f = 0.3),
       vertex.color = "#4A90D9",
       main = sprintf("Similarity Network (%s, %s)", x$method, x$measure),
       ...)
  invisible(x)
}
