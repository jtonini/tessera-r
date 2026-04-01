#' Detect Community Structure
#'
#' Step 3a of the TESSERA pipeline. Identifies communities (clusters) in the
#' similarity network using Leiden, spectral clustering, or Neighbor-Joining.
#' Community assignments encode the regime structure that divergence analysis
#' subsequently quantifies.
#'
#' @param x A `tessera_network` object.
#' @param algorithm Character. Community detection algorithm:
#'   `"leiden"` (default, requires the `leiden` package),
#'   `"louvain"` (igraph built-in),
#'   `"spectral"` (requires the `kernlab` package),
#'   `"nj"` (Neighbor-Joining via `ape`, produces a hierarchical tree).
#' @param resolution Numeric. Resolution parameter for Leiden/Louvain. Higher
#'   values produce more communities. Default 1.0.
#' @param n_communities Integer or `NULL`. Number of communities for spectral
#'   clustering. If `NULL`, estimated via eigengap heuristic. Default `NULL`.
#' @param max_communities Integer. Upper bound on communities for one-hot
#'   encoding in downstream GNN features (Phase 2). Default 20. Communities
#'   beyond this threshold are merged into an "other" category.
#'
#' @return A `tessera_communities` object containing:
#' \describe{
#'   \item{membership}{Integer vector of community assignments per entity}
#'   \item{n_communities}{Number of communities detected}
#'   \item{modularity}{Modularity score (for Leiden/Louvain)}
#'   \item{sizes}{Named integer vector of community sizes}
#'   \item{graph}{The igraph object with community assignments as vertex attribute}
#'   \item{tree}{For NJ: the phylo tree object from ape}
#'   \item{algorithm}{Algorithm used}
#' }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(200), nrow = 50, ncol = 4)
#' sim <- tessera_similarity(mat, measures = "cosine", k = 10)
#' net <- tessera_network(sim, method = "knn", k = 8)
#' com <- tessera_communities(net, algorithm = "louvain")
#' summary(com)
#'
#' @export
tessera_communities <- function(x,
                                 algorithm = c("leiden", "louvain", "spectral", "nj"),
                                 resolution = 1.0,
                                 n_communities = NULL,
                                 max_communities = 20L) {
  algorithm <- match.arg(algorithm)

  if (!inherits(x, "tessera_network")) {
    stop("`x` must be a tessera_network object.", call. = FALSE)
  }

  # Handle per-bin network objects
  if (!is.null(x$by_bin)) {
    results <- lapply(x$by_bin, function(bin_net) {
      tessera_communities(bin_net, algorithm = algorithm,
                          resolution = resolution,
                          n_communities = n_communities,
                          max_communities = max_communities)
    })
    names(results) <- names(x$by_bin)
    return(structure(
      list(
        by_bin = results,
        algorithm = algorithm,
        n_bins = length(results),
        params = list(resolution = resolution,
                      max_communities = max_communities)
      ),
      class = "tessera_communities"
    ))
  }

  g <- x$graph
  n <- igraph::vcount(g)
  tree <- NULL

  result <- switch(algorithm,
    leiden = {
      if (!requireNamespace("leiden", quietly = TRUE)) {
        stop("Package 'leiden' required for Leiden algorithm. ",
             "Install with: install.packages('leiden')", call. = FALSE)
      }
      adj <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
      tryCatch(
        {
          mem <- leiden::leiden(adj, resolution_parameter = resolution)
          list(mem = mem, mod = igraph::modularity(g, mem), alg = "leiden")
        },
        error = function(e) {
          warning("Leiden failed (", conditionMessage(e),
                  "). Falling back to Louvain.", call. = FALSE)
          cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight,
                                         resolution = resolution)
          list(mem = igraph::membership(cl),
               mod = igraph::modularity(cl),
               alg = "louvain")
        }
      )
    },

    louvain = {
      cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight,
                                     resolution = resolution)
      list(mem = igraph::membership(cl),
           mod = igraph::modularity(cl),
           alg = "louvain")
    },

    spectral = {
      if (!requireNamespace("kernlab", quietly = TRUE)) {
        stop("Package 'kernlab' required for spectral clustering. ",
             "Install with: install.packages('kernlab')", call. = FALSE)
      }
      adj <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)

      if (is.null(n_communities)) {
        deg_mat <- diag(rowSums(adj))
        laplacian <- deg_mat - adj
        deg_inv_sqrt <- diag(1 / sqrt(pmax(diag(deg_mat), 1e-10)))
        norm_lap <- deg_inv_sqrt %*% laplacian %*% deg_inv_sqrt
        eigs <- eigen(norm_lap, symmetric = TRUE, only.values = TRUE)$values
        eigs <- sort(eigs)
        max_check <- min(20L, length(eigs) - 1L)
        gaps <- diff(eigs[seq_len(max_check + 1L)])
        n_communities <- which.max(gaps)
        n_communities <- max(2L, min(n_communities, max_communities))
      }

      sc <- kernlab::specc(as.matrix(adj), centers = as.integer(n_communities))
      list(mem = as.integer(sc@.Data),
           mod = NA_real_,
           alg = "spectral")
    },

    nj = {
      adj <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
      adj[adj == 0] <- NA
      mean_sim <- mean(adj, na.rm = TRUE)
      adj[is.na(adj)] <- mean_sim
      dist_mat <- 1 - adj
      diag(dist_mat) <- 0
      dist_mat <- pmax(dist_mat, 0)

      nj_tree <- ape::nj(as.dist(dist_mat))

      if (is.null(n_communities)) n_communities <- min(10L, n %/% 5L)
      n_communities <- max(2L, n_communities)
      mem <- cutree(stats::hclust(as.dist(dist_mat), method = "average"),
                    k = n_communities)
      list(mem = mem, mod = NA_real_, alg = "nj", tree = nj_tree)
    }
  )

  membership <- result$mem
  mod <- result$mod
  algorithm <- result$alg
  if (!is.null(result$tree)) tree <- result$tree

  # Cap community count
  if (max(membership) > max_communities) {
    sizes <- sort(table(membership), decreasing = TRUE)
    keep <- as.integer(names(sizes)[seq_len(max_communities)])
    membership[!membership %in% keep] <- max_communities + 1L
    # Relabel sequentially
    old <- sort(unique(membership))
    membership <- match(membership, old)
    message(sprintf("Capped communities from %d to %d (merged small communities).",
                    length(sizes), max_communities))
  }

  # Store membership as vertex attribute
  igraph::V(g)$community <- membership

  sizes <- table(membership)
  names(sizes) <- paste0("C", names(sizes))

  structure(
    list(
      membership = membership,
      n_communities = max(membership),
      modularity = mod,
      sizes = sizes,
      graph = g,
      tree = tree,
      algorithm = algorithm,
      n_entities = n,
      params = list(resolution = resolution,
                    max_communities = max_communities)
    ),
    class = "tessera_communities"
  )
}


#' @export
print.tessera_communities <- function(x, ...) {
  cat("TESSERA Communities\n")
  if (!is.null(x$by_bin)) {
    cat(sprintf("  Bins: %d | Algorithm: %s\n", x$n_bins, x$algorithm))
    for (bname in names(x$by_bin)) {
      bc <- x$by_bin[[bname]]
      cat(sprintf("  %s: %d communities", bname, bc$n_communities))
      if (!is.na(bc$modularity)) {
        cat(sprintf(" (Q=%.3f)", bc$modularity))
      }
      cat("\n")
    }
  } else {
    cat(sprintf("  Algorithm: %s\n", x$algorithm))
    cat(sprintf("  Communities: %d | Entities: %d\n",
                x$n_communities, x$n_entities))
    if (!is.na(x$modularity)) {
      cat(sprintf("  Modularity: %.4f\n", x$modularity))
    }
    cat(sprintf("  Sizes: %s\n",
                paste(sprintf("%s=%d", names(x$sizes), x$sizes),
                      collapse = ", ")))
  }
  invisible(x)
}


#' @export
summary.tessera_communities <- function(object, ...) {
  print(object)
  if (is.null(object$by_bin)) {
    cat("\nCommunity size distribution:\n")
    s <- as.integer(object$sizes)
    cat(sprintf("  Mean: %.1f | SD: %.1f | Min: %d | Max: %d\n",
                mean(s), sd(s), min(s), max(s)))
    if (!is.null(object$tree)) {
      cat("  NJ tree available ($tree)\n")
    }
  }
  invisible(object)
}


#' @export
plot.tessera_communities <- function(x, ...) {
  if (!is.null(x$by_bin)) {
    message("Plotting first bin. Use $by_bin[[n]] to plot specific bins.")
    x <- x$by_bin[[1]]
  }

  g <- x$graph
  n_com <- x$n_communities
  palette <- if (n_com <= 8) {
    # Okabe-Ito colorblind-safe
    c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#999999")
  } else {
    grDevices::hcl.colors(n_com, "Set2")
  }

  if (igraph::vcount(g) < 200) {
    l <- igraph::layout_with_fr(g)
  } else {
    l <- igraph::layout_with_drl(g)
  }

  plot(g, layout = l,
       vertex.size = 4,
       vertex.label = NA,
       vertex.color = palette[x$membership],
       edge.width = 0.5,
       edge.color = grDevices::adjustcolor("grey60", alpha.f = 0.2),
       main = sprintf("Communities (%s, n=%d)", x$algorithm, n_com),
       ...)
  invisible(x)
}
