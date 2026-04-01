#' Compute Multi-Measure Entity Similarity
#'
#' Step 2a of the TESSERA pipeline. Computes pairwise entity similarity using
#' one or more of seven measures, returning a sparse kNN representation by
#' default. This is the computational core of TESSERA — the similarity
#' structure encodes topological information that community detection and
#' regime divergence subsequently exploit.
#'
#' @param x Either a `tessera_characterize` object (from [tessera_characterize()])
#'   or a numeric data.frame/matrix of entity features.
#' @param measures Character vector specifying which similarity measures to
#'   compute. Options: `"simpson"`, `"cosine"`, `"bray_curtis"`, `"jaccard"`,
#'   `"pearson"`, `"euclidean"`, `"mahalanobis"`, or `"all"`. Default is
#'   `"simpson"`. Use `"adaptive"` for the learned convex combination (Phase 2).
#' @param k Integer. Number of nearest neighbors for sparse representation.
#'   Default is 15. Only the top-k most similar entities per row are retained.
#' @param bin Integer or `NULL`. If `x` is a `tessera_characterize` object,
#'   compute similarity within this specific bin. If `NULL`, compute for each
#'   bin and return a list.
#' @param sparse Logical. Whether to return sparse kNN matrices. Default `TRUE`.
#'   Set to `FALSE` for small datasets where dense matrices are acceptable.
#' @param n_quantile_bins Integer. Number of quantile bins for Simpson similarity.
#'   Default is 4. Features are discretized into this many quantile-based bins
#'   before computing Simpson overlap.
#'
#' @return A `tessera_similarity` object containing:
#' \describe{
#'   \item{matrices}{Named list of similarity matrices (sparse or dense), one
#'     per measure. If `bin = NULL` and input is a characterize object, this is
#'     a nested list: `matrices[[bin]][[measure]]`.}
#'   \item{measures}{Character vector of measures computed}
#'   \item{k}{The k parameter used for sparse representation}
#'   \item{n_entities}{Number of entities}
#'   \item{params}{List of parameters}
#' }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(200), nrow = 50, ncol = 4)
#' sim <- tessera_similarity(mat, measures = c("cosine", "euclidean"), k = 10)
#' summary(sim)
#'
#' @export
tessera_similarity <- function(x,
                                measures = "simpson",
                                k = 15L,
                                bin = NULL,
                                sparse = TRUE,
                                n_quantile_bins = 4L) {

  # Extract feature matrix from characterize object or raw data
  if (inherits(x, "tessera_characterize")) {
    if (is.null(bin)) {
      # Compute per-bin
      bins <- sort(unique(x$data$.bin))
      results <- lapply(bins, function(b) {
        tessera_similarity(x, measures = measures, k = k, bin = b,
                           sparse = sparse, n_quantile_bins = n_quantile_bins)
      })
      names(results) <- paste0("bin_", bins)
      return(structure(
        list(
          by_bin = results,
          measures = results[[1]]$measures,
          k = k,
          n_bins = length(bins),
          params = list(sparse = sparse, n_quantile_bins = n_quantile_bins)
        ),
        class = "tessera_similarity"
      ))
    }
    idx <- x$data$.bin == bin
    feat_mat <- as.matrix(x$data[idx, x$features, drop = FALSE])
  } else {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
    feat_mat <- x
  }

  n <- nrow(feat_mat)
  k <- min(as.integer(k), n - 1L)

  if ("all" %in% measures) {
    measures <- c("simpson", "cosine", "bray_curtis", "jaccard",
                  "pearson", "euclidean", "mahalanobis")
  }
  valid_measures <- c("simpson", "cosine", "bray_curtis", "jaccard",
                      "pearson", "euclidean", "mahalanobis")
  bad <- setdiff(measures, c(valid_measures, "adaptive"))
  if (length(bad) > 0L) {
    stop("Unknown measures: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # Compute each similarity measure
  sim_list <- list()

  for (m in measures) {
    sim_mat <- switch(m,
      simpson   = .simpson_similarity(feat_mat, n_quantile_bins),
      cosine    = .cosine_similarity(feat_mat),
      bray_curtis = .bray_curtis_similarity(feat_mat),
      jaccard   = .jaccard_similarity(feat_mat, n_quantile_bins),
      pearson   = .pearson_similarity(feat_mat),
      euclidean = .euclidean_similarity(feat_mat),
      mahalanobis = .mahalanobis_similarity(feat_mat),
      stop("Not implemented: ", m)
    )

    # Sparsify via kNN
    if (sparse && n > k + 1L) {
      sim_mat <- .sparsify_knn(sim_mat, k)
    }

    sim_list[[m]] <- sim_mat
  }

  structure(
    list(
      matrices = sim_list,
      measures = measures,
      k = k,
      n_entities = n,
      bin = bin,
      params = list(sparse = sparse, n_quantile_bins = n_quantile_bins)
    ),
    class = "tessera_similarity"
  )
}


# --- Internal similarity functions ---

#' Simpson similarity: proportion of shared quantile-bin memberships
#' @noRd
.simpson_similarity <- function(feat_mat, n_bins = 4L) {
  n <- nrow(feat_mat)
  d <- ncol(feat_mat)

  # Discretize each feature into quantile bins
  binned <- apply(feat_mat, 2, function(col) {
    breaks <- unique(stats::quantile(col, probs = seq(0, 1, length.out = n_bins + 1L)))
    as.integer(cut(col, breaks = breaks, include.lowest = TRUE, labels = FALSE))
  })

  # Simpson overlap: fraction of features where two entities share a bin
  # Vectorized via matrix operations on one-hot encoded bins
  sim <- matrix(0, n, n)
  for (j in seq_len(d)) {
    for (b in seq_len(n_bins)) {
      indicator <- as.numeric(binned[, j] == b)
      sim <- sim + tcrossprod(indicator)
    }
  }
  sim <- sim / d
  diag(sim) <- 1
  sim
}

#' Cosine similarity
#' @noRd
.cosine_similarity <- function(feat_mat) {
  norms <- sqrt(rowSums(feat_mat^2))
  norms[norms == 0] <- 1  # avoid division by zero
  normed <- feat_mat / norms
  sim <- tcrossprod(normed)
  # Clamp to [0, 1] for consistency
  sim <- pmax(sim, 0)
  diag(sim) <- 1
  sim
}

#' Bray-Curtis similarity (1 - dissimilarity)
#' Uses vegan if available, otherwise manual computation
#' @noRd
.bray_curtis_similarity <- function(feat_mat) {
  # Shift to non-negative (Bray-Curtis requires non-negative values)
  shifted <- feat_mat
  col_mins <- apply(feat_mat, 2, min, na.rm = TRUE)
  for (j in seq_len(ncol(shifted))) {
    if (col_mins[j] < 0) {
      shifted[, j] <- shifted[, j] - col_mins[j]
    }
  }

  if (requireNamespace("vegan", quietly = TRUE)) {
    d <- as.matrix(vegan::vegdist(shifted, method = "bray"))
  } else {
    # Manual Bray-Curtis
    n <- nrow(shifted)
    d <- matrix(0, n, n)
    row_sums_mat <- rowSums(shifted)
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        num <- sum(abs(shifted[i, ] - shifted[j, ]))
        denom <- row_sums_mat[i] + row_sums_mat[j]
        d[i, j] <- d[j, i] <- if (denom > 0) num / denom else 0
      }
    }
  }
  1 - d
}

#' Jaccard similarity on quantile-binned features
#' @noRd
.jaccard_similarity <- function(feat_mat, n_bins = 4L) {
  n <- nrow(feat_mat)

  # Discretize and create binary presence/absence per bin
  binned <- apply(feat_mat, 2, function(col) {
    breaks <- unique(stats::quantile(col, probs = seq(0, 1, length.out = n_bins + 1L)))
    as.integer(cut(col, breaks = breaks, include.lowest = TRUE, labels = FALSE))
  })

  # Create binary feature matrix (one-hot across all feature-bin combos)
  bin_ids <- paste0(
    rep(seq_len(ncol(binned)), each = nrow(binned)),
    "_",
    as.vector(binned)
  )
  unique_ids <- unique(bin_ids)
  binary_mat <- matrix(0L, n, length(unique_ids))
  for (i in seq_len(n)) {
    entity_ids <- paste0(seq_len(ncol(binned)), "_", binned[i, ])
    binary_mat[i, match(entity_ids, unique_ids)] <- 1L
  }

  # Jaccard = intersection / union
  intersection <- tcrossprod(binary_mat)
  row_sums <- rowSums(binary_mat)
  union_mat <- outer(row_sums, row_sums, "+") - intersection
  sim <- intersection / pmax(union_mat, 1)
  diag(sim) <- 1
  sim
}

#' Pearson correlation similarity: (1 + r) / 2 mapped to [0, 1]
#' @noRd
.pearson_similarity <- function(feat_mat) {
  r <- stats::cor(t(feat_mat))
  r[is.na(r)] <- 0
  (1 + r) / 2
}

#' Euclidean similarity: 1 / (1 + d)
#' @noRd
.euclidean_similarity <- function(feat_mat) {
  d <- as.matrix(stats::dist(feat_mat, method = "euclidean"))
  1 / (1 + d)
}

#' Mahalanobis-based similarity
#' Uses pooled covariance. Falls back to Euclidean if covariance is singular.
#' @noRd
.mahalanobis_similarity <- function(feat_mat) {
  cov_mat <- stats::cov(feat_mat)
  cov_inv <- tryCatch(
    solve(cov_mat),
    error = function(e) {
      # Regularize with ridge
      solve(cov_mat + diag(1e-6, ncol(feat_mat)))
    }
  )

  n <- nrow(feat_mat)
  d <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    diff_mat <- sweep(feat_mat[(i + 1L):n, , drop = FALSE], 2, feat_mat[i, ])
    mah_d <- rowSums((diff_mat %*% cov_inv) * diff_mat)
    mah_d <- pmax(mah_d, 0)  # numerical guard
    d[i, (i + 1L):n] <- sqrt(mah_d)
    d[(i + 1L):n, i] <- d[i, (i + 1L):n]
  }

  1 / (1 + d)
}


#' Sparsify a similarity matrix by retaining only top-k neighbors per row
#' @noRd
.sparsify_knn <- function(sim_mat, k) {
  n <- nrow(sim_mat)
  i_idx <- integer(0)
  j_idx <- integer(0)
  x_vals <- numeric(0)

  for (row in seq_len(n)) {
    row_vals <- sim_mat[row, ]
    row_vals[row] <- -Inf  # exclude self
    top_k <- order(row_vals, decreasing = TRUE)[seq_len(k)]
    i_idx <- c(i_idx, rep(row, k))
    j_idx <- c(j_idx, top_k)
    x_vals <- c(x_vals, row_vals[top_k])
  }

  # Symmetrize: include edge if EITHER direction has it in top-k
  sparse_mat <- Matrix::sparseMatrix(
    i = i_idx, j = j_idx, x = x_vals,
    dims = c(n, n), symmetric = FALSE
  )
  # Symmetrize by taking pairwise max
  sparse_sym <- pmax(sparse_mat, Matrix::t(sparse_mat))
  sparse_sym
}


# --- S3 methods ---

#' @export
print.tessera_similarity <- function(x, ...) {
  cat("TESSERA Similarity\n")
  if (!is.null(x$by_bin)) {
    cat(sprintf("  Bins: %d\n", x$n_bins))
    cat(sprintf("  Measures: %s\n", paste(x$measures, collapse = ", ")))
    cat(sprintf("  k (sparse): %d\n", x$k))
  } else {
    cat(sprintf("  Entities: %d\n", x$n_entities))
    cat(sprintf("  Measures: %s\n", paste(x$measures, collapse = ", ")))
    cat(sprintf("  k (sparse): %d\n", x$k))
    cat(sprintf("  Bin: %s\n", if (is.null(x$bin)) "all" else x$bin))
    for (m in names(x$matrices)) {
      mat <- x$matrices[[m]]
      if (inherits(mat, "sparseMatrix")) {
        nnz <- Matrix::nnzero(mat)
        density <- nnz / (x$n_entities^2)
        cat(sprintf("  %s: sparse (%.1f%% density, %.1f KB)\n",
                    m, 100 * density,
                    object.size(mat) / 1024))
      } else {
        cat(sprintf("  %s: dense (%.1f MB)\n",
                    m, object.size(mat) / 1024^2))
      }
    }
  }
  invisible(x)
}


#' @export
summary.tessera_similarity <- function(object, ...) {
  cat("TESSERA Similarity Summary\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")

  if (!is.null(object$by_bin)) {
    for (bname in names(object$by_bin)) {
      cat(sprintf("\n--- %s ---\n", bname))
      summary(object$by_bin[[bname]])
    }
  } else {
    cat(sprintf("Entities: %d | Bin: %s\n",
                object$n_entities,
                if (is.null(object$bin)) "all" else object$bin))
    for (m in names(object$matrices)) {
      mat <- object$matrices[[m]]
      vals <- if (inherits(mat, "sparseMatrix")) {
        mat@x[mat@x > 0]
      } else {
        mat[upper.tri(mat)]
      }
      cat(sprintf("\n  %s:\n", m))
      cat(sprintf("    Range: [%.4f, %.4f]\n", min(vals), max(vals)))
      cat(sprintf("    Mean:  %.4f\n", mean(vals)))
      cat(sprintf("    SD:    %.4f\n", sd(vals)))
    }
  }
  invisible(object)
}


#' @export
plot.tessera_similarity <- function(x, measure = NULL, ...) {
  if (!is.null(x$by_bin)) {
    message("Plotting first bin. Use $by_bin[[n]] to plot specific bins.")
    x <- x$by_bin[[1]]
  }

  if (is.null(measure)) measure <- x$measures[1]
  mat <- x$matrices[[measure]]
  if (inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
  }

  image(mat, main = paste("Similarity:", measure),
        xlab = "Entity", ylab = "Entity",
        col = grDevices::hcl.colors(50, "viridis"))
  invisible(x)
}
