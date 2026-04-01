#' Quantify Regime Divergence
#'
#' Step 3b of the TESSERA pipeline. Computes regime divergence (RD) scores by
#' comparing community structure across temporal bins. RD captures how the
#' topology of entity relationships shifts over time — the core analytical
#' signal that TESSERA exploits.
#'
#' @param x A `tessera_communities` object with per-bin results (i.e., created
#'   from a pipeline that processed multiple temporal bins).
#' @param method Character. Divergence measure:
#'   `"nmi"` (default) — Normalized Mutual Information between bin community
#'     assignments;
#'   `"vi"` — Variation of Information;
#'   `"composition"` — community composition shift (proportion of entities
#'     changing community membership between adjacent bins).
#' @param reference Character. Reference strategy for divergence:
#'   `"adjacent"` (default) — compare each bin to the previous bin;
#'   `"first"` — compare each bin to the first bin;
#'   `"all_pairs"` — compute all pairwise bin divergences.
#'
#' @return A `tessera_divergence` object containing:
#' \describe{
#'   \item{scores}{Numeric vector of divergence scores (one per bin transition
#'     for adjacent/first reference, or a matrix for all_pairs)}
#'   \item{entity_divergence}{Per-entity divergence contribution}
#'   \item{method}{Divergence method used}
#'   \item{reference}{Reference strategy used}
#' }
#'
#' @details
#' Regime divergence is the operational measurement of "ecotone" — a concept
#' borrowed from biogeography where it describes the transitional zone between
#' distinct ecological communities. In TESSERA's framework, regime boundaries
#' are detected where community topology shifts significantly between temporal
#' bins.
#'
#' The three failure geometries (boundary-concentrated, isolated cluster,
#' diffuse) manifest as distinct divergence patterns:
#' \describe{
#'   \item{Boundary-concentrated}{High RD in bins near the transition, low elsewhere}
#'   \item{Isolated cluster}{Sudden spike in RD at one bin, stable elsewhere}
#'   \item{Diffuse}{Gradually increasing RD across all bins}
#' }
#'
#' @examples
#' # Requires a per-bin communities object; see tessera() for full pipeline
#'
#' @export
tessera_divergence <- function(x,
                                method = c("nmi", "vi", "composition"),
                                reference = c("adjacent", "first", "all_pairs")) {
  method <- match.arg(method)
  reference <- match.arg(reference)

  if (!inherits(x, "tessera_communities")) {
    stop("`x` must be a tessera_communities object.", call. = FALSE)
  }
  if (is.null(x$by_bin)) {
    stop("Divergence requires per-bin community results. ",
         "Run the pipeline with temporal bins.", call. = FALSE)
  }

  bin_names <- names(x$by_bin)
  n_bins <- length(bin_names)

  if (n_bins < 2L) {
    stop("At least 2 bins are required for divergence analysis.", call. = FALSE)
  }

  # Extract membership vectors per bin
  memberships <- lapply(x$by_bin, function(bc) bc$membership)

  # Compute pairwise divergence between two membership vectors
  # Entities may differ between bins — align by shared entities or use

  # co-occurrence based methods
  compute_div <- function(mem_a, mem_b, method) {
    # For entities that exist in both bins, compare community assignments
    # If different sizes, use the shorter one (entities present in both)
    n <- min(length(mem_a), length(mem_b))
    ma <- mem_a[seq_len(n)]
    mb <- mem_b[seq_len(n)]

    switch(method,
      nmi = {
        # Normalized Mutual Information
        # 1 - NMI so higher = more divergent
        1 - .nmi(ma, mb)
      },
      vi = {
        # Variation of Information
        .variation_of_information(ma, mb)
      },
      composition = {
        # Fraction of entities that changed community
        sum(ma != mb) / n
      }
    )
  }

  # Compute divergence based on reference strategy
  if (reference == "adjacent") {
    scores <- numeric(n_bins - 1L)
    names(scores) <- paste0(bin_names[-n_bins], " -> ", bin_names[-1])
    for (i in seq_len(n_bins - 1L)) {
      scores[i] <- compute_div(memberships[[i]], memberships[[i + 1]], method)
    }
  } else if (reference == "first") {
    scores <- numeric(n_bins - 1L)
    names(scores) <- paste0(bin_names[1], " -> ", bin_names[-1])
    for (i in 2:n_bins) {
      scores[i - 1L] <- compute_div(memberships[[1]], memberships[[i]], method)
    }
  } else {
    # all_pairs
    scores <- matrix(0, n_bins, n_bins, dimnames = list(bin_names, bin_names))
    for (i in seq_len(n_bins - 1L)) {
      for (j in (i + 1L):n_bins) {
        d <- compute_div(memberships[[i]], memberships[[j]], method)
        scores[i, j] <- d
        scores[j, i] <- d
      }
    }
  }

  # Per-entity divergence (adjacent bins, contribution to regime shift)
  entity_div <- NULL
  if (reference == "adjacent") {
    # For each entity, count how many bin transitions changed its community
    max_n <- max(vapply(memberships, length, integer(1)))
    entity_div <- numeric(max_n)
    for (i in seq_len(n_bins - 1L)) {
      n_shared <- min(length(memberships[[i]]), length(memberships[[i + 1]]))
      changed <- memberships[[i]][seq_len(n_shared)] !=
                 memberships[[i + 1]][seq_len(n_shared)]
      entity_div[seq_len(n_shared)] <- entity_div[seq_len(n_shared)] +
                                        as.numeric(changed)
    }
    entity_div <- entity_div / (n_bins - 1L)  # normalize to [0, 1]
  }

  structure(
    list(
      scores = scores,
      entity_divergence = entity_div,
      method = method,
      reference = reference,
      n_bins = n_bins,
      bin_names = bin_names,
      params = list()
    ),
    class = "tessera_divergence"
  )
}


# --- Internal NMI computation ---

#' Normalized Mutual Information
#' @noRd
.nmi <- function(a, b) {
  n <- length(a)
  if (n == 0) return(0)

  # Contingency table
  tab <- table(a, b)
  p_ab <- tab / n
  p_a <- rowSums(tab) / n
  p_b <- colSums(tab) / n

  # Mutual Information
  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (p_ab[i, j] > 0) {
        mi <- mi + p_ab[i, j] * log(p_ab[i, j] / (p_a[i] * p_b[j]))
      }
    }
  }

  # Entropies
  h_a <- -sum(p_a[p_a > 0] * log(p_a[p_a > 0]))
  h_b <- -sum(p_b[p_b > 0] * log(p_b[p_b > 0]))

  if (h_a + h_b == 0) return(1)
  unname(2 * mi / (h_a + h_b))
}


#' Variation of Information
#' @noRd
.variation_of_information <- function(a, b) {
  n <- length(a)
  if (n == 0) return(0)

  tab <- table(a, b)
  p_ab <- tab / n
  p_a <- rowSums(tab) / n
  p_b <- colSums(tab) / n

  h_a <- -sum(p_a[p_a > 0] * log(p_a[p_a > 0]))
  h_b <- -sum(p_b[p_b > 0] * log(p_b[p_b > 0]))

  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (p_ab[i, j] > 0) {
        mi <- mi + p_ab[i, j] * log(p_ab[i, j] / (p_a[i] * p_b[j]))
      }
    }
  }

  unname(h_a + h_b - 2 * mi)
}


# --- S3 methods ---

#' @export
print.tessera_divergence <- function(x, ...) {
  cat("TESSERA Regime Divergence\n")
  cat(sprintf("  Method: %s | Reference: %s | Bins: %d\n",
              x$method, x$reference, x$n_bins))
  if (is.matrix(x$scores)) {
    cat("  Pairwise divergence matrix:\n")
    print(round(x$scores, 4))
  } else {
    cat("  Scores:\n")
    for (i in seq_along(x$scores)) {
      cat(sprintf("    %s: %.4f\n", names(x$scores)[i], x$scores[i]))
    }
    cat(sprintf("  Mean RD: %.4f | Max RD: %.4f\n",
                mean(x$scores), max(x$scores)))
  }
  invisible(x)
}


#' @export
summary.tessera_divergence <- function(object, ...) {
  print(object)
  if (!is.null(object$entity_divergence)) {
    ed <- object$entity_divergence[object$entity_divergence > 0]
    if (length(ed) > 0) {
      cat(sprintf("\nEntity divergence (non-zero): %d/%d entities (%.1f%%)\n",
                  length(ed), length(object$entity_divergence),
                  100 * length(ed) / length(object$entity_divergence)))
      cat(sprintf("  Mean: %.4f | Max: %.4f\n", mean(ed), max(ed)))
    }
  }
  invisible(object)
}


#' @export
plot.tessera_divergence <- function(x, ...) {
  if (is.matrix(x$scores)) {
    # Heatmap for all-pairs
    image(x$scores,
          main = sprintf("Regime Divergence (%s)", x$method),
          xlab = "Bin", ylab = "Bin",
          col = grDevices::hcl.colors(50, "Reds"))
  } else {
    # Line/bar plot for sequential divergence
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      barplot(x$scores,
              main = sprintf("Regime Divergence (%s, %s)", x$method, x$reference),
              ylab = "Divergence",
              col = "#D55E00", border = NA)
      return(invisible(x))
    }

    df <- data.frame(
      transition = factor(names(x$scores), levels = names(x$scores)),
      divergence = x$scores
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$transition,
                                           y = .data$divergence)) +
      ggplot2::geom_col(fill = "#D55E00", alpha = 0.85) +
      ggplot2::geom_point(size = 3, color = "#D55E00") +
      ggplot2::labs(x = "Bin Transition", y = "Regime Divergence",
                    title = sprintf("Regime Divergence (%s)", x$method)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    print(p)
  }
  invisible(x)
}
