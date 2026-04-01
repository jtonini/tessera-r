#' Characterize Entities via Temporal Binning
#'
#' Step 1 of the TESSERA pipeline. Assigns entities to temporal bins using
#' quantile-based or fixed-width boundaries, validates features, and prepares
#' the data structure for downstream similarity computation.
#'
#' @param data A data.frame or matrix of entity features. Rows are entities,
#'   columns are features.
#' @param time_col Character string naming the column that contains temporal
#'   information (e.g., year, timestamp). If `NULL`, all entities are treated
#'   as a single bin.
#' @param entity_col Character string naming the column that identifies
#'   individual entities. If `NULL`, row indices are used.
#' @param n_bins Integer. Number of temporal bins. Default is 4.
#' @param bin_method Character. Either `"quantile"` (default, recommended) or
#'   `"fixed"` for equal-width bins. Quantile-based binning ensures balanced
#'   entity counts across bins, which is critical for datasets with uneven
#'   temporal distributions.
#' @param feature_cols Character vector of column names to use as features.
#'   If `NULL`, all numeric columns except `time_col` and `entity_col` are used.
#' @param min_entities Integer. Minimum number of entities required per bin.
#'   Bins with fewer entities are merged with adjacent bins. Default is 30.
#' @param scale Logical. Whether to standardize features to zero mean and unit
#'   variance within each bin. Default is `TRUE`.
#'
#' @return A `tessera_characterize` object (S3) containing:
#' \describe{
#'   \item{data}{The processed data.frame with bin assignments}
#'   \item{bins}{A data.frame describing bin boundaries and entity counts}
#'   \item{features}{Character vector of feature column names used}
#'   \item{time_col}{Name of the time column}
#'   \item{entity_col}{Name of the entity column}
#'   \item{params}{List of parameters used}
#' }
#'
#' @examples
#' # Synthetic example
#' set.seed(42)
#' df <- data.frame(
#'   year = sample(2000:2020, 500, replace = TRUE),
#'   temp = rnorm(500, 25, 3),
#'   depth = runif(500, 0, 50),
#'   salinity = rnorm(500, 35, 1)
#' )
#' char <- tessera_characterize(df, time_col = "year", n_bins = 4)
#' summary(char)
#'
#' @export
tessera_characterize <- function(data,
                                  time_col = NULL,
                                  entity_col = NULL,
                                  n_bins = 4L,
                                  bin_method = c("quantile", "fixed"),
                                  feature_cols = NULL,
                                  min_entities = 30L,
                                  scale = TRUE) {
  bin_method <- match.arg(bin_method)

  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  n_bins <- as.integer(n_bins)
  if (n_bins < 1L) stop("`n_bins` must be >= 1.", call. = FALSE)


  # Identify feature columns
  exclude_cols <- c(time_col, entity_col)
  if (is.null(feature_cols)) {
    num_cols <- vapply(data, is.numeric, logical(1))
    feature_cols <- setdiff(names(data)[num_cols], exclude_cols)
  }
  if (length(feature_cols) == 0L) {
    stop("No numeric feature columns found.", call. = FALSE)
  }

  # Validate feature columns exist
  missing_cols <- setdiff(feature_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Feature columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Handle missing values in features
  feature_data <- data[, feature_cols, drop = FALSE]
  complete_rows <- complete.cases(feature_data)
  n_dropped <- sum(!complete_rows)
  if (n_dropped > 0L) {
    message(sprintf("Dropping %d rows with missing feature values (%.1f%%).",
                    n_dropped, 100 * n_dropped / nrow(data)))
    data <- data[complete_rows, , drop = FALSE]
  }

  # Entity IDs
  if (!is.null(entity_col)) {
    if (!entity_col %in% names(data)) {
      stop("Entity column '", entity_col, "' not found.", call. = FALSE)
    }
    entity_ids <- data[[entity_col]]
  } else {
    entity_ids <- seq_len(nrow(data))
    entity_col <- ".entity_id"
    data[[entity_col]] <- entity_ids
  }

  # Temporal binning
  if (!is.null(time_col)) {
    if (!time_col %in% names(data)) {
      stop("Time column '", time_col, "' not found.", call. = FALSE)
    }
    time_vals <- data[[time_col]]
    if (!is.numeric(time_vals)) {
      time_vals <- as.numeric(time_vals)
      if (any(is.na(time_vals))) {
        stop("Could not coerce `time_col` to numeric.", call. = FALSE)
      }
    }

    if (bin_method == "quantile") {
      breaks <- unique(stats::quantile(time_vals, probs = seq(0, 1, length.out = n_bins + 1L)))
      if (length(breaks) < 3L) {
        warning("Quantile binning produced fewer bins than requested. ",
                "Falling back to fixed-width bins.", call. = FALSE)
        breaks <- seq(min(time_vals), max(time_vals), length.out = n_bins + 1L)
      }
    } else {
      breaks <- seq(min(time_vals), max(time_vals), length.out = n_bins + 1L)
    }

    # Assign bins (right-closed, include lowest)
    data$.bin <- as.integer(cut(time_vals, breaks = breaks,
                                 include.lowest = TRUE, labels = FALSE))
    actual_n_bins <- max(data$.bin, na.rm = TRUE)

    # Merge small bins
    bin_counts <- table(data$.bin)
    small_bins <- as.integer(names(bin_counts)[bin_counts < min_entities])
    if (length(small_bins) > 0L && length(small_bins) < actual_n_bins) {
      for (sb in small_bins) {
        # Merge with adjacent bin (prefer right neighbor, fall back to left)
        if (sb < actual_n_bins) {
          data$.bin[data$.bin == sb] <- sb + 1L
        } else {
          data$.bin[data$.bin == sb] <- sb - 1L
        }
      }
      # Re-label bins sequentially
      old_bins <- sort(unique(data$.bin))
      new_labels <- seq_along(old_bins)
      data$.bin <- new_labels[match(data$.bin, old_bins)]
      message(sprintf("Merged %d small bins (< %d entities). Final bin count: %d.",
                      length(small_bins), min_entities, max(data$.bin)))
    }

    # Build bin summary
    bin_summary <- do.call(rbind, lapply(sort(unique(data$.bin)), function(b) {
      idx <- data$.bin == b
      data.frame(
        bin = b,
        n_entities = sum(idx),
        time_min = min(time_vals[idx]),
        time_max = max(time_vals[idx]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    # Single bin
    data$.bin <- 1L
    bin_summary <- data.frame(bin = 1L, n_entities = nrow(data),
                              time_min = NA_real_, time_max = NA_real_,
                              stringsAsFactors = FALSE)
  }

  # Scale features within bins
  if (scale) {
    for (b in unique(data$.bin)) {
      idx <- data$.bin == b
      for (fc in feature_cols) {
        vals <- data[idx, fc]
        s <- sd(vals, na.rm = TRUE)
        if (!is.na(s) && s > 0) {
          data[idx, fc] <- (vals - mean(vals, na.rm = TRUE)) / s
        } else {
          data[idx, fc] <- 0
        }
      }
    }
  }

  structure(
    list(
      data = data,
      bins = bin_summary,
      features = feature_cols,
      time_col = time_col,
      entity_col = entity_col,
      params = list(n_bins = n_bins, bin_method = bin_method,
                    min_entities = min_entities, scale = scale,
                    n_dropped = n_dropped)
    ),
    class = "tessera_characterize"
  )
}


#' @export
print.tessera_characterize <- function(x, ...) {
  cat("TESSERA Characterization\n")
  cat(sprintf("  Entities: %d across %d bin(s)\n",
              nrow(x$data), nrow(x$bins)))
  cat(sprintf("  Features: %d (%s)\n",
              length(x$features),
              if (x$params$scale) "scaled" else "raw"))
  if (x$params$n_dropped > 0L) {
    cat(sprintf("  Dropped: %d rows with missing values\n", x$params$n_dropped))
  }
  cat(sprintf("  Binning: %s (%d requested)\n",
              x$params$bin_method, x$params$n_bins))
  invisible(x)
}


#' @export
summary.tessera_characterize <- function(object, ...) {
  cat("TESSERA Characterization Summary\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  cat(sprintf("Total entities: %d\n", nrow(object$data)))
  cat(sprintf("Features (%d): %s\n",
              length(object$features),
              paste(object$features, collapse = ", ")))
  cat(sprintf("Scaling: %s\n", if (object$params$scale) "yes" else "no"))
  cat("\nBin distribution:\n")
  print(object$bins, row.names = FALSE)
  invisible(object)
}


#' @export
plot.tessera_characterize <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    barplot(x$bins$n_entities, names.arg = x$bins$bin,
            xlab = "Bin", ylab = "Entity count",
            main = "Entity Distribution Across Temporal Bins",
            col = "#4A90D9")
    return(invisible(x))
  }

  p <- ggplot2::ggplot(x$bins, ggplot2::aes(x = factor(.data$bin),
                                              y = .data$n_entities)) +
    ggplot2::geom_col(fill = "#4A90D9", alpha = 0.8) +
    ggplot2::labs(x = "Temporal Bin", y = "Entity Count",
                  title = "Entity Distribution Across Temporal Bins") +
    ggplot2::theme_minimal()
  print(p)
  invisible(x)
}
