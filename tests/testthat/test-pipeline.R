# tests/testthat/test-pipeline.R

test_that("tessera_characterize handles basic input", {
  set.seed(42)
  df <- data.frame(
    year = sample(2000:2020, 200, replace = TRUE),
    x1 = rnorm(200),
    x2 = rnorm(200),
    x3 = runif(200)
  )

  char <- tessera_characterize(df, time_col = "year", n_bins = 4)
  expect_s3_class(char, "tessera_characterize")
  expect_equal(length(char$features), 3)
  expect_true(nrow(char$bins) >= 1)
  expect_true(all(c(".bin") %in% names(char$data)))
  expect_equal(nrow(char$data), 200)
})

test_that("tessera_characterize uses quantile binning by default", {
  set.seed(42)
  # Skewed temporal distribution
  df <- data.frame(
    year = c(rep(2000, 100), rep(2010, 50), rep(2020, 50)),
    x1 = rnorm(200)
  )

  char <- tessera_characterize(df, time_col = "year", n_bins = 3)
  # Quantile bins should balance counts better than fixed-width
  counts <- char$bins$n_entities
  expect_true(max(counts) / min(counts) < 5)
})

test_that("tessera_characterize drops rows with NAs", {
  df <- data.frame(
    year = 1:100,
    x1 = c(rnorm(90), rep(NA, 10)),
    x2 = rnorm(100)
  )

  expect_message(
    char <- tessera_characterize(df, time_col = "year", n_bins = 2),
    "Dropping 10 rows"
  )
  expect_equal(nrow(char$data), 90)
})

test_that("tessera_characterize scales features within bins", {
  set.seed(42)
  df <- data.frame(
    year = rep(c(1, 2), each = 50),
    x1 = c(rnorm(50, 0, 1), rnorm(50, 100, 1))
  )

  char <- tessera_characterize(df, time_col = "year", n_bins = 2, scale = TRUE)
  # After within-bin scaling, each bin's features should be ~N(0,1)
  for (b in unique(char$data$.bin)) {
    vals <- char$data$x1[char$data$.bin == b]
    expect_equal(mean(vals), 0, tolerance = 1e-10)
    expect_equal(sd(vals), 1, tolerance = 1e-10)
  }
})


test_that("tessera_similarity computes multiple measures", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 50, ncol = 4)

  sim <- tessera_similarity(mat, measures = c("cosine", "euclidean"), k = 10)
  expect_s3_class(sim, "tessera_similarity")
  expect_equal(length(sim$matrices), 2)
  expect_true("cosine" %in% names(sim$matrices))
  expect_true("euclidean" %in% names(sim$matrices))
})

test_that("tessera_similarity returns sparse by default", {
  set.seed(42)
  mat <- matrix(rnorm(400), nrow = 100, ncol = 4)

  sim <- tessera_similarity(mat, measures = "cosine", k = 10)
  expect_true(inherits(sim$matrices$cosine, "sparseMatrix"))
})

test_that("similarity values are in expected range", {
  set.seed(42)
  mat <- matrix(runif(200), nrow = 50, ncol = 4)

  sim <- tessera_similarity(mat, measures = "cosine", k = 10, sparse = FALSE)
  vals <- sim$matrices$cosine[upper.tri(sim$matrices$cosine)]
  expect_true(all(vals >= 0 & vals <= 1))

  sim_e <- tessera_similarity(mat, measures = "euclidean", k = 10, sparse = FALSE)
  vals_e <- sim_e$matrices$euclidean[upper.tri(sim_e$matrices$euclidean)]
  expect_true(all(vals_e >= 0 & vals_e <= 1))
})


test_that("tessera_network builds from similarity", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 50, ncol = 4)
  sim <- tessera_similarity(mat, measures = "cosine", k = 10)

  net <- tessera_network(sim, method = "knn", k = 8)
  expect_s3_class(net, "tessera_network")
  expect_true(igraph::is_igraph(net$graph))
  expect_equal(igraph::vcount(net$graph), 50)
  expect_true(igraph::ecount(net$graph) > 0)
})

test_that("tessera_network threshold method works", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 50, ncol = 4)
  sim <- tessera_similarity(mat, measures = "cosine", k = 49, sparse = FALSE)

  net <- tessera_network(sim, method = "threshold", threshold = 0.8)
  expect_s3_class(net, "tessera_network")
  # All edge weights should be >= threshold
  w <- igraph::E(net$graph)$weight
  if (length(w) > 0) {
    expect_true(all(w >= 0.8))
  }
})


test_that("tessera_communities detects structure with louvain", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 50, ncol = 4)
  sim <- tessera_similarity(mat, measures = "cosine", k = 10)
  net <- tessera_network(sim, method = "knn", k = 8)

  com <- tessera_communities(net, algorithm = "louvain")
  expect_s3_class(com, "tessera_communities")
  expect_equal(length(com$membership), 50)
  expect_true(com$n_communities >= 1)
  expect_false(is.na(com$modularity))
})

test_that("community membership is capped at max_communities", {
  set.seed(42)
  # Large enough dataset to get many communities
  mat <- matrix(rnorm(2000), nrow = 200, ncol = 10)
  sim <- tessera_similarity(mat, measures = "cosine", k = 5)
  net <- tessera_network(sim, method = "knn", k = 5)

  com <- tessera_communities(net, algorithm = "louvain",
                              resolution = 5.0,  # force many communities
                              max_communities = 5L)
  expect_true(com$n_communities <= 5)
})


test_that("full pipeline runs with tessera() wrapper", {
  set.seed(42)
  df <- data.frame(
    year = sample(2000:2020, 200, replace = TRUE),
    x1 = rnorm(200),
    x2 = rnorm(200),
    x3 = runif(200)
  )

  result <- tessera(df, time_col = "year", verbose = FALSE)
  expect_s3_class(result, "tessera_result")
  expect_s3_class(result$characterize, "tessera_characterize")
  expect_s3_class(result$similarity, "tessera_similarity")
  expect_s3_class(result$network, "tessera_network")
  expect_s3_class(result$communities, "tessera_communities")
  expect_s3_class(result$divergence, "tessera_divergence")
})

test_that("ecology preset uses expected parameters", {
  set.seed(42)
  df <- data.frame(
    year = sample(2000:2020, 200, replace = TRUE),
    x1 = rnorm(200),
    x2 = rnorm(200)
  )

  # Use louvain fallback since leiden may not be installed
  result <- tessera(df, time_col = "year", preset = "ecology",
                     community_algorithm = "louvain",  # override leiden dep
                     verbose = FALSE)
  expect_equal(result$similarity$measures, c("bray_curtis", "simpson"))
})


test_that("NMI is 1 for identical partitions", {
  a <- c(1, 1, 2, 2, 3, 3)
  expect_equal(tessera:::.nmi(a, a), 1.0)
})

test_that("NMI is 0 for independent partitions", {
  # Approximately 0 for large enough independent partitions
  set.seed(42)
  a <- sample(1:5, 1000, replace = TRUE)
  b <- sample(1:5, 1000, replace = TRUE)
  expect_true(tessera:::.nmi(a, b) < 0.05)
})

test_that("variation of information is 0 for identical partitions", {
  a <- c(1, 1, 2, 2, 3, 3)
  expect_equal(tessera:::.variation_of_information(a, a), 0)
})
