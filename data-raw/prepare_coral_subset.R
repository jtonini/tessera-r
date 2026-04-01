# data-raw/prepare_coral_subset.R
#
# Prepares a representative subset of the Global Coral Bleaching Database
# (van Woesik et al. 2022) for inclusion as example data in the tessera
# R package.
#
# Source: van Woesik, R. et al. (2022). Global Coral Bleaching Database.
#         figshare. https://doi.org/10.6084/m9.figshare.20078993
#         License: CC-BY 4.0
#
# Run once on badenpowell:
#   cd ~/tessera-r
#   Rscript data-raw/prepare_coral_subset.R
#
# Output: data/coral_example.rda

library(DBI)
library(RSQLite)

# --- Configuration ---
SQLITE_PATH <- "~/tessera/data/Global_Coral_Bleaching_Database_SQLite_11_24_21.db"
TARGET_N <- 500L
SEED <- 42L

stopifnot(file.exists(SQLITE_PATH))
cat("Reading from:", SQLITE_PATH, "\n")

con <- dbConnect(SQLite(), SQLITE_PATH)

# --- Join the three key tables on Sample_ID ---
query <- "
SELECT
  se.Sample_ID,
  se.Date_Year   AS year,
  se.Depth_m     AS depth_m,
  e.ClimSST      AS clim_sst,
  e.Temperature_Mean AS sst_mean,
  e.SSTA         AS ssta,
  e.SSTA_DHW     AS ssta_dhw,
  e.SSTA_Frequency AS ssta_freq,
  e.Windspeed    AS windspeed,
  e.TSA          AS tsa,
  e.TSA_DHW      AS tsa_dhw,
  b.Percent_Bleached AS percent_bleached
FROM Sample_Event_tbl se
INNER JOIN Environmental_tbl e ON se.Sample_ID = e.Sample_ID
INNER JOIN Bleaching_tbl b     ON se.Sample_ID = b.Sample_ID
WHERE se.Date_Year IS NOT NULL
  AND se.Date_Year > 0
"

full_data <- dbGetQuery(con, query)
dbDisconnect(con)

cat("Joined dataset:", nrow(full_data), "rows,", ncol(full_data), "columns\n")
cat("Year range:", range(full_data$year, na.rm = TRUE), "\n")

# --- Select feature columns (exclude IDs, year, response) ---
feature_cols <- c("depth_m", "clim_sst", "sst_mean", "ssta", "ssta_dhw",
                  "ssta_freq", "windspeed", "tsa", "tsa_dhw")

# --- Remove rows with NAs in features ---
complete_mask <- complete.cases(full_data[, feature_cols])
cat("Complete cases:", sum(complete_mask), "of", nrow(full_data), "\n")
eligible <- full_data[complete_mask, ]

# --- Stratified sampling ---
set.seed(SEED)

# Binary bleaching indicator for stratification
eligible$bleached <- as.integer(!is.na(eligible$percent_bleached) &
                                  eligible$percent_bleached > 0)

# Temporal quartiles
breaks <- quantile(eligible$year, probs = seq(0, 1, length.out = 5))
eligible$time_bin <- cut(eligible$year, breaks = breaks,
                         include.lowest = TRUE, labels = FALSE)

strata <- interaction(eligible$time_bin, eligible$bleached)
strata_table <- table(strata)
cat("\nStrata distribution:\n")
print(strata_table)

# Sample per stratum
n_strata <- length(strata_table[strata_table > 0])
per_stratum <- ceiling(TARGET_N / n_strata)

sampled_idx <- integer(0)
for (s in names(strata_table)) {
  s_idx <- which(strata == s)
  if (length(s_idx) == 0) next
  n_take <- min(per_stratum, length(s_idx))
  sampled_idx <- c(sampled_idx, sample(s_idx, n_take))
}

if (length(sampled_idx) > TARGET_N) {
  sampled_idx <- sample(sampled_idx, TARGET_N)
}

# --- Build final dataset ---
keep_cols <- c("year", feature_cols, "percent_bleached")
coral_example <- eligible[sampled_idx, keep_cols]
rownames(coral_example) <- NULL

cat("\n--- Final dataset ---\n")
cat("Rows:", nrow(coral_example), "\n")
cat("Columns:", ncol(coral_example), "\n")
cat("Names:", paste(names(coral_example), collapse = ", "), "\n")
cat("Year range:", range(coral_example$year), "\n")
cat("Bleaching rate:", round(mean(coral_example$percent_bleached > 0,
                                   na.rm = TRUE) * 100, 1), "%\n")
cat("Size:", format(object.size(coral_example), units = "KB"), "\n")

str(coral_example)

# --- Save ---
dir.create("data", showWarnings = FALSE)
save(coral_example, file = "data/coral_example.rda", compress = "xz")
cat("\nSaved to data/coral_example.rda\n")
