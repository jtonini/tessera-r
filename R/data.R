#' Coral Bleaching Example Dataset
#'
#' A curated subset of the Global Coral Bleaching Database (van Woesik et al.
#' 2022) for demonstrating the TESSERA framework. Contains approximately 500
#' records stratified across temporal quartiles and bleaching outcomes, ensuring
#' balanced representation for pipeline demonstration.
#'
#' @format A data.frame with approximately 500 rows and columns including:
#' \describe{
#'   \item{year}{Year of observation (1983–2019)}
#'   \item{sst}{Sea surface temperature}
#'   \item{ssta}{Sea surface temperature anomaly}
#'   \item{depth}{Reef depth (meters)}
#'   \item{distance_to_shore}{Distance from shore (km)}
#'   \item{turbidity}{Water turbidity}
#'   \item{cyclone_frequency}{Cyclone frequency in the region}
#'   \item{windspeed}{Wind speed}
#'   \item{climsst}{Climatological sea surface temperature}
#'   \item{percent_bleached}{Percentage of coral bleached (response variable)}
#' }
#'
#' Column names are normalized to lowercase with underscores. The exact set of
#' columns depends on the source database version; use `names(coral_example)`
#' to inspect.
#'
#' @details
#' The full Global Coral Bleaching Database contains 31,199 records across 28
#' features spanning 1983–2019, with a 24.2\% bleaching rate. This subset
#' preserves the temporal and outcome distribution through stratified sampling
#' (seed = 42).
#'
#' @section Regenerating this dataset:
#' The script `data-raw/prepare_coral_subset.R` documents the exact extraction
#' procedure. It requires the original SQLite database from figshare.
#'
#' @source van Woesik, R. et al. (2022). Global Coral Bleaching Database.
#'   figshare. \doi{10.6084/m9.figshare.20078993}. License: CC-BY 4.0.
#'
#' @examples
#' \dontrun{
#' data(coral_example)
#'
#' # Quick look
#' str(coral_example)
#' summary(coral_example)
#'
#' # Run TESSERA pipeline
#' result <- tessera(coral_example,
#'                   time_col = "year",
#'                   preset = "ecology")
#' result
#' plot(result)
#' }
"coral_example"
