# paper/generate_figures.R
#
# Generates all figures for the MEE tessera paper.
# Run from the tessera-r directory:
#   Rscript paper/generate_figures.R
#
# Requires: tessera installed, full coral database accessible

library(DBI)
library(RSQLite)
library(tessera)
library(ggplot2)
library(igraph)

# --- Okabe-Ito colorblind palette ---
oi_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999",
                "#000000", "#F5C710", "#1A85FF", "#D41159")

outdir <- "paper/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================
# Load full coral dataset
# =============================================================
cat("Loading coral data...\n")
SQLITE_PATH <- "~/tessera/data/Global_Coral_Bleaching_Database_SQLite_11_24_21.db"
con <- dbConnect(SQLite(), SQLITE_PATH)

# Only select features — no response variable, no IDs
query <- "
SELECT
  se.Date_Year       AS year,
  se.Depth_m         AS depth_m,
  e.ClimSST          AS clim_sst,
  e.Temperature_Mean AS sst_mean,
  e.SSTA             AS ssta,
  e.SSTA_DHW         AS ssta_dhw,
  e.SSTA_Frequency   AS ssta_freq,
  e.Windspeed        AS windspeed,
  e.TSA              AS tsa,
  e.TSA_DHW          AS tsa_dhw
FROM Sample_Event_tbl se
INNER JOIN Environmental_tbl e ON se.Sample_ID = e.Sample_ID
WHERE se.Date_Year IS NOT NULL
  AND se.Date_Year > 0
"
coral_full <- dbGetQuery(con, query)
dbDisconnect(con)
cat("Full dataset:", nrow(coral_full), "rows\n")

# Remove NAs in features
feature_cols <- c("depth_m", "clim_sst", "sst_mean", "ssta", "ssta_dhw",
                  "ssta_freq", "windspeed", "tsa", "tsa_dhw")
complete_mask <- complete.cases(coral_full[, feature_cols])
coral_clean <- coral_full[complete_mask, ]
cat("After NA removal:", nrow(coral_clean), "rows\n")
cat("Year range:", range(coral_clean$year), "\n")

# =============================================================
# Run TESSERA pipeline
# =============================================================
cat("\nRunning TESSERA pipeline...\n")
t0 <- proc.time()

result <- tessera(coral_clean,
                  time_col = "year",
                  measures = c("bray_curtis", "simpson"),
                  n_bins = 4,
                  bin_method = "quantile",
                  k = 15,
                  network_method = "knn",
                  community_algorithm = "louvain",
                  resolution = 1.0,
                  divergence_method = "nmi",
                  divergence_reference = "adjacent",
                  verbose = TRUE)

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\nPipeline complete in %.1f seconds.\n", elapsed))

# Print summary for manuscript
cat("\n=== RESULTS FOR MANUSCRIPT ===\n")
print(result)
cat("\n--- Bin details ---\n")
print(result$characterize$bins)
cat("\n--- Divergence scores ---\n")
print(result$divergence)
cat("\n--- Community counts per bin ---\n")
for (bname in names(result$communities$by_bin)) {
  bc <- result$communities$by_bin[[bname]]
  cat(sprintf("  %s: %d communities, modularity=%.4f\n",
              bname, bc$n_communities, bc$modularity))
}

# =============================================================
# FIGURE 1: Pipeline diagram
# =============================================================
cat("\nGenerating Figure 1 (pipeline)...\n")

pdf(file.path(outdir, "fig1_pipeline.pdf"), width = 10, height = 3)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 2))

boxes <- data.frame(
  x = c(1, 3, 5, 7, 9),
  label = c("Characterize", "Similarity", "Network", "Communities", "Divergence"),
  func = c("tessera_characterize()", "tessera_similarity()", "tessera_network()",
           "tessera_communities()", "tessera_divergence()"),
  step = c("Step 1", "Step 2a", "Step 2b", "Step 3a", "Step 3b"),
  col = c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00")
)

for (i in seq_len(nrow(boxes))) {
  rect(boxes$x[i] - 0.8, 0.4, boxes$x[i] + 0.8, 1.6,
       col = adjustcolor(boxes$col[i], alpha.f = 0.2),
       border = boxes$col[i], lwd = 2)
  text(boxes$x[i], 1.2, boxes$label[i], cex = 0.9, font = 2)
  text(boxes$x[i], 0.85, boxes$func[i], cex = 0.55, family = "mono")
  text(boxes$x[i], 0.55, boxes$step[i], cex = 0.6, col = "grey40")
  if (i < nrow(boxes)) {
    arrows(boxes$x[i] + 0.85, 1.0, boxes$x[i + 1] - 0.85, 1.0,
           length = 0.1, lwd = 2, col = "grey40")
  }
}

text(5, 1.85, "tessera() convenience wrapper", cex = 0.8, font = 3, col = "grey30")
rect(0.05, 0.25, 9.95, 1.75, border = "grey70", lty = 2, lwd = 1)
dev.off()
cat("  -> fig1_pipeline.pdf\n")


# =============================================================
# FIGURE 2: Regime divergence across temporal bins
# =============================================================
cat("Generating Figure 2 (divergence)...\n")

div <- result$divergence
div_df <- data.frame(
  transition = factor(names(div$scores), levels = names(div$scores)),
  divergence = div$scores
)

bins <- result$characterize$bins

p2 <- ggplot(div_df, aes(x = transition, y = divergence)) +
  geom_col(fill = "#D55E00", alpha = 0.85, width = 0.6) +
  geom_point(size = 3, color = "#D55E00") +
  labs(x = "Bin Transition",
       y = "Regime Divergence (1 - NMI)",
       title = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(outdir, "fig2_divergence.pdf"), p2,
       width = 7, height = 4.5, dpi = 300)
ggsave(file.path(outdir, "fig2_divergence.png"), p2,
       width = 7, height = 4.5, dpi = 300)
cat("  -> fig2_divergence.pdf/.png\n")


# =============================================================
# FIGURE 3: Community structure for two contrasting bins
# =============================================================
cat("Generating Figure 3 (communities)...\n")

bin_names <- names(result$communities$by_bin)

pdf(file.path(outdir, "fig3_communities.pdf"), width = 12, height = 5.5)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

for (i in c(1, length(bin_names))) {
  bname <- bin_names[i]
  bc <- result$communities$by_bin[[bname]]
  g <- bc$graph
  n_com <- bc$n_communities
  palette <- oi_palette[seq_len(min(n_com, length(oi_palette)))]

  if (vcount(g) < 500) {
    set.seed(42)
    l <- layout_with_fr(g, niter = 500)
  } else {
    set.seed(42)
    l <- layout_with_drl(g)
  }

  w <- E(g)$weight
  w_scaled <- 0.3 + 1.5 * (w - min(w)) / max(max(w) - min(w), 1e-10)

  bin_idx <- as.integer(gsub("bin_", "", bname))
  yr_min <- bins$time_min[bin_idx]
  yr_max <- bins$time_max[bin_idx]

  plot(g, layout = l,
       vertex.size = 2,
       vertex.label = NA,
       vertex.color = palette[bc$membership],
       vertex.frame.color = NA,
       edge.width = w_scaled,
       edge.color = adjustcolor("grey60", alpha.f = 0.15),
       main = paste0("(", c("a", "b")[ifelse(i == 1, 1, 2)], ") ",
                     yr_min, "-", yr_max, "\n",
                     n_com, " communities, Q=",
                     round(bc$modularity, 3)))
}
dev.off()
cat("  -> fig3_communities.pdf\n")


# =============================================================
# FIGURE 4: Per-entity divergence distribution
# =============================================================
cat("Generating Figure 4 (entity divergence)...\n")

if (!is.null(div$entity_divergence)) {
  ed <- div$entity_divergence
  ed_df <- data.frame(divergence = ed[ed > 0])

  p4 <- ggplot(ed_df, aes(x = divergence)) +
    geom_histogram(fill = "#0072B2", alpha = 0.8, bins = 30,
                   color = "white", linewidth = 0.3) +
    geom_vline(xintercept = mean(ed_df$divergence), linetype = "dashed",
               color = "#D55E00", linewidth = 0.8) +
    annotate("text", x = mean(ed_df$divergence) + 0.02,
             y = Inf, vjust = 2, hjust = 0,
             label = paste0("mean = ", round(mean(ed_df$divergence), 3)),
             color = "#D55E00", size = 4) +
    labs(x = "Per-Entity Divergence Score",
         y = "Count",
         title = NULL) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(outdir, "fig4_entity_divergence.pdf"), p4,
         width = 6, height = 4, dpi = 300)
  ggsave(file.path(outdir, "fig4_entity_divergence.png"), p4,
         width = 6, height = 4, dpi = 300)
  cat("  -> fig4_entity_divergence.pdf/.png\n")

  cat(sprintf("\n  Entity divergence stats:\n"))
  cat(sprintf("    Total entities: %d\n", length(ed)))
  cat(sprintf("    Non-zero: %d (%.1f%%)\n",
              sum(ed > 0), 100 * mean(ed > 0)))
  cat(sprintf("    Mean (non-zero): %.4f\n", mean(ed[ed > 0])))
  cat(sprintf("    Max: %.4f\n", max(ed)))
} else {
  cat("  Per-entity divergence not available.\n")
}

cat("\n=== All figures generated in", outdir, "===\n")
