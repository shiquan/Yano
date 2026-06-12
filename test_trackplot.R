# Test updated TrackPlot function
library(Yano)

# Load GTF database
message("Loading GTF database...")
gtf <- gtf2db("/projects/conco/qshi/shiquan.github.io/gencode.v44.annotation.gtf.gz")

bamfile <- "/projects/conco/qshi/shiquan.github.io/Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam"

# ---- Test 1: SERF2 gene with highlights (no cell.group) ----
message("Test 1: SERF2 with highlights, no cell.group")
png("/projects/conco/qshi/Yano/test_trackplot_SERF2.png", width = 12, height = 6, units = "in", res = 150)
TrackPlot(
  bamfile = bamfile,
  gtf = gtf,
  gene = "SERF2",
  highlights = c(43801711, 43804427),
  max.depth = 2,
  display.genes = "SERF2"
)
dev.off()
message("Saved test_trackplot_SERF2.png")

# ---- Test 2: MYL6 gene with junction and highlights (no cell.group) ----
message("Test 2: MYL6 with junction and highlights")
png("/projects/conco/qshi/Yano/test_trackplot_MYL6.png", width = 12, height = 8, units = "in", res = 150)
TrackPlot(
  bamfile = bamfile,
  gtf = gtf,
  gene = "MYL6",
  junc = TRUE,
  highlights = list(c(56160320, 56161387), c(56161387, 56161465))
)
dev.off()
message("Saved test_trackplot_MYL6.png")

# ---- Test 3: SERF2 with default max.depth (no capping) ----
message("Test 3: SERF2 default depth, no highlights")
png("/projects/conco/qshi/Yano/test_trackplot_SERF2_default.png", width = 12, height = 6, units = "in", res = 150)
TrackPlot(
  bamfile = bamfile,
  gtf = gtf,
  gene = "SERF2",
  display.genes = "SERF2"
)
dev.off()
message("Saved test_trackplot_SERF2_default.png")

message("All tests completed.")
