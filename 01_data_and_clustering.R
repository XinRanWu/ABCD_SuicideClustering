# Manuscript: sample, 34 clustering variables and K-means.
# Run from the repository root; data must remain outside this public repository.
data_dir <- Sys.getenv("ABCD_DATA_DIR", unset = "../ABCD")
prepared_dir <- Sys.getenv("ABCD_PREPARED_DIR", unset = "../Data_and_Scripts")
prepared_file <- file.path(prepared_dir, "abcd_SuicideData.RData")
if (!file.exists(prepared_file)) stop("Set ABCD_PREPARED_DIR to the directory containing abcd_SuicideData.RData")
load(prepared_file)
stopifnot(ncol(abcd_trait.bl.clust) >= 35, !anyDuplicated(abcd_trait.bl.clust$subjectkey))
# This saved object contains the original group labels; do not overwrite it when reproducing tables.
print(table(abcd_trait.bl.clust$clust, useNA="ifany"))
# Original clustering call, extracted from abcd_Suicide_Part1.R.
# Requires factoextra; randomness and K selection must be reconciled with the original run.
# Run the following block interactively when recalculation is required.
if (identical(Sys.getenv("RUN_CLUSTER"), "1")) {
  if (!requireNamespace("factoextra", quietly=TRUE)) stop("Install factoextra first")
  clustering_new <- factoextra::eclust(scale(abcd_trait.bl.clust[, 2:35]), nboot=100, FUNcluster="kmeans", iter.max=50)
  print(table(clustering_new$cluster))
}
# Original split-half validation: K is selected by eclust separately in each half.
# First 812 and remaining participants follow the original row order.
if (identical(Sys.getenv("RUN_SPLIT_HALF"), "1")) {
  if (!requireNamespace("factoextra", quietly=TRUE)) stop("Install factoextra first")
  if (nrow(abcd_trait.bl.clust) != 1624) stop("Original split requires 1,624 rows")
  split_one <- factoextra::eclust(scale(abcd_trait.bl.clust[1:812,2:35]), nboot=100, FUNcluster="kmeans", iter.max=50)
  split_two <- factoextra::eclust(scale(abcd_trait.bl.clust[813:1624,2:35]), nboot=100, FUNcluster="kmeans", iter.max=50)
  print(cor(t(split_one$centers), t(split_two$centers)))
  # Align labels and verify results against the original supplementary table.
}
# Source: Data_and_Scripts/abcd_Suicide_Part1.R, sample construction lines 499-591.
# Input sample and labels here come from the prepared workspace, not a fresh ABCD download.
