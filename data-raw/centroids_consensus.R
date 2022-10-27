## code to prepare `centroids_consensus` dataset goes here
library(consensusMIBC)
library(dplyr)

og_centroids <- consensusMIBC::centroids |>
  as_tibble()

colnames(og_centroids) <- c("lum_p", "lum_ns", "lum_u", "stroma_rich", "ba_sq", "ne_like", "hgnc", "entrez", "ensembl")

centroids_consensus <- relocate(og_centroids, ensembl, .before = "entrez")

usethis::use_data(centroids_consensus, overwrite = TRUE)
