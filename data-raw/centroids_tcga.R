## code to prepare `centroids_tcga` dataset goes here
library(org.Hs.eg.db)
library(BLCAsubtyping)
library(dplyr)
library(EnsDb.Hsapiens.v86)

og_centroids <- BLCAsubtyping:::tcga.centroids |>
  as_tibble()

colnames(og_centroids) <- c(
  "hgnc", "entrez", "Basal_squamous", "Luminal", "Luminal_infiltrated",
  "Luminal_papillary", "Neuronal"
)

centroids_tcga <- dplyr::relocate(og_centroids, hgnc, entrez, .after = last_col())

# The supplied IDs are entrez - converting to Ensembl is essentially trying to
# go from less specific to more specific - I can't spin information out of
# nothing.

# Instead, I'll plan to convert any provided Ensembl IDs down to entrez IDs.

usethis::use_data(centroids_tcga, overwrite = TRUE)
