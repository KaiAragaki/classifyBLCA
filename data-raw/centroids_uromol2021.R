library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(classifyNMIBC)
library(dplyr)

og_centroids <- classifyNMIBC::centroids |>
  as_tibble()

annot_eg <- AnnotationDbi::select(
  org.Hs.eg.db,
  columns = c("ENTREZID", "SYMBOL", "ENSEMBL"),
  keys = keys(org.Hs.eg.db)
)

annot_ez <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  columns = c("ENTREZID", "SYMBOL", "GENEID"),
  keys = keys(EnsDb.Hsapiens.v86)
)

# First try to map ensembl to entrez...
og_centroids_w_entrez <- left_join(og_centroids, annot_eg, by = c("ensembl_gene_ID" = "ENSEMBL"))

duplicated_ids <- og_centroids_w_entrez[duplicated(og_centroids_w_entrez$ensembl_gene_ID),]$ensembl_gene_ID
# There are 13 duplicate ids

duplicated_rows <- og_centroids_w_entrez[og_centroids_w_entrez$ensembl_gene_ID %in% duplicated_ids,]

# If there are duplicates, prefer those that match the hgnc of the og dataset:
rm_dupes <- og_centroids_w_entrez |>
  dplyr::filter(!ensembl_gene_ID %in% duplicated_ids | SYMBOL == hgnc_symbol)

still_missing <- rm_dupes[is.na(rm_dupes$ENTREZID),]
# 111 rows

# Try matching by provided symbols
by_provided_symbols <- left_join(still_missing, annot_eg, by = c("hgnc_symbol" = "SYMBOL"))

# There's one duplicate - CTAGE8
# The 'correct' ENSG appears to be ENSG00000289604
by_provided_symbols <- dplyr::filter(by_provided_symbols, hgnc_symbol != "CTAGE8" | ENSEMBL == "ENSG00000289604")

# Format to rejoin with old data
by_provided_symbols <- dplyr::select(by_provided_symbols, Class_1, Class_2a,
                                     Class_2b, Class_3, hgnc_symbol, ensembl_gene_ID,
                                     ENTREZID = ENTREZID.y, SYMBOL)

no_missing <- rm_dupes[!is.na(rm_dupes$ENTREZID),]

together <- rbind(no_missing, by_provided_symbols)

still_missing <- together[is.na(together$ENTREZID),]
# 5 rows

# These genes require manual curation :/
# https://useast.ensembl.org/Homo_sapiens/Gene/Idhistory?db=core;g=ENSG00000288825;r=1:149842218-149842750;t=ENST00000369159
still_missing[1, 5] <- "H2AC18"
still_missing[1, 6] <- "ENSG00000288825"
still_missing[1, 7] <- "8337"

still_missing[3, 7] <- "494150"

still_missing[4, 5] <- "SMC5-DT"
still_missing[4, 7] <- "100507299"

still_missing[5, 5] <- "H3C2"
still_missing[5, 6] <- "ENSG00000286522"
still_missing[5, 7] <- "8358"

# There are no current external IDs for LINC01336

no_missing <- together[!is.na(together$ENTREZID),]

final <- rbind(no_missing, still_missing)

centroids_uromol2021 <- dplyr::select(
  final, Class_1, Class_2a, Class_2b, Class_3,
  hgnc = hgnc_symbol, ensembl = ensembl_gene_ID,
  entrez = ENTREZID
)

usethis::use_data(centroids_uromol2021, overwrite = TRUE)
