#' Calculate centroid correlation for samples
#'
#' @param df A `data.frame` containing RNA expression with unique genes as rows
#'   and samples as columns. RNA-seq data must be log-transformed.
#' @param gene_id Character specifying the type of gene identifiers used for the
#'   row names or first column of `df`.
#' @param classifier Classifier to be used.
#' @param tidy Logical. If TRUE, assumes the first column contains the gene
#'   identifiers. Otherwise, assumes IDs are row names
#'
#' @details If using the TCGA classifier, you must provide either Entrez IDs or
#'   HGNC symbols. The original centroids provided for the TCGA classifier only
#'   included Entrez and HGNC, and the conversion from Entrez to Ensembl is
#'   ambiguous. Therefore, conversion and selection of unique genes must be
#'   performed by the user on a case-by-case basis, perhaps using the expression
#'   value of the gene as a guide for which genes to keep.
#'
#' @return A `tibble` containing: \describe{ \item{estimate}{Pearson correlation
#'   of a given `sample` to the given `class` centroid} \item{conf.low,
#'   conf.high}{low and high end of 95% confidence interval}
#'   \item{nearest}{centroid for which the sample has the highest correlation
#'   to} \item{statistic}{t-statistic} \item{parameter}{degrees of freedom}
#'   \item{sep_lvl}{(Highest Correlation - 2nd Highest
#'   Correlation)/median(Distance to highest correlation)} }
#' @author Aurelie Kamoun
#' @importFrom rlang .data
#' @export

classify_blca <- function(df,
                          gene_id = c("entrez", "ensembl", "hgnc"),
                          classifier = c("tcga", "uromol2021", "consensus"),
                          tidy = FALSE) {

  gene_id <- rlang::arg_match(gene_id)
  classifier <- rlang::arg_match(classifier)

  if (gene_id == "ensembl" && classifier == "tcga") {
    rlang::abort(
      "Ensembl IDs cannot be used with the TCGA classifier.",
      "Convert to either unique HGNC symbols or Entrez IDs."
    )
  }

  if (classifier == "tcga") {
    centroids <- classifyBLCA::centroids_tcga
  } else if (classifier == "uromol2021") {
    centroids <- classifyBLCA::centroids_uromol2021
  } else {
    centroids <- classifyBLCA::centroids_consensus
  }

  # If tidy, assume first col is gene_id. Make it row names and drop the col.
  if (tidy) {
    df <- as.data.frame(df)
    rownames(df) <- df[, 1]
    df <- df[, -1]
  }

  # Order genes the same in dataset as in centroids. Drop non-centroid genes.
  y <- df[match(centroids[[gene_id]], rownames(df)), ]

  if (nrow(y) == 0) {
    rlang::abort(
      c("Empty intersection between profiled genes and genes used for classification",
        "Ensure that gene names correspond to identifier specified by gene_id"),
    )
  }

  if (nrow(y) < 0.6 * nrow(centroids)) {
    warning(
      "Input data includes <60% of genes used for classification. Results may not be relevant"
    )
  }

  array_cor_test <- function(df, class) {
    df |>
      apply(2, \(x) stats::cor.test(x, centroids[[class]]) |> broom::tidy()) |>
      dplyr::bind_rows() |>
      dplyr::mutate(sample = colnames(df), class = class)
  }

  classes <- setdiff(colnames(centroids), c("hgnc", "ensembl", "entrez"))

  cor_res <-
    lapply(classes, array_cor_test, df = y) |>
    dplyr::bind_rows() |>
    dplyr::group_by(.data$sample) |>
    dplyr::arrange(dplyr::desc(.data$estimate)) |>
    dplyr::mutate(
      nearest = .data$class[1],
      nearest_cor = .data$estimate[1],
      dist_to_sec = .data$nearest_cor - .data$estimate[2],
      delta_med = stats::median(.data$nearest_cor - .data$estimate),
      sep_lvl = .data$dist_to_sec / .data$delta_med
    ) |>
    dplyr::ungroup() |>
    dplyr::select(c("sample", "class", "estimate", "conf.low", "conf.high", "nearest", "statistic", "parameter", "sep_lvl")) |>
    dplyr::arrange(.data$sample, .data$class)

  cor_res
}
