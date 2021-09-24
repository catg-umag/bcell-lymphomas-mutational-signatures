library(cluster)
library(deconstructSigs)
library(nnls)
library(tidyverse)


c_cosmic_signatures_urls <- "https://raw.githubusercontent.com/CATG-UMAG/bcell-lymphomas-mutational-signatures/main/data/cosmic_signatures_urls.csv"


#' Makes mutation matrix with motiffs in rows and some group in columns
#'
#' @param data Dataframe containing list of mutations
#' @param by Column containing the groups to separate the data
#' @param normalize Make each column sum 1
make_mutation_table <- function(data, by, normalize = TRUE) {
  by <- dplyr::enquo(by)

  table <- data %>%
    dplyr::count(substitution, context, !!by)

  if (normalize) {
    table <- table %>%
      dplyr::group_by(!!by) %>%
      dplyr::mutate(n = n / sum(n))
  }

  table <- table %>%
    tidyr::spread(!!by, n) %>%
    as.data.frame() %>%
    replace(is.na(.), 0) %>%
    dplyr::arrange(substitution, context)

  return(table)
}


#' Gets COSMIC signatures from internet
#'
#' @param version COSMIC version (2.0, 3.0, 3.1, 3.2)
#' @param genome Genome (GRCh37, GRCh38, mm9, mm10, rn6)
#' @param type Signature type (SBS, DBS, ID)
get_cosmic_signatures <- function(version = 3.2, genome = "GRCh38", type = "SBS") {
  urls_df <- read.csv(c_cosmic_signatures_urls)

  selected_signature <- urls_df %>%
    dplyr::filter(version == !!version & genome == !!genome & type == !!type)

  if (nrow(selected_signature) == 1) {
    signatures <- read.table(selected_signature$url, sep = "\t", header = TRUE) %>%
      dplyr::select(where(function(x) any(!is.na(x)))) %>%
      dplyr::rename(mutation = Type) %>%
      dplyr::mutate(X = sub("^(.)\\[(.>.)\\](.)$", "\\2 \\1.\\3", mutation)) %>%
      tidyr::separate(X, c("substitution", "context"), " ") %>%
      dplyr::select(substitution, context, mutation, dplyr::everything()) %>%
      dplyr::arrange(substitution, context)
  } else {
    stop(glue::glue("Could not find signatures for version={version}, genome={genome}, type={type}"))
  }

  return(signatures)
}


#' Formats substitution and context to X[X>Y]X format
#'
#' @param substitution Substitution (format X>Y)
#' @param context Substitution context (format X.X or XXX)
format_mutation <- function(substitution, context) {
  mutation <- gsub("(.)>(.) (.).(.)", "\\3[\\1>\\2]\\4", paste(substitution, context))
  return(mutation)
}


#' Fits samples to a group of signature using deconstructSigs
#'
#' @param samples Dataframe containing signatures to recreate with fitting columns: [mutation, samples...]
#' @param reference_signatures Dataframe containing reference signatures columns: [mutation, samples...]
fit_signatures <- function(samples, reference_signatures) {
  #' used to prepare the correct input for deconstructSigs
  prepare_for_fitting <- function(df) {
    df <- df %>%
      tibble::column_to_rownames("mutation") %>%
      dplyr::select(-substitution, -context) %>%
      t() %>%
      as.data.frame()

    return(df)
  }

  samples_table <- prepare_for_fitting(samples)
  signatures_table <- prepare_for_fitting(reference_signatures)

  weights <- do.call(rbind, lapply(rownames(samples_table), function(x) {
    ws <- whichSignatures(
      tumor.ref = samples_table,
      signatures.ref = signatures_table,
      signature.cutoff = 0.0,
      sample.id = x
    )
    return(ws$weights)
  }))

  table <- weights %>%
    dplyr::select_if(~ sum(.) != 0) %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather("signature", "contribution", -sample)

  return(table)
}


#' Hierarchical clustering of samples fitted using DIANA
#'
#' @param data Dataframe containing fitting information columns: [sample, signature, contribution]
fitting_clustering <- function(data) {
  fitting_matrix <- data %>%
    tidyr::spread(signature, contribution) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix()

  result <- cluster::diana(fitting_matrix, stand = FALSE)

  return(result)
}


#' Use a NNLS optimization to reconstruct a signature using combination of reference signatures
#'
#' @param signature Signature to reconstruct
#' @param reference_signatures Dataframe with reference signatures (matrix format)
#' @param n Number of signatures to combine for the reconstruction
reconstruct_signatures <- function(signature, reference_signatures, n = 2) {
  combinations <- combn(
    names(signatures_cosmic) %>%
      purrr::discard(~ . %in% c("substitution", "context", "mutation")),
    n
  )

  results <- apply(combinations, n, function(x) {
    selected_signatures <- as.matrix(reference_signatures[x])
    # solve NNLS
    result <- nnls(selected_signatures, signature)
    # get coefficients and normalize them (cosine similarity won't change)
    coefficients <- result$x / sum(result$x)

    name <- paste(x, collapse = "+")
    proportions <- paste(sapply(1:length(x), function(y) {
      glue::glue("{x[y]}={round(coefficients[y], 2)}")
    }), collapse = ";")
    # get cosine similarity for original and reconstructed signature
    similarity <- cosine_similarity(as.numeric(selected_signatures %*% coefficients), signature)
    return(c(name = name, proportions = proportions, similarity = similarity))
  })

  df <- data.frame(t(results)) %>%
    dplyr::arrange(desc(similarity))

  return(df)
}


#' Calculates cosine distance between two vectors
cosine_similarity <- function(x, y) {
  return(as.numeric(x %*% y / sqrt(x %*% x * y %*% y)))
}


#' Loads value from environment variable, with a default value
get_var <- function(name, default) {
  return(if (Sys.getenv(name) != "") {
    Sys.getenv(name)
  } else {
    default
  })
}