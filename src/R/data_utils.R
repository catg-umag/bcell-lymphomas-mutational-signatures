library(dplyr)
library(SomaticCancerAlterations)


#' Calculates the confidence interval from standard deviation (sd) and number of objects (n)
ci <- function(vsd, vn, c = 0.95) {
  return(qt(1 - ((1 - c) / 2), vn - 1) * vsd / sqrt(vn))
}


#' Calculates cosine distance between two vectors
cosine_dist <- function(x, y) {
  return(1 - as.numeric(x %*% y / sqrt(x %*% x * y %*% y)))
}


#' Reads a dataframe with variants data and produces counts for each kind of mutation
get_mut_type_occurrences <- function(data) {
  mutation_order <- c(
    "C>T at CpG", "C>T other", "C>T",
    "T>C", "C>A", "C>G", "T>A", "T>G"
  )

  table <- data %>%
    dplyr::mutate(sub2 = ifelse(substitution == "C>T",
      ifelse(grepl("..G", context), "C>T at CpG", "C>T other"),
      substitution
    )) %>%
    dplyr::count(lymph, sample, sub2) %>%
    tidyr::spread(key = sub2, value = n)

  # total mutations by sample to calculate percentages later
  totals <- table %>%
    dplyr::mutate(total = rowSums(.[-c(1:2)])) %>%
    dplyr::select(lymph, sample, total)

  table <- table %>%
    dplyr::mutate(`C>T` = `C>T at CpG` + `C>T other`) %>%
    tidyr::gather(key = mutation_type, value = n, -c(lymph, sample)) %>%
    dplyr::left_join(totals, by = c("lymph", "sample")) %>%
    dplyr::mutate(perc = n / total * 100) %>%
    dplyr::select(-total) %>%
    dplyr::mutate(mutation_type = factor(mutation_type, levels = mutation_order)) %>%
    dplyr::arrange(lymph, sample, mutation_type)

  return(table)
}


#' Makes mutation matrix with motiffs rows and sample in columns
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


#' Find the rows were the tables are located in a file with multiple tables
#' How? Taking groups of meaningful rows (not empty or commentary)
find_tables <- function(filepath) {
  con <- file(filepath, "r")
  lines <- readLines(con)
  close(con)

  l <- list()
  lidx <- 1
  in_table <- FALSE
  for (i in seq_len(length(lines))) {
    if ((nchar(lines[i]) == 0 || startsWith(lines[i], "#"))) {
      if (in_table) {
        l[[lidx]] <- c(t_begin, i - 1)
        lidx <- lidx + 1
        in_table <- FALSE
      }
    }
    else if (!in_table) {
      t_begin <- i
      in_table <- TRUE
    }
  }

  return(l)
}


#' Load HsMetrics from different files and return them in a dataframe.
#' Only certain columns are selected
load_hs_metrics <- function(directory) {
  filenames <- list.files(directory, full.names = TRUE)

  df <- do.call(rbind, lapply(filenames, function(f) {
    tp <- find_tables(f)

    t <- read.table(f, sep = "\t", header = TRUE, skip = tp[[1]][1] - 1, nrows = tp[[1]][2] - tp[[1]][1])
    t <- t %>%
      dplyr::select(
        TOTAL_READS,
        PCT_PF_UQ_READS_ALIGNED,
        PCT_SELECTED_BASES,
        MEAN_BAIT_COVERAGE,
        PCT_TARGET_BASES_10X,
        TARGET_TERRITORY,
        ON_TARGET_BASES
        ) %>%
      dplyr::mutate(Case = str_extract(f, "[:digit:]+_(CLL|NO|FL|FB|Fibro|PBL|BM|GL|MBL)([:digit:]?)")) %>%
      dplyr::mutate(type = gsub("[0-9]+_", "", Case)) %>%
      dplyr::mutate(group = ifelse(type %in% c("NO", "FB", "Fibro", "GL"), "normal", "tumor"))

    return(t)
  }))

  return(df)
}


#' Summarizes dataframe of HsMetrics by a certain columns, calculating mean, sd and ci for mean
#'
#' @param data data
#' @param val_cols Number of columns with the values to summarize. It's assumed that they are the first columns.
summarize_hs_metrics <- function(data, val_cols = 5, by = group) {
  by <- dplyr::enquo(by)

  result <- data %>%
    dplyr::select(1:val_cols, !!by) %>%
    tidyr::gather(var, val, -!!by) %>%
    dplyr::group_by(!!by, var) %>%
    dplyr::summarise(
      mean = mean(val),
      sd = sd(val),
      n = n()
    ) %>%
    dplyr::mutate(
      low_ci = mean - ci(sd, n),
      high_ci = mean + ci(sd, n)
    ) %>%
    as.data.frame()

  return(result)
}


load_cosmic_signatures <- function(path = NA) {
  if (is.na(path)) {
    path <- "https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  }
  signatures <- read.table(path, sep = "\t", header = TRUE)
  signatures <- signatures[, 1:33] # discards some empty columns in the file
  signatures <- signatures %>%
    dplyr::rename(substitution = Substitution.Type, context = Trinucleotide, mutation = Somatic.Mutation.Type) %>%
    dplyr::mutate(context = sub("(.)(.)(.)", "\\1.\\3", context)) %>% # context
    dplyr::arrange(substitution, context) # standard order

  return(signatures)
}


get_cross_similarity <- function(sigs_a, sigs_b) {
  dists <- sapply(seq_len(ncol(sigs_a)), function(x) {
    sapply(seq_len(ncol(sigs_b)), function(y) {
      1 - cosine_dist(sigs_a[, x], sigs_b[, y])
    })
  })

  colnames(dists) <- colnames(sigs_a)
  rownames(dists) <- colnames(sigs_b)

  return(dists)
}


# Get the names of the closest cosmic signatures to our signature
get_top_reference_signatures <- function(reference_signatures, similarities, signature_name, n = 5) {
  top_names <- rownames(similarities)[order(similarities[, signature_name], decreasing = TRUE)[1:n]]

  top <- reference_signatures[, c("substitution", "context", top_names)]

  return(top)
}


count_and_get_perc <- function(data, count_vars, group_vars) {
  counts <- data %>%
    dplyr::count(!!!count_vars) %>%
    dplyr::group_by(!!!group_vars) %>%
    dplyr::mutate(perc = n * 100 / sum(n)) %>%
    as.data.frame()

  return(counts)
}


count_aid_motifs <- function(data) {
  motif_counts <- count_and_get_perc(
    data,
    quos(lymph, aid_motif),
    quos(lymph)
  ) %>%
    dplyr::mutate(aid_motif = factor(aid_motif,
      levels = c("WRCY", "WA", "RCG", "None")
    ))

  return(motif_counts)
}


count_location <- function(data) {
  return(count_and_get_perc(
    data,
    quos(location, aid_motif),
    quos(location)
  ))
}


count_location_pathway <- function(data) {
  return(count_and_get_perc(
    data,
    quos(pathway, location, aid_motif),
    quos(location)
  ))
}


count_genes <- function(data) {
  return(count_and_get_perc(
    data,
    quos(lymph, pathway, geneName),
    quos(lymph, pathway)
  ))
}


count_sample <- function(data) {
  return(count_and_get_perc(
    data,
    quos(lymph, sample),
    quos(lymph)
  ))
}


count_sample_pathway <- function(data) {
  return(count_and_get_perc(
    data,
    quos(pathway, lymph, sample),
    quos(lymph)
  ))
}


count_lymph <- function(data) {
  return(count_and_get_perc(
    data,
    quos(lymph),
    quos(lymph)
  ))
}


count_pathway <- function(data) {
  return(count_and_get_perc(
    data,
    quos(pathway, lymph),
    quos(lymph)
  ))
}


get_sca_data <- function(dataset_name) {
  gr_data <- SomaticCancerAlterations::scaLoadDatasets(dataset_name, merge = TRUE)

  # Create dataframe
  df <- data.frame(
    sample = gr_data$Patient_ID,
    lymph = gr_data$Dataset,
    type = gr_data$Variant_Type,
    chrom = paste0("chr", seqnames(gr_data)),
    pos = start(ranges(gr_data)),
    ref = gr_data$Reference_Allele,
    alt = gr_data$Tumor_Seq_Allele2
  )

  # Filter autosomes chr and SNP
  df <- df %>%
    filter(chrom %in% paste0("chr", 1:22)) %>%
    filter(type %in% c("SNP"))

  return(df)
}


get_signature_contrib_bylymph <- function(contribution) {
  df <- as.data.frame(contribution) %>%
    tibble::rownames_to_column(var = "names") %>%
    dplyr::mutate(group = ifelse(grepl("FL_.*", names), "FL", "CLL/MBL")) %>%
    dplyr::select(-names) %>%
    dplyr::group_by(group) %>%
    tidyr::gather(signature, value, -group) %>%
    dplyr::group_by(group, signature) %>%
    dplyr::summarise(value = sum(value)) %>%
    dplyr::mutate(value = value * 100 / sum(value)) %>%
    tidyr::spread(signature, value)

  return(df)
}


get_mutation_repetitions <- function(data) {
  df <- data %>%
    dplyr::count(lymph, chrom, pos, ref, alt) %>%
    dplyr::count(lymph, n) %>%
    tidyr::complete(lymph, n) %>%
    replace(is.na(.), 0) %>%
    dplyr::group_by(lymph) %>%
    dplyr::mutate(cs = cumsum(nn), cs_n = cumsum(nn) / sum(nn) * 100)

  return(df)
}

load_mut_mat <- function(MUTATIONAL_MATRIX) {
  mut_mat <- read.csv(MUTATIONAL_MATRIX) %>%
    arrange(X) %>%
    mutate(X = gsub("(.)(.) (.)\\.(.)", "\\3[\\1>\\2]\\4", X)) %>%
    # normalize
    mutate_at(colnames(.)[-1], ~ (. / sum(.))) %>%
    column_to_rownames("X") %>%
    as.matrix()

  return(mut_mat)
}

load_mut_mat_mx <- function(data) {
  mut_mat <- data %>%
    as.data.frame() %>%
    rownames_to_column(var = "X") %>%
    mutate(X = gsub("(.)(.) (.)\\.(.)", "\\3[\\1>\\2]\\4", X)) %>%
    # normalize
    mutate_at(colnames(.)[-1], ~ (. / sum(.))) %>%
    column_to_rownames("X") %>%
    as.matrix()

  return(mut_mat)
}

tidying_for_ds <- function(data) {
  mut_mat_ds <- data %>%
    t() %>%
    as.data.frame()

  return(mut_mat_ds)
}

reference_sig <- function(REFERENCE_SIGNATURES, REFERENCE_SUBSET_DS) {
  ref_signatures <- read.csv(REFERENCE_SIGNATURES) %>%
    unite("X", Type, SubType, sep = " ") %>%
    arrange(X) %>%
    mutate(X = gsub("(.)>(.) (.).(.)", "\\3[\\1>\\2]\\4", X)) %>%
    column_to_rownames("X") %>%
    as.matrix()


  if (exists("REFERENCE_SUBSET_DS")) {
    ref_signatures_ds <- ref_signatures[, REFERENCE_SUBSET_DS] %>%
      t() %>%
      as.data.frame()
  } else {
    ref_signatures_ds <- ref_signatures %>%
      as.data.frame()
  }

  return(ref_signatures_ds)
}

get_lymph_group <- function(names) {
  return(
    ifelse(
      grepl("^CLL(.MBL)?", names),
      "CLL/MBL",
      str_extract(names, "^[^_]*")
    )
  )
}