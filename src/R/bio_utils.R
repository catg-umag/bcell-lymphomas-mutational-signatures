library(VariantAnnotation)

# AID patterns; f/r = forward/reverse, m/c = mutation/context
aidp <- list(
  # Canonical AID signature should be C>T/G RCY
  WRCY = c(f = c(m = "C-[TGA]", c = "[AT][AG].[CT]."), r = c(m = "G-[ACT]", c = ".[AG].[CT][AT]")),
  # non-canonical according to Kasar A>C at WA
  WA = c(f = c(m = "A-[TGC]", c = ".[AT].[ACGT]."), r = c(m = "T-[CGA]", c = ".[ACGT].[AT].")),
  # signature 9 accordign to Alexandrov C>T at N.G
  RCG = c(f = c(m = "C-[TGA]", c = ".[AG].G."), r = c(m = "G-[ACT]", c = ".C.[CT]."))
)

# AID patterns for 1-nucl context by side
aidp3 <- list(
  # Canonical AID signature should be C>T/G RCY
  RCY = c(f = c(m = "C-[T]", c = "[AG].[CT]"), r = c(m = "G-[A]", c = "[AG].[CT]")),
  # non-canonical according to Kasar A>C at WA
  WA = c(f = c(m = "A-[TGC]", c = "[AT].[ACGT]"), r = c(m = "T-[CGA]", c = "[ACGT].[AT]")),
  # signature 9 accordign to Alexandrov C>T at N.G
  RCG = c(f = c(m = "C-[T]", c = "[AG].G"), r = c(m = "G-[A]", c = "C.[CT]"))
)


hs_autosomes <- function() {
  chrs <- as.character(1:22)
  return(chrs)
}

hs_linear <- function() {
  chrs <- c(1:22, "X", "Y")
  return(chrs)
}

ncbi <- function(x) {
  suppressMessages(GenomeInfoDb::seqlevelsStyle(x) <- "NCBI")
  GenomeInfoDb::genome(x) <- NA # avoid mismatches in 'genome' slots for overlaps
  return(x)
}


#' Load all the vcf from the directory passed as argument in a single dataframe.
#'
#' @param vcf_dir Directory with VCFs
#' @param anno Extra fields to extract from VCF
load_vcfs_in_dataframe <- function(vcf_dir, anno = c()) {
  vcf_files <- list.files(vcf_dir, pattern = "*.vcf", full.names = TRUE)

  data <- do.call(rbind, lapply(vcf_files, function(f) {
    vcf <- ncbi(VariantAnnotation::readVcf(f, "GRCh38"))
    GenomeInfoDb::seqlevels(vcf, pruning.mode = "coarse") <- hs_autosomes()

    name <- gsub("[.].*$", "", gsub("^.*/", "", f))

    df <- data.frame(
      sample = name,
      lymph = ifelse(grepl("FL_.*", name), "FL", "CLL/MBL"),
      chrom = gsub("chr", "", seqnames(vcf)),
      pos = start(ranges(vcf)),
      ref = ref(vcf),
      alt = as.character(unlist(alt(vcf)))
    )

    if (length(anno) > 0) {
      # posible fields to extract from VCF (df name => vcf info name)
      vcf_fields <- c(
        geneName = "Gene.refGene",
        location = "Func.refGene",
        effect = "ExonicFunc.refGene",
        SIFT = "SIFT_pred",
        LRT = "LRT_pred",
        ExAC = "ExAC_ALL"
      )


      for (field in anno) {
        if (field == "tumor_vaf") {
          df[["tumor_vaf"]] <- as.numeric(gsub("%", "", VariantAnnotation::geno(vcf)$FREQ[, "TUMOR"]))
        } else if (field == "SPV") {
          df[["SPV"]] <- vcf@info$SPV
        } else if (field == "1000genomes") {
          df$`1000genomes` <- vcf@info$`1000g2015aug_all`
        } else if (field %in% names(vcf_fields)) {
          df[[field]] <- unlist(vcf@info[[vcf_fields[[field]]]], use.names = FALSE)
          # replaces '\x3b' by ';'
          df[[field]] <- gsub("\\\\x3b", ";", df[[field]])
        }
      }
    }

    return(df)
  }))

  return(data)
}


#' Add context to mutation dataframe and convert mutations to pyrimidines (the standard).
#'
#' @param data Dataframe with mutations, must have $chrom, $pos, $ref and $alt
#' @param genome Genome used to extract contexts
get_context_and_convert_to_ct <- function(data, genome) {
  # Get contexts from reference (this one takes a bit)
  data$context2 <- as.character(BSgenome::getSeq(genome, data$chrom, data$pos - 2, data$pos + 2))
  data$context2 <- sub("(.{2})(.)(.{2})", "\\1.\\3", data$context2)
  data$context <- substr(data$context2, 2, 4)

  # Create substitution field and convert everything to C/T>
  data$substitution <- paste0(data$ref, "-", data$alt)

  ag_mutation <- data$ref %in% c("A", "G")
  data[ag_mutation, "substitution"] <- Biostrings::complement(
    Biostrings::DNAStringSet(data[ag_mutation, "substitution"])
  )
  data[ag_mutation, "context"] <- Biostrings::reverseComplement(
    Biostrings::DNAStringSet(data[ag_mutation, "context"])
  )
  data[ag_mutation, "context2"] <- Biostrings::reverseComplement(
    Biostrings::DNAStringSet(data[ag_mutation, "context2"])
  )

  data$substitution <- gsub("-", ">", data$substitution)

  return(data)
}


#' Check the presence of an AID spot on a given mutation
#'
#' @param mutation Mutation in format <ref>-<alt>
#' @param context Context in format CC.CC (one or two nucletides by each side depending of the pattern)
#' @param aid_patterns Patterns to search (check earlier in this file to view the format)
identify_aid_patterns <- function(mutation, context, aid_patterns = aidp) {
  matches <- names(aid_patterns)[sapply(aid_patterns, function(p) {
    (grepl(p["f.m"], mutation) & grepl(p["f.c"], context)) |
      (grepl(p["r.m"], mutation) & grepl(p["r.c"], context))
  })]

  if (length(matches) == 1) {
    return(matches)
  }
  else if (length(matches) == 0) {
    return("None")
  }
  else {
    print("Problem!")
  }
}
