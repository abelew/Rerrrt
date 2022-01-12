#' Invoke some normalization methods on the various matrices from
#' matrices_from_tables().
#'
#' @param matrices Tables provided by matrices_from_tables().
#' @param reads_per_sample Numeric vector defining how many reads survived
#'  filtering for each sample.
#' @param tags_per_sample Numeric vector defining how many tags survived
#'  filtering for each sample.
#' @return Normalized matrices
#' @export
normalize_matrices <- function(matrices, reads_per_sample, tags_per_sample) {
  ## I still need to figure out how I want to normalize these numbers.
  mat_cpm <- matrices
  mat_cpmlength <- matrices
  mat_over_counts <- matrices
  mat_over_countslength <- matrices
  normalized_by_counts <- matrices
  matrices_by_counts <- matrices
  sequence_length <- nrow(matrices[["miss_tags_by_position"]])
  mnames <- names(matrices)

  for (t in 1:length(mnames)) {
    table_name <- mnames[t]
    if (nrow(matrices[[table_name]]) == 0) {
      mat_cpm[[table_name]] <- NULL
      mat_cpmlength[[table_name]] <- NULL
      mat_over_counts[[table_name]] <- NULL
      mat_over_countslength[[table_name]] <- NULL
      normalized_by_counts[[table_name]] <- NULL
      matrices_by_counts[[table_name]] <- NULL
      next
    }

    ## str_detect checks the table name for whether it is reads vs. tags vs. sequencer.
    read_table <- stringr::str_detect(string=table_name, pattern="reads")

    ## Recast the data table as a matrix for doing cpm and such, also pull the names
    ## back out to their original purpose, rownames.
    data_columns <- 2:ncol(mat_cpm[[table_name]])
    mtrx <- as.matrix(mat_cpm[[table_name]][, 2:ncol(mat_cpm[[table_name]])])
    rownames(mtrx) <- mat_cpm[[table_name]][["names"]]

    ## For a few tables, the control sample is all 0s, so let us take that into account.
    usable_columns <- colSums(mtrx) > 0
    mtrx <- mtrx[, usable_columns]
    ## Now copy in place the new data
    mat_cpm[[table_name]] <- mtrx
    mat_cpmlength[[table_name]] <- mtrx
    mat_over_counts[[table_name]] <- mtrx
    mat_over_countslength[[table_name]] <- mtrx
    normalized_by_counts[[table_name]] <- mtrx
    matrices_by_counts[[table_name]] <- mtrx

    ## Now perform some normalizations
    if (isTRUE(read_table)) {
      mat_cpm[[table_name]] <- try(edgeR::cpm(y=mat_cpm[[table_name]],
                                              lib.sizes=reads_per_sample), silent=TRUE)
      used_columns <- colnames(mat_over_counts[[table_name]])
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / reads_per_sample[used_columns]
    } else {
      mat_cpm[[table_name]] <- try(edgeR::cpm(y=mat_cpm[[table_name]],
                                              lib.sizes=tags_per_sample), silent=TRUE)
      used_columns <- colnames(mat_over_counts[[table_name]])
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / tags_per_sample[used_columns]
    }
    mat_over_countslength[[table_name]] <- mat_over_counts[[table_name]] / sequence_length
    if (class(mat_cpm[[table_name]])[1] == "try-error") {
      message("Normalization failed for table: ", table_name, ", setting it to null.")
      mat_cpm[[table_name]] <- NULL
      mat_cpmlength[[table_name]] <- NULL
    } else {
      mat_cpmlength[[table_name]] <- mat_cpmlength[[table_name]] / sequence_length
    }

    ## Final step, earlier we removed the names from the tables,
    ## We are going to need them for plotting
    mat_cpm[[table_name]] <- mat_cpm[[table_name]] %>%
      as.data.frame() %>%
      tibble::add_column(names = rownames(mat_cpm[[table_name]]),
                         .before = 1)
    mat_cpmlength[[table_name]] <- mat_cpmlength[[table_name]] %>%
      as.data.frame() %>%
      tibble::add_column(names = rownames(mat_cpmlength[[table_name]]),
                         .before = 1)
    mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] %>%
      as.data.frame() %>%
      tibble::add_column(names = rownames(mat_over_counts[[table_name]]),
                         .before = 1)
    mat_over_countslength[[table_name]] <- mat_over_countslength[[table_name]] %>%
      as.data.frame() %>%
      tibble::add_column(names = rownames(mat_over_countslength[[table_name]]),
                         .before = 1)    
  } ## End iterating over the tables

  ## Since we removed the names earlier, put them back.
 
  retlist <- list(
    "mat_cpm" = mat_cpm,
    "mat_cpmlength" = mat_cpmlength,
    "mat_over_counts" = mat_over_counts,
    "mat_over_countslength" = mat_over_countslength)
  return(retlist)
}
