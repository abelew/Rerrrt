#' Given the many tables provided by quantify, merge them into matrices of error
#' types.
#'
#' @param tables The big pile of error tables.
#' @param samples The set of samples in the data.
#' @param verbose Print some information while running.
#' @return Matrices with error types as rows, samples as columns.
#' @export
matrices_from_tables <- function(tables, samples, verbose=FALSE) {
  matrices <- list()
  indexes_per_sample <- c()
  reads_per_sample <- c()
  num_tables <- length(tables)
  for (t in 1:num_tables) {
    table_name <- tables[t]
    col_names <- c()
    num_samples <- length(samples)
    for (s in 1:num_samples) {
      sample_name <- names(samples)[s]
      col_names <- c(col_names, sample_name)
      table <- samples[[sample_name]][[table_name]]
      last_col <- ncol(table)
      column <- table[[last_col]]
      if (is.null(matrices[[table_name]])) {
        matrices[[table_name]] <- data.table::data.table("names"=table[[1]], "s"=column)
        colnames(matrices[[table_name]]) <- c("names", sample_name)
      } else {
        tmpdt <- data.table::data.table(keep.rownames=TRUE, "names"=table[[1]], "s"=column)
        matrices[[table_name]] <- merge(matrices[[table_name]], tmpdt, by="names", all=TRUE)
        colnames(matrices[[table_name]])[s + 1] <- sample_name
      }
      indexes_per_sample[sample_name] <- samples[[sample_name]][["filtered_index_count"]]
      reads_per_sample[sample_name] <- sum(samples[[sample_name]][["changed_table"]][["all_reads"]]) +
        sum(samples[[sample_name]][["ident_table"]][["all_reads"]])
    } ## End iterating over every sample name.
    if (nrow(matrices[[table_name]]) > 0) {
      na_idx <- is.na(matrices[[t]])
      matrices[[t]][na_idx] <- 0
      numeric_rownames <- stringr::str_detect(string=rownames(matrices[[t]]),
                                              pattern="[[:digit:]]")
      if (sum(numeric_rownames) == length(numeric_rownames)) {
        numeric_order <- suppressWarnings(order(as.numeric(rownames(matrices[[t]]))))
        matrices[[t]] <- matrices[[t]][numeric_order, ]
      }
    }
  }
  retlist <- list(
    "matrices" = matrices,
    "indexes_per_sample" = indexes_per_sample,
    "reads_per_sample" = reads_per_sample)
  return(retlist)
}

#' Invoke some normalization methods on the various matrices from
#' matrices_from_tables().
#'
#' @param matrices Tables provided by matrices_from_tables().
#' @param reads_per_sample Numeric vector defining how many reads survived
#'  filtering for each sample.
#' @param indexes_per_sample Numeric vector defining how many indexes survived
#'  filtering for each sample.
#' @return Normalized matrices
#' @export
normalize_matrices <- function(matrices, reads_per_sample, indexes_per_sample) {
  ## I still need to figure out how I want to normalize these numbers.
  mat_cpm <- matrices
  mat_cpmlength <- matrices
  mat_over_counts <- matrices
  mat_over_countslength <- matrices
  normalized_by_counts <- matrices
  matrices_by_counts <- matrices
  sequence_length <- nrow(matrices[["miss_indexes_by_position"]])
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

    ## str_detect checks the table name for whether it is reads vs. indexes vs. sequencer.
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
                                              lib.sizes=indexes_per_sample), silent=TRUE)
      used_columns <- colnames(mat_over_counts[[table_name]])
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / indexes_per_sample[used_columns]
    }
    mat_over_countslength[[table_name]] <- mat_over_counts[[table_name]] / sequence_length
    if (class(mat_cpm[[table_name]])[1] == "try-error") {
      message("Normalization failed for table: ", table_name, ", setting it to null.")
      mat_cpm[[table_name]] <- NULL
      mat_cpmlength[[table_name]] <- NULL
    } else {
      mat_cpmlength[[table_name]] <- mat_cpmlength[[table_name]] / sequence_length
    }
  }

  retlist <- list(
    "mat_cpm" = mat_cpm,
    "mat_cpmlength" = mat_cpmlength,
    "mat_over_counts" = mat_over_counts,
    "mat_over_countslength" = mat_over_countslength)
  return(retlist)
}


add_density_tables <- function(first, second) {
  first_df <- as.data.frame(first)
  colnames(first_df) <- c("num", "first")
  second_df <- as.data.frame(second)
  colnames(second_df) <- c("num", "second")
  both_df <- merge(first_df, second_df, by="num", all=TRUE)
  rownames(both_df) <- both_df[["num"]]
  both_df[["num"]] <- NULL
  sum <- rowSums(both_df, na.rm=TRUE)
  idx <- order(as.numeric(names(sum)), decreasing=FALSE)
  sum <- sum[idx]
  return(sum)
}

make_frequency_matrix <- function(mtrx) {
  new <- mtrx
  csums <- colSums(mtrx)
  for (col in colnames(mtrx)) {
    new[, col] <- new[, col] / csums[col]
  }
  multiplier <- 1 / max(new)
  new <- new * multiplier
  return(new)
}

#' Given the large table of index counts for one sample, bring them together into a small matrix.
#'
#' @param composite The matrix to build.
#' @param piece The piece to add.
#' @param sname Sample name.
make_index_table <- function(composite, piece, sname) {
  piece_table <- table(piece[["index_count"]])
  composite <- merge_index_table(composite, piece_table, sname)
  return(composite)
}

merge_index_table <- function(composite, piece_table, sname) {
  sname <- as.character(sname)
  piece_dt <- data.table("num" = names(piece_table), sname=as.numeric(piece_table))
  colnames(piece_dt) <- c("num", sname)
  if (ncol(composite) == 0) {
    composite <- piece_dt
  } else {
    composite <- merge(composite, piece_dt, by="num", all=TRUE)
  }
  return(composite)
}

#' Extract the numbers from the growing density matrix for this specific sample.
#'
#' Then add them and make a new matrix with these numbers.
#'
#' @param composite Table to build up.
#' @param first First piece to add to the composite.
#' @param second Second piece to add to the composite.
#' @param sname Sample name.
add_index_densities <- function(composite, first, second, sname) {
  sname <- as.character(sname)
  first_piece <- first[, c("num", sname), with=FALSE]
  colnames(first_piece) <- c("num", "first")
  second_piece <- second[, c("num", sname), with=FALSE]
  colnames(second_piece) <- c("num", "second")
  both <- merge(first_piece, second_piece, by="num", all=TRUE)
  num_names <- as.numeric(both[["num"]])
  both[["num"]] <- NULL
  sum_table <- rowSums(both, na.rm=TRUE)
  names(sum_table) <- num_names
  composite <- merge_index_table(composite, sum_table, sname)
  return(composite)
}
