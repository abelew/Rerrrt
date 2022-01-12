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
  tags_per_sample <- c()
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
      tags_per_sample[sample_name] <- samples[[sample_name]][["filtered_tag_count"]]
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
    "tags_per_sample" = tags_per_sample,
    "reads_per_sample" = reads_per_sample)
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

#' Given the large table of tag counts for one sample, bring them together into a small matrix.
#'
#' @param composite The matrix to build.
#' @param piece The piece to add.
#' @param sname Sample name.
make_tag_table <- function(composite, piece, sname) {
  piece_table <- table(piece[["tag_count"]])
  composite <- merge_tag_table(composite, piece_table, sname)
  return(composite)
}

merge_tag_table <- function(composite, piece_table, sname) {
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
add_tag_densities <- function(composite, first, second, sname) {
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
  composite <- merge_tag_table(composite, sum_table, sname)
  return(composite)
}
